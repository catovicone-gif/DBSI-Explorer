"""
dbsi_config.py -- read/write the DBSI GUI's input parameters from/to a
plain-text INI-format configuration file (Python's stdlib configparser,
no extra dependency), so a full parameter set can be prepared in a text
editor or scripted, rather than re-typed into the GUI by hand every run.

See example_config.ini in this folder for a working example (Golden, CO,
this paper's own defaults). Run `python dbsi_config.py` to regenerate that
example file from DEFAULT_PARAMS below.

FILE FORMAT: standard INI, five sections -- [orbital], [location], [time],
[atmosphere], [bird]. Every key is optional; any key left out falls back
to the default noted in its comment in example_config.ini. Unknown keys
are ignored (not an error), so a config file only needs to state the
parameters the user actually wants to change.
"""
import configparser
import os

DEFAULT_PARAMS = {
    'orbital': {
        'eccentricity': 0.0167141,
        'obliquity_deg': 23.44,
        'omega_tilde_deg': 103.1609,
        # Matches bin/equation_of_time.py's own DEFAULT_TP_JD (2459583.437188867)
        # exactly -- this is the reference perihelion epoch every figure in the
        # paper is actually calibrated against. Since this is a MEAN-elements
        # model, any epoch works in principle (the mean anomaly is exactly
        # periodic), but using a DIFFERENT one here would silently shift the
        # phase relative to the paper's own validated results unless offset by
        # a whole number of orbital periods -- so this default is pinned, not
        # arbitrary. Changing it is fine; just know it moves the phase.
        'perihelion_datetime': '2022-01-03T22:29:33.118110',
        'orbital_period_days': 365.2421748385884,
    },
    'location': {
        'name': 'Golden',
        'latitude_deg': 39.7423,
        'longitude_deg': -105.221,
        'altitude_m': 1730.0,
    },
    'time': {
        'start_datetime': '2024-01-01T00:00:00',
        'end_datetime': '2025-01-01T00:00:00',
        'samples_per_day': 24,
    },
    'atmosphere': {
        'tsi_wm2': 1361.0,
        'sigma_ext_m2kg': 0.00250,
        'atm_thickness_km': 100.0,
        'atm_density_kgm3': 1.22,
    },
    'bird': {
        'include_bird': False,
        'tau_a038': 0.039,
        'tau_a05': 0.030,
        'ozone_atmcm': 0.35,
        'water_cm': 1.4,
    },
}

# (section, key) -> python type, used both to parse on load and to decide
# str()-formatting on save.
_TYPES = {
    ('orbital', 'eccentricity'): float,
    ('orbital', 'obliquity_deg'): float,
    ('orbital', 'omega_tilde_deg'): float,
    ('orbital', 'perihelion_datetime'): str,
    ('orbital', 'orbital_period_days'): float,
    ('location', 'name'): str,
    ('location', 'latitude_deg'): float,
    ('location', 'longitude_deg'): float,
    ('location', 'altitude_m'): float,
    ('time', 'start_datetime'): str,
    ('time', 'end_datetime'): str,
    ('time', 'samples_per_day'): int,
    ('atmosphere', 'tsi_wm2'): float,
    ('atmosphere', 'sigma_ext_m2kg'): float,
    ('atmosphere', 'atm_thickness_km'): float,
    ('atmosphere', 'atm_density_kgm3'): float,
    ('bird', 'include_bird'): bool,
    ('bird', 'tau_a038'): float,
    ('bird', 'tau_a05'): float,
    ('bird', 'ozone_atmcm'): float,
    ('bird', 'water_cm'): float,
}

_COMMENTS = {
    'orbital': "Earth's mean heliocentric orbital elements (paper Sec. 2 / Appendix E).",
    'location': "Observer location. 'name' is a label only -- lat/lon/alt below are what's used.",
    'time': "Study period and sampling resolution (samples_per_day: 24=hourly, 1440=1-minute).",
    'atmosphere': "Solar constant and this paper's own single-extinction-cross-section transmittance model.",
    'bird': "Optional Bird & Hulstrom (1981) comparison. include_bird=false skips this model entirely.",
}


def load_config(path):
    """Reads an INI config file and returns a nested dict matching
    DEFAULT_PARAMS's shape, with any key not present in the file filled in
    from DEFAULT_PARAMS. Raises FileNotFoundError if path doesn't exist,
    and ValueError (naming the offending key) if a value can't be parsed
    as its expected type."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Config file not found: {path}")

    cp = configparser.ConfigParser()
    cp.read(path)

    params = {sec: dict(vals) for sec, vals in DEFAULT_PARAMS.items()}
    for sec in cp.sections():
        if sec not in params:
            continue  # unknown section, ignore
        for key, raw in cp.items(sec):
            if (sec, key) not in _TYPES:
                continue  # unknown key, ignore
            typ = _TYPES[(sec, key)]
            try:
                if typ is bool:
                    params[sec][key] = cp.getboolean(sec, key)
                else:
                    params[sec][key] = typ(raw)
            except ValueError as exc:
                raise ValueError(f"Config file {path}: can't parse [{sec}] {key} = {raw!r} "
                                  f"as {typ.__name__}") from exc
    return params


def save_config(path, params):
    """Writes a nested dict (same shape as DEFAULT_PARAMS -- missing keys
    are filled from DEFAULT_PARAMS) out to an INI file at path, with a
    short explanatory comment above each section."""
    cp = configparser.ConfigParser()
    lines = ["# DBSI GUI configuration file",
             "# Companion tool to ArxivPaperR04InsolationExpansion.tex",
             "# Any key omitted falls back to this paper's own default value.",
             ""]
    for sec, defaults in DEFAULT_PARAMS.items():
        vals = {**defaults, **params.get(sec, {})}
        cp[sec] = {k: str(v) for k, v in vals.items()}
        lines.append(f"# {_COMMENTS.get(sec, '')}")
        lines.append(f"[{sec}]")
        for k, v in vals.items():
            lines.append(f"{k} = {v}")
        lines.append("")
    with open(path, 'w') as f:
        f.write("\n".join(lines) + "\n")


if __name__ == '__main__':
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'example_config.ini')
    save_config(out_path, DEFAULT_PARAMS)
    print(f"Wrote {out_path}")
    # Round-trip check
    reloaded = load_config(out_path)
    assert reloaded == {sec: {k: v for k, v in vals.items()} for sec, vals in DEFAULT_PARAMS.items()}, \
        "Round-trip load_config(save_config(...)) did not reproduce DEFAULT_PARAMS!"
    print("Round-trip load/save check: PASS")
