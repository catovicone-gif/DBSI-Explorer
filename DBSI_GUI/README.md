# DBSI0 / DBSI2 Explorer

A small desktop GUI companion to *"Parametric Expansion of Clear-Sky Beam
Solar Irradiance into Trigonometric Time-Series"* (`ArxivPaperR03InsolationExpansion.tex`),
for exploring the paper's closed-form Direct Beam Solar Irradiance
approximations, **DBSI0** (eq. 3.1) and **DBSI2** (eq. 3.2), for any
location, any date/time range, and any set of orbital elements -- no
ephemeris file, no TSI dataset, nothing to download. Everything this tool
computes comes directly from the same validated modules
(`bin/*.py`, vendored into this repo -- see the note below) the paper's own
figures and validation tables were built from, so it can never numerically
drift from what's published.

This repository is a focused, "coarse demonstration" companion tool, not
the full research codebase behind the paper (which includes the dozens of
validation/sensitivity scripts, figures, and datasets used to write it).
It vendors only the eight `bin/` modules `dbsi_core.py` actually calls into.

Optionally, it also computes the same two quantities using the Bird &
Hulstrom (1981) broadband transmittance model in place of the paper's own
single-extinction-cross-section model, for direct side-by-side comparison
(matching Table 3 / Figure 5 of the paper).

## Requirements

- Python 3.9+
- `numpy`, `matplotlib`
- Tkinter (`tkinter` is bundled with the standard Windows/macOS Python
  installers; on Linux it's usually a separate OS package, e.g.
  `sudo apt install python3-tk` on Debian/Ubuntu)

No other files from this repository are needed at runtime beyond the
`bin/` folder at the repo root, which this folder imports via a relative
path -- keep `DBSI_GUI/` and `bin/` as siblings (their layout as cloned)
or update the `sys.path.insert(...)` line at the top of `dbsi_core.py` if
you move either one.

## Running it

```
cd PythonCode/DBSI_GUI
python dbsi_gui_app.py
```

This opens a window with input panels (orbital elements, location, time
range, atmosphere, and an optional Bird & Hulstrom panel), a **Run**
button, and an embedded plot. Two more buttons export the current run:
**Export CSV...** (the full computed time series) and **Save Plot PNG...**.

## Inputs

| Group | Fields | Notes |
|---|---|---|
| Orbital elements | eccentricity *e*, obliquity *ε*, longitude of perihelion *ω̃*, perihelion date/time, orbital period *Pa* | Defaults are Earth's 2022 mean values, the same ones used throughout the paper. |
| Location | preset dropdown (Toronto / Alderville / Golden / Izana / Custom) or manual latitude / longitude / altitude | Picking a preset auto-fills lat/lon/alt, which stay editable afterward. |
| Time range | start, end (UTC, ISO-8601), resolution | Set start/end 24h apart for a daily profile, a year apart for an annual one, or anything else. Resolution is samples/day (4, 24, 96, 288, 1440); pick a finer one for a shorter range. |
| Atmosphere | TSI / solar constant, σ_ext, *H_a*, *ρ₀* | This paper's own spherical-shell, single-extinction-cross-section transmittance model (Appendix D). |
| Bird & Hulstrom (1981) | on/off, turbidity at 0.38 µm and 0.5 µm, ozone column, precipitable water | When on, DBSI0/DBSI2 are ALSO computed with Bird & Hulstrom's transmittance in place of the paper's own, and both are plotted/exported side by side. Defaults are the paper's own Golden-tuned values (Appendix F.2); see that section for why they're tuned against Golden only, not pooled with Izana. |

### A scope note on "any planet"

The orbital *period* (*Pa*) is a genuine free parameter -- it only ever
enters as `360/Pa * (t - t_p)` in the mean-anomaly phase, so setting it to
another planet's period is physically meaningful. The *rotation* rate is
not adjustable the same way: the underlying hour-angle formula derives the
Greenwich hour angle directly from the Julian Date's fractional part,
which is Earth's mean solar day by definition. Changing the time
resolution here only changes how finely a day is sampled, not the
rotation period itself. So treat any other-planet exploration as "Earth's
rotation rate, another planet's revolution period and orbital geometry" --
not a complete other-planet model.

## Config files

Rather than re-typing every field by hand, a full parameter set can be
loaded from (or saved to) a plain-text `.ini` file via **Load Config...**
/ **Save Config...** (or the File menu). See `example_config.ini` in this
folder for a working example (Golden, CO, this paper's own defaults) and
`dbsi_config.py`'s docstring for the exact key list. Any key left out of
a config file falls back to that default, so a config only needs to state
what it's actually changing, e.g.:

```ini
[location]
name = Izana
latitude_deg = 28.30935
longitude_deg = -16.49926
altitude_m = 2372.9

[time]
start_datetime = 2010-05-16T00:00:00
end_datetime = 2010-05-17T00:00:00
samples_per_day = 1440
```

## Output

- **Plot**: DBSI0 and DBSI2 over the selected range (top), and DBSI2-DBSI0
  -- the contribution of the higher-order S2/S3/T1/T2 harmonic terms --
  underneath. If Bird & Hulstrom is enabled, its DBSI0/DBSI2 curves are
  overlaid (dashed) on the top panel.
- **CSV**: one row per time step, with Julian Date, calendar date/time,
  solar longitude, declination, equation of time, local hour angle,
  cos(zenith), transmissivity, DBSI0, DBSI2, and (if enabled) the Bird &
  Hulstrom transmissivity/DBSI0/DBSI2 columns.

## Files in this folder

- `dbsi_core.py` -- the computational engine (`run_dbsi()`); no GUI
  dependency, callable from a plain script or notebook. Running it directly
  (`python dbsi_core.py`) prints a friendly notice and exits cleanly here,
  since its bit-for-bit self-test needs `apprDBSIv00_driver.py`, which lives
  in the paper's full research repository, not this standalone demo -- it
  does not affect `run_dbsi()` itself.
- `dbsi_config.py` -- config file load/save (`load_config()`,
  `save_config()`), plus `DEFAULT_PARAMS`.
- `dbsi_gui_app.py` -- the Tkinter application.
- `example_config.ini` -- a working example config (regenerate it with
  `python dbsi_config.py`).

## Known limitation

This tool was developed and its computational core (`dbsi_core.py`) was
verified numerically in an environment without a display, so the GUI
layout itself has not been visually confirmed to look right end-to-end.
The underlying physics is bit-for-bit tested against the paper's own
pipeline (see above); if a widget looks misplaced or a control behaves
oddly, that's a layout issue, not a physics one -- please report it (or
fix it directly; the layout code in `dbsi_gui_app.py` is a straightforward
`ttk` grid/pack layout with no exotic dependencies).
