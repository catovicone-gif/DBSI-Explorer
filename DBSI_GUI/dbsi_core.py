"""
dbsi_core.py -- self-contained-parameter wrapper around the paper's own
validated DBSI0/DBSI2 (and, optionally, Bird & Hulstrom 1981) pipeline,
with EVERY orbital, location, time and atmosphere input exposed as a
plain function argument instead of being read off a module-level default
constant (as PythonCode/apprDBSIv00_driver.py's appr_dbsi_pipeline() does).

This is the computational engine behind dbsi_gui_app.py. It deliberately
does NOT re-derive the physics from scratch: every formula below calls
straight into PythonCode/bin/*.py, the same modules the paper's own
validation scripts import, so this tool can never numerically drift from
what's actually published -- it only changes WHICH values feed those
formulas, from hardcoded defaults to caller-supplied parameters.

Run standalone (`python dbsi_core.py`) to self-test: reproduces
apprDBSIv00_driver.py's own default-parameter output bit-for-bit (to
float64 precision) as a regression check that this wrapper hasn't
silently diverged from the pipeline it wraps.

SCOPE NOTE (read before treating this as a general "any planet" tool):
Pa (orbital period) is a genuine free parameter -- it only ever appears
as 360/Pa*(t-t_p) in the mean-anomaly phase, so setting it to another
planet's period is meaningful. The rotation-rate side is NOT free the
same way: local_mean_hour_angle_deg() (bin/dbsi_common.py) derives the
Greenwich hour angle directly from the Julian Date's fractional part,
which is BY DEFINITION one Earth mean solar day -- there is no separate
"Pd" the way there is a "Pa". Changing "samples per day" in this tool's
time grid only changes sampling resolution, not the rotation period
actually used by the physics. Treat any other-planet exploration here as
"Earth's rotation rate, another planet's revolution period and orbital
geometry" -- not a full other-planet model.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin'))

import numpy as np

from gdt2jd import gdt2jd
from trig_utils import sind, cosd
from dbsi_common import (
    Location, LOCATIONS, TimeGrid, build_time_grid, wrap_to_pm180,
    cosz_from_declination, dbsi_series_terms, atmos_surface_transmittance,
    approx_solar_longitude_deg, local_mean_hour_angle_deg,
    ECCENTRICITY, OMEGA_TILDE_DEG, OBLIQUITY_DEG, ANOMALISTIC_YEAR_DAYS,
    SOLAR_CONSTANT_WM2, ATM_THICKNESS_KM, ATM_SURFACE_DENSITY_KGM3,
    ATM_EXTINCTION_M2_PER_KG,
)
from ampl_phase import ampl_phase
from equation_of_time import eot_formula_minutes, DEFAULT_TP_JD
from bird_hulstrom_1981 import bird_direct_transmittance, isa_surface_pressure_mb


def run_dbsi(loc, time_grid, *,
             eccentricity=ECCENTRICITY,
             obliquity_deg=OBLIQUITY_DEG,
             omega_tilde_deg=OMEGA_TILDE_DEG,
             t_p_jd=DEFAULT_TP_JD,
             orbital_period_days=ANOMALISTIC_YEAR_DAYS,
             tsi_wm2=SOLAR_CONSTANT_WM2,
             atm_thickness_km=ATM_THICKNESS_KM,
             atm_density_kgm3=ATM_SURFACE_DENSITY_KGM3,
             sigma_ext_m2kg=ATM_EXTINCTION_M2_PER_KG,
             include_bird=False,
             bird_tau_a038=0.039,
             bird_tau_a05=0.030,
             bird_ozone_atmcm=0.35,
             bird_water_cm=1.4):
    """
    Computes DBSI0 (eq. 3.1) and DBSI2 (eq. 3.2) over time_grid for one
    Location, using this paper's closed-form, ephemeris-free, TSI-file-free
    approximation chain (approx solar longitude eq. 2.6, equation of time
    eq. C.8, true local hour angle eq. 2.5, declination/zenith geometry,
    the paper's own single-extinction-cross-section spherical-shell
    transmittance) -- exactly PythonCode/apprDBSIv00_driver.py's own
    appr_dbsi_pipeline(), but with every orbital/atmosphere/TSI constant
    it reads from bin/dbsi_common.py module defaults instead exposed here
    as an explicit keyword argument.

    If include_bird is True, ALSO computes DBSI0_bird/DBSI2_bird using the
    Bird & Hulstrom (1981) broadband direct-beam transmittance in place of
    the paper's own transmittance model (same geometry, same S/T-series --
    only Theta(z) changes), for direct side-by-side comparison, exactly as
    Table 3/Figure 5 of the paper compare the two.

    Returns a dict of arrays (all the same length as time_grid.JD) plus a
    few scalar echo fields; see the bottom of this function for the full
    key list.
    """
    JD = time_grid.JD
    e = float(eccentricity)
    eps = float(obliquity_deg)
    omg = float(omega_tilde_deg)
    Pa = float(orbital_period_days)

    # --- Approximate solar longitude, eq. (2.6) ---
    l_sun_approx = approx_solar_longitude_deg(JD, t_p_jd, Pa=Pa, e=e, omg_deg=omg)

    # --- Equation of time (eq. C.8, post-round-25 appendix numbering) and
    # true local hour angle (eq. 2.5) ---
    EOT_min = eot_formula_minutes(JD, t_p_jd=t_p_jd, Pa=Pa, e=e, eps_deg=eps, omg_deg=omg)
    EOT_deg = EOT_min / 4.0
    Tsun_mean_local_deg = local_mean_hour_angle_deg(JD, loc.lon_deg)
    Tsun_true_local_deg = wrap_to_pm180(Tsun_mean_local_deg + EOT_deg)

    # --- Approximate declination and zenith geometry ---
    sd = sind(eps) * sind(l_sun_approx)
    decl_approx_deg = np.degrees(np.arcsin(np.clip(sd, -1, 1)))
    cosz_approx = cosz_from_declination(loc.lat_deg, decl_approx_deg, Tsun_true_local_deg)
    cosz_approx = np.clip(cosz_approx, 0.0, None)  # night -> 0
    zz_approx_deg = np.degrees(np.arccos(np.clip(cosz_approx, -1, 1)))

    # --- This paper's own transmittance model (spherical-shell, single
    # effective extinction cross section), at the site's own altitude ---
    z_norm, t_s_surf = atmos_surface_transmittance(
        alt_km=loc.alt_m / 1000.0, Hatm_km=atm_thickness_km,
        rho0=atm_density_kgm3, sigma_ext=sigma_ext_m2kg)
    transmissivity = np.interp(zz_approx_deg, z_norm, t_s_surf)

    # --- S/T-series terms, and DBSI0 (eq 3.1) / DBSI2 (eq 3.2) ---
    amp_l, phase_l, amp_t, phase_t = ampl_phase(eps, e, omg)
    terms = dbsi_series_terms(amp_l, phase_l, amp_t, phase_t, l_sun_approx, Tsun_true_local_deg)

    I00 = float(tsi_wm2)
    S01 = terms['S0'] + terms['S1']
    S23 = terms['S2'] + terms['S3']
    T0 = terms['T0']
    T12 = terms['T1'] + terms['T2']

    DBSI0_toa = I00 * (sind(loc.lat_deg) * S01 + cosd(loc.lat_deg) * T0)
    DBSI2_toa = DBSI0_toa + I00 * (sind(loc.lat_deg) * S23 + cosd(loc.lat_deg) * T12)

    DBSI0 = np.where(cosz_approx <= 0, 0.0, transmissivity * DBSI0_toa)
    DBSI2 = np.where(cosz_approx <= 0, 0.0, transmissivity * DBSI2_toa)

    out = dict(JD=JD, GD=time_grid.GD, l_sun_approx_deg=l_sun_approx,
               decl_approx_deg=decl_approx_deg, EOT_min=EOT_min,
               Tsun_true_local_deg=Tsun_true_local_deg, cosz_approx=cosz_approx,
               zz_approx_deg=zz_approx_deg, transmissivity=transmissivity,
               DBSI0=DBSI0, DBSI2=DBSI2, include_bird=include_bird)

    # --- Optional Bird & Hulstrom (1981) comparison, same geometry ---
    if include_bird:
        P_mb = isa_surface_pressure_mb(loc.alt_m)
        transmissivity_bird = bird_direct_transmittance(
            zz_approx_deg, P_mb, Uo_cm=bird_ozone_atmcm, Uw_cm=bird_water_cm,
            tau_a038=bird_tau_a038, tau_a05=bird_tau_a05)
        DBSI0_bird = np.where(cosz_approx <= 0, 0.0, transmissivity_bird * DBSI0_toa)
        DBSI2_bird = np.where(cosz_approx <= 0, 0.0, transmissivity_bird * DBSI2_toa)
        out.update(transmissivity_bird=transmissivity_bird,
                   DBSI0_bird=DBSI0_bird, DBSI2_bird=DBSI2_bird)

    return out


def make_time_grid(start, end, samples_per_day=24):
    """Thin convenience wrapper over bin/dbsi_common.build_time_grid() using
    a plain "samples per day" resolution knob instead of the (Pd, nstep)
    pair the underlying function takes (Pd=24h fixed -- see this module's
    SCOPE NOTE docstring).

    build_time_grid()'s own step size is (1.0/Pd)/nstep days -- i.e. nstep
    subdivisions of ONE Pd-th of a day (Pd=24 means Pd-ths are hours, so
    nstep=4 there means 4 subdivisions per HOUR, 96/day, matching its own
    docstring example "nstep=4 -> 15-minute steps"). That is a per-HOUR
    step count, not a per-DAY one, so it must be divided by Pd=24 here to
    turn a "samples_per_day" request into the nstep build_time_grid()
    actually wants -- passing samples_per_day straight through as nstep
    would silently produce 24x too many points."""
    return build_time_grid(start, end, Pa=None, Pd=24.0, nstep=float(samples_per_day) / 24.0)


def datetime_to_jd(dt_str_or_val):
    """Converts an ISO-8601 string (or numpy datetime64 / python datetime)
    to a Julian Date float, via the same gdt2jd() every driver in this
    project uses -- so a perihelion date typed into the GUI as plain text
    becomes the t_p_jd run_dbsi() needs."""
    return float(gdt2jd(np.datetime64(dt_str_or_val)))


if __name__ == '__main__':
    # Self-test: reproduce apprDBSIv00_driver.py's own default-parameter
    # output for its own default location/time-grid, bit-for-bit, as a
    # regression check that this wrapper's re-parameterized pipeline has
    # not drifted from the one actually used to produce the paper's figures.
    #
    # apprDBSIv00_driver.py itself lives in the full research repository,
    # not in this standalone demo (which vendors only the bin/ modules
    # run_dbsi() actually needs) -- so this bit-for-bit cross-check simply
    # skips itself here rather than failing. It does NOT affect run_dbsi()
    # or the GUI; it's a regression check for the paper's own dev repo.
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    try:
        from apprDBSIv00_driver import appr_dbsi_pipeline
    except ImportError:
        print("Skipping bit-for-bit self-test: apprDBSIv00_driver.py (from the paper's "
              "full research repository) is not part of this standalone demo repo.")
        print("run_dbsi() itself is unaffected -- see README.md to run the GUI instead.")
        sys.exit(0)

    from dbsi_common import SIDEREAL_DAY_HOURS

    loc = LOCATIONS['Alderville']
    tg = build_time_grid('2015-01-01', '2016-01-01',
                          Pa=ANOMALISTIC_YEAR_DAYS, Pd=SIDEREAL_DAY_HOURS, nstep=4)

    ref = appr_dbsi_pipeline(loc, tg, verbose=False)
    mine = run_dbsi(loc, tg)

    max_diff0 = np.max(np.abs(ref['DBSI0'] - mine['DBSI0']))
    max_diff2 = np.max(np.abs(ref['DBSI2'] - mine['DBSI2']))
    print(f"Self-test vs apprDBSIv00_driver.py, default params, Alderville, 2015:")
    print(f"  max|DBSI0 diff| = {max_diff0:.3e} W/m^2")
    print(f"  max|DBSI2 diff| = {max_diff2:.3e} W/m^2")
    assert max_diff0 < 1e-9 and max_diff2 < 1e-9, "dbsi_core.run_dbsi has DRIFTED from appr_dbsi_pipeline()!"
    print("  PASS -- bit-for-bit match.")

    # Bird & Hulstrom path, cross-checked against apprDBSIv00_driver's own
    # 'bird1981' transmittance_model option, same defaults.
    ref_bird = appr_dbsi_pipeline(loc, tg, verbose=False, transmittance_model='bird1981')
    mine_bird = run_dbsi(loc, tg, include_bird=True)
    max_diffb0 = np.max(np.abs(ref_bird['DBSI0'] - mine_bird['DBSI0_bird']))
    max_diffb2 = np.max(np.abs(ref_bird['DBSI2'] - mine_bird['DBSI2_bird']))
    print(f"\nSelf-test vs appr_dbsi_pipeline(transmittance_model='bird1981'):")
    print(f"  max|DBSI0 diff| = {max_diffb0:.3e} W/m^2")
    print(f"  max|DBSI2 diff| = {max_diffb2:.3e} W/m^2")
    assert max_diffb0 < 1e-9 and max_diffb2 < 1e-9, "Bird & Hulstrom path has DRIFTED!"
    print("  PASS -- bit-for-bit match.")
