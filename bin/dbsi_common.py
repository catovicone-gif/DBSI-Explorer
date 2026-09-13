"""
Common Functions shared by trueDBSIv00_driver.py and apprDBSIv00_driver.py,
per "High Level Description of Common Functions to both drivers" in
Instruction_for_DBSI_August25_2026.docx:

    (A) Location (name, longitude, time-zone, latitude, altitude)
    (B) Local Mean Time Array (JD array + matching Gregorian array, built
        from start/end date-time, Earth's revolution period, and Earth's
        rotation period)
    (C) Orbital elements (a, e, omega~) and other constants (obliquity,
        Sun's radius, Earth's radius, atmosphere thickness, atmospheric
        surface density, atmospheric extinction cross section) -- each
        below is tagged with its source or flagged as an assumption.

Both drivers import from here instead of duplicating this boilerplate, so
the upcoming benchmarking/validation scripts (next phase, per project plan)
get the same building blocks rather than a third copy of it.
"""
from dataclasses import dataclass
import numpy as np

from gdt2jd import gdt2jd
from jd2gdt import jd2gdt
from trig_utils import sind, cosd
from atmos_transmittance import atmos_transmittance


# =====================================================================
# (A) Location
# =====================================================================

@dataclass
class Location:
    """
    name:     human-readable label
    lon_deg:  geographic longitude, degrees, EAST positive (matches the
              sign convention already used throughout MatlabCode/DBSI00.m,
              e.g. Toronto l_geo=-79.631)
    lat_deg:  geographic latitude, degrees, north positive
    tz_hours: civil time-zone offset from UTC, hours -- informational only;
              the DBSI pipeline itself runs entirely in mean solar time
              (via lon_deg, see local_mean_hour_angle_deg below), not civil
              clock/zone time, so tz_hours is not consumed by the physics.
              It's carried on Location per instruction (A) for labeling
              output and for any future civil-time-facing I/O.
    alt_m:    altitude above sea level, meters. Consumed by
              atmos_surface_transmittance() below to evaluate the
              atmospheric transmittance at this site's actual elevation
              rather than assuming a sea-level surface.
    """
    name: str
    lon_deg: float
    lat_deg: float
    tz_hours: float
    alt_m: float


# Predefined locations, values carried over unchanged from MatlabCode/DBSI00.m
LOCATIONS = {
    'Toronto':    Location('Toronto',     -79.631,   43.651070, -5.0,  173.0),
    'Alderville': Location('Alderville',  -78.081,   44.190939, -5.0,  259.0),
    'Golden':     Location('Golden',     -105.221,   39.742300, -7.0, 1730.0),
    # Izana Atmospheric Observatory, Tenerife, Canary Islands (BSRN station
    # code "iza") -- added round 16, 2026-09-02, as a third validate_dbsiv00_
    # driver.py site chosen specifically for its elevation (~2,373 m) to
    # exercise the altitude-aware spherical-shell transmittance model at a
    # much higher elevation than Golden's ~1,730 m. Lat/lon/elevation here
    # are read directly from the real downloaded BSRN/PANGAEA data file's
    # own "Coverage"/"Event(s)" metadata (LATITUDE: 28.309350, LONGITUDE:
    # -16.499260, ELEVATION: 2372.9 m) -- more precise than, and superseding,
    # an earlier web-search-derived estimate (28.30, -16.48, 2367.0) used
    # before the real file was in hand. Civil time zone is WET (UTC+0
    # standard / UTC+1 DST, i.e. Canary Islands time, NOT mainland Spain's
    # CET/CEST) -- informational only, per the Location docstring; not
    # consumed by the physics.
    'Izana':      Location('Izana',       -16.499260, 28.309350,  0.0, 2372.9),
}


# =====================================================================
# (B) Local Mean Time Array
# =====================================================================

@dataclass
class TimeGrid:
    JD: np.ndarray        # Julian Date array, mean (UTC) time scale
    GD: np.ndarray         # matching Gregorian datetime64[us] array
    t: np.ndarray            # days elapsed since JD[0] (t = JD - JD[0])
    step_days: float


def build_time_grid(start, end, Pa=None, Pd=24.0, nstep=4):
    """
    (B) Local Mean Time Array.

    start, end: start/end date-times of the study period (numpy
        datetime64, python datetime, or ISO-8601 string). Interpreted as
        UTC/mean time, no timezone handling -- matches gdt2jd()'s own
        convention.
    Pa: Period of Earth's revolution, days. Unused in the grid construction
        itself (the grid is built from Pd, the rotation period, which sets
        the timestep) -- accepted here only so callers can pass ANOMALISTIC_YEAR_DAYS
        for symmetry with the instructions' "Period of Earth's revolution"
        requirement; kept for signature completeness / documentation, not
        because the time array needs it.
    Pd: Period of Earth's rotation, HOURS (24.0 = one mean solar day).
    nstep: timesteps per rotation period Pd (nstep=4 -> 15-minute steps,
        matching every existing driver script in this project).

    Returns a TimeGrid: JD (Julian Date, mean/UTC time scale -- required by
    the instructions to be one of the two output arrays), GD (the matching
    Gregorian Date/Time array -- the other required output array, built via
    jd2gdt.py's jd2gdt(), the corrected/canonical JD->Gregorian function --
    see jd2gdt.py's own docstring for the fix), and t (days since JD[0],
    the convenience "simple time scale" used throughout the existing bin/
    functions such as se_orbit()).
    """
    jd_start = float(gdt2jd(np.datetime64(start)))
    jd_end = float(gdt2jd(np.datetime64(end)))
    step = (1.0 / Pd) / nstep
    JD = np.arange(jd_start, jd_end + step / 2.0, step)
    GD = jd2gdt(JD)
    t = JD - JD[0]
    return TimeGrid(JD=JD, GD=GD, t=t, step_days=step)


# =====================================================================
# (C) Orbital elements and other constants
# =====================================================================
# Each constant below carries its source; where the existing MATLAB/Python
# codebase does not document a source, that is stated plainly rather than
# invented -- per the instructions' "all with supported sources or stated
# as assumptions."

# --- Earth's mean heliocentric orbital elements ---
# Values as used throughout MatlabCode/InsolApprox.m, DBSI00.m and
# ArxivPaperR03InsolationExpansion.tex Sec. 2 / Table D.1: annual mean
# values for calendar year 2022 (per the paper's own Table 1/2 caption).
ECCENTRICITY = 0.0167141              # e
OMEGA_TILDE_DEG = 103.1609            # omega~, longitude of perihelion, deg
OBLIQUITY_DEG = 23.44                 # eps, ecliptic obliquity, deg
ANOMALISTIC_YEAR_DAYS = 365.2421748385884   # Pa, Earth's revolution period, days
SIDEREAL_DAY_HOURS = 24.0             # Pd, Earth's rotation period (mean solar day), hours

# Semi-major axis / astronomical unit -- RESOLVED (2026-08-27): earlier
# revisions of this module carried a separate SEMI_MAJOR_AXIS_M = 149600438471.4 m
# (from MatlabCode/InsolApprox.m), distinct from and ~2,568 km larger than
# the IAU-defined astronomical unit (149597870700 m). The user confirmed
# the value in GreenEclSol.mat (i.e. the DE431-sourced value, below) is
# the correct one -- SEMI_MAJOR_AXIS_M is removed; AU_KM is now the single
# value used everywhere in this project's v00 codebase for this quantity
# (it was already the only one actually consumed by trueDBSIv00_driver.py/
# apprDBSIv00_driver.py -- SEMI_MAJOR_AXIS_M had been defined but unused).
# NOTE: MatlabCode/InsolApprox.m and the legacy Python parity port
# insol_approx_driver.py (built specifically to reproduce InsolApprox.m's
# published a1TSI/a2TSI numbers exactly) still use the old, unsourced
# 149600438471.4 m value by design -- changing those is a separate,
# user-confirm-first decision (see DBSIv00_open_issues.md), not folded in
# here automatically, since it would change what those legacy scripts are
# for (exact parity with previously-validated historical results).

# Astronomical unit, IAU-defined exact value (used for the Sun-Earth
# distance factor with DE431-derived r -- CORRECTED 2026-09-02, was
# mislabeled "DE441-derived" here and in read_green_sol.py; see that
# module's docstring for the correction -- which is expressed in AU in
# GreenEclSol.mat/.npz before that file's own AU field converts it to km,
# see read_green_sol.py's read_ecl_sol()). Matches GreenEclSol's own
# stored AU value exactly.
AU_KM = 149597870.7

# --- Sun ---
RSUN_M = 695700e3                     # Sun's radius, m (MatlabCode/InsolApprox.m)
SOLAR_CONSTANT_WM2 = 1361.0           # nominal TSI at 1 AU, W/m^2.
# SOURCE: ArxivPaperR03InsolationExpansion.tex Sec. 3 ("Currently, the
# solar constant value is 1361 W/m^2 [2]"). Used as-is (not re-derived via
# Stefan-Boltzmann law from an assumed effective temperature, unlike the
# older InsolApprox.m/insol_approx_driver.py, which back into ~1361 via
# Lo=Sigma*TempSun^4, TempSun=5777.257 K tuned to match -- using the
# paper's literal cited constant directly is simpler and matches eq.
# (3.1)-(3.2)'s own text calling this factor "the solar constant").

# --- Earth / atmosphere ---
# Values as passed to AtmosTransmittance.m / atmos_transmittance() in every
# existing driver: AtmosTransmittance(6375, 100, 1.22, 0.00350, z_norm).
REARTH_KM = 6375.0                    # mean Earth radius, km. ASSUMPTION:
# a simple sphere -- no documented source in the existing code for this
# specific radius value (it's close to, but not exactly, the IUGG mean
# radius 6371.0 km); flagging as an unsourced round-number assumption.
ATM_THICKNESS_KM = 100.0              # Hatm, effective atmosphere thickness, km.
# ASSUMPTION: a nominal, unsourced effective height for the exponential
# atmosphere model (roughly the altitude of the Karman line); not tied to
# a specific standard-atmosphere reference in the existing code.
ATM_SURFACE_DENSITY_KGM3 = 1.22       # rho_0, sea-level air density, kg/m^3.
# Close to the ICAO/US Standard Atmosphere sea-level value (1.225 kg/m^3);
# treated as an assumption here since the existing code does not cite the
# standard-atmosphere table directly.
ATM_EXTINCTION_M2_PER_KG = 0.00250    # sigma_ext, broadband extinction
# cross section, m^2/kg. RETUNED 2026-09-13 (was 0.00295, itself tuned
# 2026-08-28 from an original 0.00350): least-squares best fit against a
# POOLED Golden, CO + Izana, Tenerife clear-sky objective (Golden:
# Global-Diffuse, diffuse-free target; Izana: DIR*cos(z), direct-beam
# target) at the current Hatm/rho_0 -- see
# investigate_extinction_izana_extension.py, which extends the original
# Golden-only investigate_extinction_sensitivity.py by adding Izana (a
# high-altitude, low-turbidity site) to the objective, so the single
# shared sigma_ext is no longer tuned against one mid-altitude site alone.
# Standalone optima differ considerably by site: Golden alone still wants
# 0.00295 (RMSE 22.3 W/m^2), while Izana alone wants 0.00200 (RMSE 11.7
# W/m^2) -- reflecting real, unmodeled site-to-site turbidity/aerosol
# differences that this single non-spectral parameter cannot capture
# simultaneously. The pooled Golden+Izana objective (equal per-sample
# weight) is minimized at 0.00250 (pooled RMSE 27.0 W/m^2, vs. 37.2 W/m^2
# if the old Golden-only 0.00295 value were kept); at 0.00250, Golden's own
# RMSE worsens to 36.0 W/m^2 (bias +22.0) while Izana's improves to 21.3
# W/m^2 (bias -14.9) -- a deliberate compromise, not an improvement for
# either site individually. Alderville, Ontario was evaluated but NOT
# included in this objective (see that script's docstring): its target
# (G1) includes diffuse irradiance, which is not fully separable from the
# extinction-driven direct-beam signal being tuned here; its standalone
# optimum (0.00200, informational only) is broadly consistent with
# Izana's. For the joint Hatm/sigma_ext check (Golden-only, predates the
# Izana extension), see investigate_hatm_sigext_joint_sweep.py (that joint
# sweep found a lower-RMSE combination near Hatm=300 km, but it was NOT
# adopted -- it pushes Hatm far past its own stated ~100 km physical
# rationale for only a marginal further RMSE gain, along a broad/shallow,
# near-degenerate ridge; see that script's docstring). NOT a universal
# physical constant -- like Hatm and rho_0 below, this is a single,
# non-spectral, TUNABLE stand-in for the whole spectrum and all
# atmospheric constituents at once (see the paper's own Sec. 3: "no
# spectral phenomena were simulated for this benchmark"); the paper
# explicitly defers a detailed transmittance derivation to a separate
# paper (Appendix F). A different site/day would likely warrant a
# different value.


# =====================================================================
# Shared atmospheric-transmittance helper (site-altitude-aware)
# =====================================================================

def atmos_surface_transmittance(alt_km=0.0, Rearth_km=REARTH_KM, Hatm_km=ATM_THICKNESS_KM,
                                 rho0=ATM_SURFACE_DENSITY_KGM3,
                                 sigma_ext=ATM_EXTINCTION_M2_PER_KG, z_norm=None):
    """
    Returns (z_norm, t_s_surf): spherical-shell atmospheric transmittance
    as a function of zenith distance z_norm (degrees, 0..90, with a padded
    91st entry at z=90 -- see below), evaluated at a site alt_km above sea
    level (default 0.0, sea level).

    atmos_transmittance()'s own altitude grid ha runs from ha=0 (top of
    atmosphere) to ha=Hatm_km (sea-level surface) in 0.5 km steps. A site
    alt_km above sea level sees a thinner atmosphere below it than a
    sea-level site does -- equivalent to evaluating the transmittance at
    ha = Hatm_km - alt_km instead of ha = Hatm_km. This replaces the ad
    hoc "t_ss(:,end-4)"-style column offsets used for site-elevation
    adjustment in MatlabCode/InsolApproxB.m (Golden, alt ~1.73 km,
    approximated there by a fixed 2 km/4-column shift on a coarser
    Hatm=70 km grid, alongside a separately re-tuned extinction cross
    section, 0.005 vs Alderville's 0.0035) with a direct, altitude-value
    lookup on the ONE common atmosphere model (Hatm_km, rho0, sigma_ext --
    the same for every location) -- i.e. this deliberately does NOT
    re-tune the atmospheric physics per site, only the geometry (how much
    atmosphere sits between the site and space).

    z=90 is padded with the same *0.9 convention used throughout this
    project (sec(90) = infinity, so z_norm itself only runs 0..89 inside
    atmos_transmittance(); this appends one more point for z=90).
    """
    if z_norm is None:
        z_norm = np.arange(0, 91)
    _t_pp, t_ss, ha, _za = atmos_transmittance(Rearth_km, Hatm_km, rho0, sigma_ext, z_norm)
    target_ha = Hatm_km - alt_km
    col = int(np.clip(round(target_ha / 0.5), 0, len(ha) - 1))
    t_s_surf = t_ss[:, col].copy()
    t_s_surf = np.append(t_s_surf, t_s_surf[-1] * 0.9)
    return z_norm, t_s_surf


# =====================================================================
# Shared 1 AU un-normalization / inverse-square-law helper
# =====================================================================

def au_distance_factor(r_km, au_km=AU_KM):
    """
    (AU/r)^2 -- the inverse-square-law factor between an irradiance value
    expressed at exactly 1 AU (as every solar irradiance CDR product this
    project has used is: NRLSSI2's TSI and NRLSSI3's SSI are both natively
    1 AU-normalized -- see NRLSSIv03/build_coarse_nrlssi3.py's docstring
    and the project status doc's "1 AU normalization" write-up) and the
    true irradiance actually received at Sun-Earth distance r_km on a
    given day.

    This is the SAME factor trueDBSIv00_driver.py's true_dbsi_pipeline()
    already applies to NormCoarseNRLSSI2.mat's 1 AU TSI (previously written
    inline as `(AU_KM / r_true_km) ** 2`) -- pulled out here, at the one
    place r(t) is already computed from GreenEclSol, so any other 1
    AU-referenced dataset (e.g. NRLSSI v3's SSI, once wired in) gets
    un-normalized the same way, by calling this function again at its own
    point of use, rather than needing a separate un-normalization pass
    over the raw file itself (see round 13 of the project status doc for
    the fuller reasoning).
    """
    return (au_km / np.asarray(r_km, dtype=float)) ** 2


# =====================================================================
# Shared approximate Earth/solar ecliptic longitude helper (eq. 2.6/A.6)
# =====================================================================

def approx_earth_longitude_deg(JD, t_p_jd, Pa=ANOMALISTIC_YEAR_DAYS,
                                e=ECCENTRICITY, omg_deg=OMEGA_TILDE_DEG):
    """
    Approximate Earth's heliocentric ecliptic longitude l_earth, from mean
    orbital elements only, eq. (2.6)/(A.6):

        M       = 360/Pa * (JD - t_p_jd)            (mean anomaly, deg)
        l_earth = M + omega~ + 2 e sin(M)            (deg, first order in e)

    Returns l_earth wrapped to [0, 360). JD, t_p_jd: Julian Date(s) and
    reference perihelion epoch (JD); Pa, e, omg_deg default to this
    project's mean orbital elements above but can be overridden.

    This is the single shared implementation of eq. (2.6)'s formula, used
    both for the Sun's apparent geocentric ecliptic longitude l_sun (=
    l_odot, see approx_solar_longitude_deg() below, feeding the S/T-series
    eq. 2.3/2.4) and, unchanged, for the obliquity term of the equation of
    time (eq. D.7, bin/equation_of_time.py) -- previously two independent,
    textually-duplicated copies of this same formula.
    """
    JD = np.asarray(JD, dtype=float)
    M = 360.0 / Pa * (JD - t_p_jd)
    l_earth = M + omg_deg + (180.0 / np.pi) * 2 * e * sind(M)
    return np.mod(l_earth, 360.0)


def approx_solar_longitude_deg(JD, t_p_jd, Pa=ANOMALISTIC_YEAR_DAYS,
                                e=ECCENTRICITY, omg_deg=OMEGA_TILDE_DEG):
    """
    Sun's apparent geocentric ecliptic longitude l_sun (= l_odot, eq. 2.6),
    the antipode of approx_earth_longitude_deg(): l_sun = l_earth - 180,
    wrapped to [0, 360).
    """
    l_earth = approx_earth_longitude_deg(JD, t_p_jd, Pa=Pa, e=e, omg_deg=omg_deg)
    return np.mod(l_earth - 180.0, 360.0)


# =====================================================================
# Shared hour-angle helper (longitude-aware)
# =====================================================================

def local_mean_hour_angle_deg(JD, lon_deg):
    """
    Mean (no equation-of-time correction) LOCAL solar hour angle, in
    degrees, wrapped to (-180, 180], for a site at longitude lon_deg
    (degrees, east positive).

    JD's fractional part is 0 at Greenwich mean noon and 0.5 at Greenwich
    mean midnight (standard astronomical convention -- see gdt2jd.py /
    jd2gdt.py docstrings), so (JD mod 1)*360 - 180... i.e. the existing
    project convention wrap_to_pm180((JD mod 1)*360) -- is the GREENWICH
    mean hour angle. Adding the site's longitude converts that to the
    local mean hour angle: this is eq. (2.5)'s "+ l" term (see paper Sec.
    2), which every driver before this common module implicitly set to
    zero (see hour_angle_eot_validation.py's own docstring, "this
    project's own Tsun construction has no explicit longitude term
    either").
    """
    JD = np.asarray(JD, dtype=float)
    Tsun_greenwich = np.mod(JD, 1.0) * 360.0
    Tsun_greenwich = np.where(Tsun_greenwich > 180, Tsun_greenwich - 360, Tsun_greenwich)
    Tsun_local = Tsun_greenwich + lon_deg
    return wrap_to_pm180(Tsun_local)


def wrap_to_pm180(deg):
    """Wrap an angle (degrees) to (-180, 180]."""
    d = np.mod(np.asarray(deg, dtype=float) + 180.0, 360.0) - 180.0
    return np.where(d == -180.0, 180.0, d)


def cosz_from_declination(lat_deg, decl_deg, hour_angle_deg):
    """
    cos(zenith distance) from latitude, solar declination and local (true
    or mean, as supplied) solar hour angle, all in degrees -- the standard
    spherical-astronomy relation used throughout this project.
    """
    return (sind(lat_deg) * sind(decl_deg)
            + cosd(lat_deg) * cosd(decl_deg) * cosd(hour_angle_deg))


# =====================================================================
# Individually-addressable S/T-series terms (for DBSI0/DBSI2 per paper
# eq. 3.1/3.2, and any higher-order combination benchmarking/validation
# scripts want later)
# =====================================================================

def dbsi_series_terms(amp_l, phase_l, amp_t, phase_t, Lsun, Tsun):
    """
    Returns the individual S0..S3 and T0..T4 terms (each an array, same
    shape as Lsun/Tsun) making up the S_l_odot (eq 2.3) and T_l_odot,t_odot
    (eq 2.4) series -- NOT summed, unlike bin/approx_dir_irr_r01.py, whose
    DI0 already sums all four S-terms together and whose DI1 already sums
    T0+T1 together. That's fine for a cumulative truncation-order study
    (what approx_dir_irr_r01.py was built for), but eq. (3.1)'s DBSI_0
    needs S0+S1 WITHOUT S2/S3, and T0 WITHOUT T1 -- so this function keeps
    every term separate and lets the caller group them however eq.
    (3.1)/(3.2) (or any other combination) requires.

    amp_l, phase_l, amp_t, phase_t: from ampl_phase(eps, e, omg) --
    amp_l=[L0,L1,L2,L3] (=S0..S3), phase_l=[l1,l2,l3], amp_t=[C0..C4]
    (=T0..T4), phase_t=[PHI1..PHI4].

    Term signs match bin/approx_dir_irr_r01.py's (validated) convention
    exactly (T1, T3 carry a minus sign; T0, T2, T4 a plus sign -- an
    alternating pattern, not arbitrary; see e.g. the analogous alternating
    sign in equation_of_time.py's E_obl series). NOTE: this alternating
    sign does not literally match how ArxivPaperR03InsolationExpansion.tex
    eq. (3.1)/(3.2) is typeset (all terms shown with a bare '+'), and the
    paper's own printed eq. (3.2) also has an inconsistency of its own --
    the T1 term's second phase is printed as 103.496 instead of matching
    the first phase's 103.718 (they should be equal per the general T1
    term's own symmetric form, cos(t-l+Phi1)+cos(t+l-Phi1)). Rather than
    transcribe those literal (and inconsistent) printed numbers, this
    function uses the general, already-validated formulas that exactly
    reproduce Table 1/Table 2's amplitudes -- worth reconciling with the
    paper text/erratum separately.

    Returns a dict: {'S0':.., 'S1':.., 'S2':.., 'S3':..,
                      'T0':.., 'T1':.., 'T2':.., 'T3':.., 'T4':..}
    """
    L0, L1, L2, L3 = amp_l
    l1, l2, l3 = phase_l
    C0, C1, C2, C3, C4 = amp_t
    PHI1, PHI2, PHI3, PHI4 = phase_t

    S0 = np.full_like(np.asarray(Lsun, dtype=float), L0)
    S1 = L1 * sind(Lsun - l1)
    S2 = L2 * sind(2 * Lsun - l2)
    S3 = L3 * sind(3 * Lsun - l3)

    T0 = C0 * cosd(Tsun)
    T1 = -C1 * (cosd(Tsun - Lsun + PHI1) + cosd(Tsun + Lsun - PHI1))
    T2 = C2 * (cosd(Tsun - 2 * Lsun + PHI2) + cosd(Tsun + 2 * Lsun - PHI2))
    T3 = -C3 * (cosd(Tsun - 3 * Lsun + PHI3) + cosd(Tsun + 3 * Lsun - PHI3))
    T4 = C4 * (cosd(Tsun - 4 * Lsun + PHI4) + cosd(Tsun + 4 * Lsun - PHI4))

    return {'S0': S0, 'S1': S1, 'S2': S2, 'S3': S3,
            'T0': T0, 'T1': T1, 'T2': T2, 'T3': T3, 'T4': T4}
