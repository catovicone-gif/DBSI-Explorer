"""
Analytic equation of time E_odot, following eq. (D.8) / Appendix D of
ArxivPaperR03InsolationExpansion.tex ("the approximation formula for EOT
developed by Claude.ai" per the AI Agents for physics research chat).

This replaces MatlabCode/EquationOfTime.m's approach -- which loads a
precomputed DE441 lookup table (GreenEquSol.mat, external to this
project) and linearly interpolates -- with a closed-form expression
using the same fixed MEAN orbital elements used throughout the rest of
this model (e, obliquity, longitude of perihelion), exactly as
InsolApprox.m/insol_approx_driver.py do for everything else. No external
ephemeris file is needed at call time.

    E_odot = E_ecc - E_obl                                      (D.1)
    E_ecc  = 2 e sin(M) + (5/4) e^2 sin(2M)                      (D.3)
    E_obl  = sum_{n=1..N} (-1)^(n+1) (y^n / n) sin(2 n l_earth)  (D.7)
             with y = tan^2(eps/2)

M is the mean anomaly and l_earth the (first-order, eq. 2.6/A.6)
approximate ecliptic longitude, computed via dbsi_common.py's
approx_earth_longitude_deg() -- the single shared implementation of eq.
(2.6)'s formula, also used by apprDBSIv00_driver.py for l_sun_approx (via
dbsi_common.approx_solar_longitude_deg(), l_earth's antipode), so the two
call sites can no longer drift apart. Using l_earth (= l_sun - 180 deg,
mod 360) instead of l_sun here is exact, not an approximation: (D.7) only
ever appears as
sin(2n * l_sun) for even n = 2,4,6,8, and sin(2n*(l_earth+180)) =
sin(2n*l_earth) identically, since 2n*180 deg is always a whole multiple
of 360 deg when n is an integer. So the +/-180 deg Earth/Sun longitude
convention (see plot_de441_driver.py's docstring) cannot matter here.

The RHS of (D.1)/(D.3)/(D.7) comes out in RADIANS (the coefficients 2e,
(5/4)e^2, y, y^2/2, ... are dimensionless prefactors on a sine, so the
whole sum is a radian-valued angle) -- confirmed numerically against
Table D.1 of the paper: 2e = 0.033428 rad * (1440/2pi) min/rad = 7.66
min, matching the well-known ~7.7-minute eccentricity amplitude of the
real equation of time; y = 0.043037 rad -> 9.86 min, matching the
well-known ~9.9-minute obliquity amplitude.

SIGN CONVENTION -- a real finding, not a guess: the paper's own D.6
remarks flag that "the sign convention for E_odot varies... equation
(D.8) should be spot-checked... before being adopted", and asks for
exactly the check done here. Implementing (D.1)=(D.3)-(D.7) completely
literally and simplifying algebraically: since l_earth = nu + omega~
(true longitude) and E_ecc = nu - M, E_obl = l_earth - alpha (with
alpha = RA), (D.1) reduces exactly to E_odot_literal = alpha - (M +
omega~) = alpha - L_mean, i.e. the NEGATIVE of the standard "mean minus
apparent" equation-of-time convention (and of eq. 2.5's own stated
intent, "positive when the true Sun leads the mean Sun, e.g. early
November"). Verified two independent ways: (1) this algebraic
simplification, and (2) numerically -- the literal (D.1) formula gives
-14.2 min on 2022-02-11 and +16.4 min on 2022-11-03 relative to real
almanac values of -14.2 and +16.4 (i.e. exactly backwards), and matches
the DE441-derived truth in eot_de441_comparison.py only after negating.
This function therefore returns E_obl - E_ecc (D.7 minus D.3, the
negative of the literal (D.1)/(D.8) as printed) so its sign matches the
real world and the paper's own stated convention -- flagging this for
the paper text/eq. (2.5) or (D.1) to be corrected for self-consistency,
rather than silently leaving the discrepancy unreported.
"""
import numpy as np
from trig_utils import sind
from gdt2jd import gdt2jd
from dbsi_common import approx_earth_longitude_deg

# Mean orbital elements used throughout this project (paper Table D.1)
E_ECC = 0.0167141          # eccentricity
EPS_DEG = 23.44            # obliquity, degrees
OMG_DEG = 103.1609         # longitude of perihelion, omega~, degrees
PA_DAYS = 365.2421748385884  # anomalistic year, days

# Reference perihelion epoch (JD), the mean of DE441_HeliCenter.txt's Tp
# column over calendar year 2022 (see eot_de441_comparison.py) -- the
# same "mean elements" epoch used for the DE441 comparison below. Only
# sets the phase (where M=0); override t_p_jd for a different reference.
DEFAULT_TP_JD = 2459583.437188867


def eot_formula_minutes(t_jd, t_p_jd=DEFAULT_TP_JD, Pa=PA_DAYS, e=E_ECC,
                         eps_deg=EPS_DEG, omg_deg=OMG_DEG, n_obliquity_terms=4):
    """
    Core formula, eq. (D.8) with the sign corrected per this module's
    docstring (E_obl - E_ecc, not the literal D.1's E_ecc - E_obl).
    t_jd: Julian Date(s), scalar or array.
    Returns EOT in MINUTES, positive = apparent Sun ahead of mean Sun
    (standard convention; matches real almanac values and the DE441
    comparison in eot_de441_comparison.py).
    """
    t_jd = np.asarray(t_jd, dtype=float)
    M = 360.0 / Pa * (t_jd - t_p_jd)  # mean anomaly, degrees

    E_ecc = 2 * e * sind(M) + 1.25 * e ** 2 * sind(2 * M)  # radians, eq. (D.3)

    # eq. (2.6)/(A.6), via dbsi_common's shared implementation (also used
    # by apprDBSIv00_driver.py for l_sun_approx -- see that function's
    # docstring for why using l_earth directly here, rather than l_sun, is
    # exact for eq. (D.7)'s sin(2n*l_odot) terms)
    l_earth = approx_earth_longitude_deg(t_jd, t_p_jd, Pa=Pa, e=e, omg_deg=omg_deg)

    y = np.tan(np.radians(eps_deg / 2.0)) ** 2
    E_obl = 0.0
    sign = 1.0
    for n in range(1, n_obliquity_terms + 1):
        E_obl = E_obl + sign * (y ** n / n) * sind(2 * n * l_earth)  # eq. (D.7)
        sign = -sign

    E_odot_rad = E_obl - E_ecc  # (D.1) with corrected sign -- see docstring
    return np.degrees(E_odot_rad) * 4.0  # rad -> deg -> minutes (1 deg = 4 min)


def equation_of_time(dt, t_p_jd=DEFAULT_TP_JD, **kwargs):
    """
    Equation of time (minutes) for any Gregorian date/time input.

    dt: numpy datetime64, Python datetime, or array-like of these.
    Any extra kwargs (Pa, e, eps_deg, omg_deg, n_obliquity_terms) are
    passed through to eot_formula_minutes() to use different elements.
    """
    jd = gdt2jd(dt)
    return eot_formula_minutes(jd, t_p_jd=t_p_jd, **kwargs)
