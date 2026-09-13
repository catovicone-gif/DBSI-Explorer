"""
bird_hulstrom_1981.py -- Python port of the "Bird Model" direct-beam
transmittance from:

    R.E. Bird and R.L. Hulstrom, "A Simplified Clear Sky Model for Direct
    and Diffuse Insolation on Horizontal Surfaces," SERI/TR-642-761,
    February 1981 (References/1981-BirdModelTR-642-761.pdf), Table 2-6
    ("Equations for Total Downward Irradiance for the Bird Model"), p.8.

Implemented 2026-08-28 at the user's request, following a discussion of
whether DBSIv00's current single-effective-extinction-cross-section
Theta(z) is adequate for this paper, or whether a more physically
structured (though still broadband, not spectral) transmittance model
would strengthen the validation section. This module is that alternative:
a decomposition into five distinct, independently-shaped attenuation
processes (Rayleigh scattering, aerosol extinction, water vapor
absorption, ozone absorption, uniformly-mixed-gas absorption) instead of
one lumped exponential -- the same structural idea as Bird & Riordan
(1984)/SPCTRAL2 (see DBSIv00_Bird1984_comparison.md), but broadband
(closed-form, no per-wavelength integration) rather than spectral, so it
is a genuinely small addition, not a sub-model requiring its own
extensive validation chapter (it inherits the 1981 report's own
validation against BRITE/SOLTRAN/Dave rigorous codes).

SCOPE: only the DIRECT-BEAM transmittance (Table 2-6's I_d equation,
divided by I_0*cosZ) is implemented here -- DBSI is a direct-beam-only
quantity, so the report's diffuse-irradiance terms (I_as, I_G, I_T) are
not needed and are not implemented.

INPUTS NOT AVAILABLE FROM THE VALIDATION DATASETS: neither
20150324_AFN12.csv (Alderville) nor N2024March09NLR.csv (Golden) records
surface pressure, precipitable water, ozone amount, or aerosol turbidity
-- the Bird model's own required inputs (Table 2-8, p.10) beyond zenith
angle. Standard/typical placeholder values are used here (documented at
each default below), site-altitude-adjusted only for pressure via the
ISA barometric formula. This mirrors, structurally, what DBSIv00's
current model already does (one shared, unmeasured-per-day parameter
set) -- the difference is that the Bird model's five parameters each
carry a specific, checkable physical meaning and a literature-typical
default, rather than one lumped, unpublished constant.
"""
import numpy as np


def isa_surface_pressure_mb(alt_m):
    """
    Standard-atmosphere (ISA) surface pressure at altitude alt_m, in
    millibars. Used as a default for Bird's P input when local measured
    pressure isn't available (as here) -- P0=1013.25 mb at sea level.
    """
    return 1013.25 * (1.0 - 2.25577e-5 * np.asarray(alt_m, dtype=float)) ** 5.25588


def bird_airmass(zenith_deg):
    """
    Kasten's relative airmass, M, as used in the Bird model (Table 2-6):
        M = [cosZ + 0.15 (93.885 - Z)^-1.25]^-1
    Z is clipped to [0, 90] deg -- the formula is not meant for Z>90
    (below horizon), and DBSI's own cos(z) factor zeroes those points out
    downstream anyway. Returns M (dimensionless), same shape as zenith_deg.
    """
    Z = np.clip(np.asarray(zenith_deg, dtype=float), 0.0, 90.0)
    cosZ = np.cos(np.radians(Z))
    return 1.0 / (cosZ + 0.15 * (93.885 - Z) ** (-1.25))


def bird_direct_transmittance(zenith_deg, P_mb, Uo_cm=0.35, Uw_cm=1.4,
                               tau_a038=0.039, tau_a05=0.030, Ba=0.84, K1=0.1):
    """
    Broadband direct-beam transmittance product from Bird & Hulstrom
    (1981), Table 2-6 -- i.e. DNI_Bird(Z) = TSI * this, and
    DBSI_Bird(Z) = TSI * cos(Z) * this, playing the same role DBSIv00's
    atmos_surface_transmittance() Theta(z) plays now.

        Theta_Bird = 0.9662 * T_R * T_o * T_UM * T_w * T_A

    Parameters
    ----------
    zenith_deg : array_like
        Solar zenith angle Z, degrees.
    P_mb : array_like or float
        Surface pressure, millibars (site-altitude-adjusted -- see
        isa_surface_pressure_mb()).
    Uo_cm : float
        Total ozone amount, atm-cm. Default 0.35: a typical global-mean
        column ozone value (Bird & Riordan 1984's own Golden comparisons
        used 0.31-0.344 atm-cm -- see DBSIv00_Bird1984_comparison.md);
        not measured for either validation date.
    Uw_cm : float
        Precipitable water vapor, cm, in a vertical column. Default 1.4:
        a typical mid-latitude value (Bird & Riordan 1984's own Golden
        August 1981 comparisons measured 1.35-1.42 cm); not measured for
        either validation date.
    tau_a038, tau_a05 : float
        Aerosol optical depth (broadband turbidity) at 0.38 and 0.5 um.
        Defaults 0.039/0.030, TUNED 2026-08-28 (were 0.13/0.10, a generic
        "moderately clear" literature guess): least-squares best fit
        against the Golden, CO clear-sky validation day (Global-Diffuse,
        diffuse-free target) -- see investigate_bird_turbidity_sensitivity.py.
        A single scalar turbidity-scale factor was swept (both values
        scaled together, ratio held fixed) rather than optimized
        independently, since these two clear-sky days can't separately
        constrain the aerosol's spectral shape. Golden RMSE 41.9->6.8
        W/m^2 at this value -- notably BELOW the current DBSIv00
        placeholder's own best achievable RMSE even after tuning (22.2
        single-parameter / 15.9 joint, see
        investigate_hatm_sigext_joint_sweep.py). tau_a05=0.030 sits below
        (not within) Bird & Riordan (1984)'s own measured range at Golden
        (0.10-0.51 across their several comparison days), suggesting this
        particular 2024-03-09 day was unusually clean rather than that an
        implausible parameter was needed to fit it. NOT a universal
        constant -- like DBSIv00's own SIGext, this is the single largest
        source of day-to-day/site-to-site uncertainty in the Bird model
        per the 1981 report's own Sec 3.1 ("the most significant
        attenuator at all zenith angles is the aerosol"); a different
        site/day would likely warrant a different value. (Alderville's
        own apparent optimum, k=0 i.e. zero aerosol, was NOT used -- it's
        an artifact of G1 including diffuse irradiance neither model
        produces, not a genuine atmospheric finding; see the sweep
        script's docstring.)
    Ba : float
        Aerosol forward-scattering ratio. Default 0.84, the report's own
        suggested value "unless good information on the aerosol is
        available" (p.9) -- note Ba only enters the (unused here) diffuse
        term, so it has no effect on this function's output; kept as a
        parameter for interface completeness/future diffuse-term use.
    K1 : float
        Aerosol-absorptance constant. Default 0.1, the report's own
        suggested default (p.9) for when the rural/urban split isn't
        known; the rural-aerosol-fitted value from the report is 0.0933.

    Returns
    -------
    Theta_Bird : ndarray, same shape as zenith_deg. Zero where the sun is
        below the horizon is NOT applied here (that's cos(Z)'s job, as
        with DBSIv00's existing Theta(z)) -- values for Z>90 are simply
        the Z=90 value (airmass is clipped, see bird_airmass()).
    """
    Z = np.asarray(zenith_deg, dtype=float)
    M = bird_airmass(Z)
    Mp = M * np.asarray(P_mb, dtype=float) / 1013.0  # pressure-corrected airmass, M'

    # --- Rayleigh scattering ---
    T_R = np.exp(-0.0903 * Mp ** 0.84 * (1.0 + Mp - Mp ** 1.01))

    # --- Ozone absorption ---
    Xo = Uo_cm * M
    T_o = (1.0
           - 0.1611 * Xo * (1.0 + 139.48 * Xo) ** (-0.3035)
           - 0.002715 * Xo / (1.0 + 0.044 * Xo + 0.0003 * Xo ** 2))

    # --- Uniformly mixed gas absorption (CO2, O2) ---
    T_UM = np.exp(-0.0127 * Mp ** 0.26)

    # --- Water vapor absorption ---
    Xw = Uw_cm * M
    T_w = 1.0 - 2.4959 * Xw / ((1.0 + 79.034 * Xw) ** 0.6828 + 6.385 * Xw)

    # --- Aerosol extinction (scattering + absorption) ---
    tau_A = 0.2758 * tau_a038 + 0.35 * tau_a05
    T_A = np.exp(-(tau_A ** 0.873) * (1.0 + tau_A - tau_A ** 0.7088) * M ** 0.9108)

    Theta_Bird = 0.9662 * T_R * T_o * T_UM * T_w * T_A
    return Theta_Bird
