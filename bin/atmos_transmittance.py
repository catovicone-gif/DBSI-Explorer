"""
Python port of MatlabCode/bin/AtmosTransmittance.m

Computes transmittance of an exponential atmosphere under two geometries:
(I) plane-parallel and (II) spherical-shell, as a function of altitude
and zenith distance.
"""
import numpy as np
from trig_utils import cosd, secd


def atmos_transmittance(Rearth, Hatm, rho_0, SIGext, zd):
    zd = np.asarray(zd, dtype=float)
    ha = np.arange(0, Hatm + 0.5, 0.5)  # 0 = top of atmosphere, Hatm = surface
    Hb = 0.0
    in_ = 4  # integration points between consecutive altitude steps

    zz = zd[zd != zd[-1]]  # drop z=90 (sec(90) -> inf)

    tau_pp_cols, tau_ss_cols = [], []
    for i, h in enumerate(ha):
        # Plane-parallel atmosphere approximation
        tau_pp_temp = rho_0 * Hatm * secd(zz) * (np.exp(-(Hatm - h) / (Hatm - Hb)) - np.exp(-1))
        tau_pp_cols.append(tau_pp_temp)

        # Spherical-shell atmosphere approximation
        intg = (i + 1) * in_
        dh = h / intg if intg > 0 else 0.0
        h_temp = np.cumsum(np.ones(intg) * dh) if intg > 0 else np.array([])
        rho_temp = rho_0 * np.exp(-(Hatm - h_temp) / (Hatm - Hb)) if intg > 0 else np.array([])
        del_temp = h_temp / Rearth if intg > 0 else np.array([])

        tau_ss_temp = np.zeros_like(zz)
        for j in range(intg):
            tau_ss_temp = tau_ss_temp + rho_temp[j] * (1 + del_temp[j]) / \
                np.sqrt(cosd(zz) ** 2 + 2 * del_temp[j] + del_temp[j] ** 2) * dh
        tau_ss_cols.append(tau_ss_temp)

    tau_pp = np.column_stack(tau_pp_cols)
    tau_ss = np.column_stack(tau_ss_cols)
    t_pp = np.exp(-SIGext * tau_pp)
    t_ss = np.exp(-SIGext * tau_ss)
    return t_pp, t_ss, ha, zz
