"""
Python port of MatlabCode/bin/AmplPhase.m

Computes the S-series (S0-S3) and T-series (T0-T4) amplitudes and phases
of the DBSI trigonometric expansion, from obliquity eps, eccentricity e,
and longitude of perihelion omg.

Validated: reproduces the paper's Table 1 and Table 2 amplitudes exactly.
"""
import numpy as np
from trig_utils import sind, cosd, atan2d


def ampl_phase(eps, e, omg):
    E0 = (1 + e ** 2 / 2) / ((1 - e ** 2) ** 2)
    E1 = 2 * e / ((1 - e ** 2) ** 2)
    E2 = e ** 2 / 2 / ((1 - e ** 2) ** 2)

    se2 = sind(eps) ** 2
    se4 = sind(eps) ** 4
    se6 = sind(eps) ** 6
    se8 = sind(eps) ** 8
    se10 = sind(eps) ** 10

    S1 = sind(eps)
    S0 = 1 - (se2 / 4 + 3 * se4 / 64 + 5 * se6 / 256 + 175 * se8 / 16384 + 441 * se10 / 65536)
    S2 = se2 / 4 + se4 / 16 + 15 * se6 / 512 + 140 * se8 / 8192 + 1470 * se10 / 131072
    S4 = se4 / 64 + 3 * se6 / 256 + 140 * se8 / 16384 + 420 * se10 / 65536

    L1 = S1 * np.sqrt(E0 ** 2 - E0 * E2 * cosd(2 * omg) + (E2 / 2) ** 2)
    L0 = -S1 * E1 / 2 * sind(omg)
    L2 = -S1 * E1 / 2
    L3 = S1 * E2 / 2
    phi_y = E2 / 2 * sind(2 * omg)
    phi_x = E0 - E2 / 2 * cosd(2 * omg)
    l1 = atan2d(phi_y, phi_x)
    l1 = l1 + 360 if phi_y < 0 else l1

    C0 = E0 * S0 + E2 * S2 / 2 * cosd(2 * omg)
    C1 = E1 / 2 * np.sqrt(S0 ** 2 + S2 ** 2 / 4 + S0 * S2 * cosd(2 * omg))
    phi_y = sind(omg) * (S0 - S2 / 2)
    phi_x = cosd(omg) * (S0 + S2 / 2)
    PHI1 = atan2d(phi_y, phi_x)
    PHI1 = PHI1 + 360 if phi_y < 0 else PHI1

    C2 = np.sqrt((E0 * S2 / 2) ** 2 + E0 * E2 * S2 * (S0 - S4 / 2) * cosd(2 * omg) / 2
                 + E2 ** 2 * (S0 ** 2 + S4 ** 2 / 4 - S0 * S4 * cosd(4 * omg)) / 4)
    phi_y = E2 / 2 * (S0 + S4 / 2) * sind(2 * omg)
    phi_x = E0 * S2 / 2 + E2 / 2 * (S0 - S4 / 2) * cosd(2 * omg)
    PHI2 = atan2d(phi_y, phi_x)
    PHI2 = PHI2 + 360 if phi_y < 0 else PHI2

    C3 = E1 / 4 * np.sqrt(S2 ** 2 + S4 ** 2 - 2 * S2 * S4 * cosd(2 * omg))
    phi_y = sind(omg) * (S2 + S4)
    phi_x = cosd(omg) * (-S2 + S4)
    PHI3 = atan2d(phi_y, phi_x)
    PHI3 = PHI3 + 360 if phi_y < 0 else PHI3

    C4 = 0.5 * np.sqrt((E0 * S4) ** 2 + (E2 * S2) ** 2 - E0 * E2 * S2 * S4 * cosd(2 * omg))
    phi_y = (E2 * S2 / 4) * sind(2 * omg)
    phi_x = -(E0 * S4 / 2) + (E2 * S2 / 4) * cosd(2 * omg)
    PHI4 = atan2d(phi_y, phi_x)
    PHI4 = PHI4 + 360 if phi_y < 0 else PHI4

    amp_l = np.array([L0, L1, L2, L3])
    phase_l = np.array([l1, omg, 2 * omg])
    amp_t = np.array([C0, C1, C2, C3, C4])
    phase_t = np.array([PHI1, PHI2, PHI3, PHI4])
    return amp_l, phase_l, amp_t, phase_t
