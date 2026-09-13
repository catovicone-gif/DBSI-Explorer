"""
Degree-based trigonometric helpers.

MATLAB provides sind/cosd/secd/atan2d as built-in functions; plain numpy
does not, so this module supplies degree-based equivalents used across
every ported function in this /bin folder. Not itself a port of any one
MatlabCode/bin/*.m file -- it stands in for MATLAB's built-ins.
"""
import numpy as np


def sind(x):
    return np.sin(np.deg2rad(x))


def cosd(x):
    return np.cos(np.deg2rad(x))


def secd(x):
    return 1.0 / np.cos(np.deg2rad(x))


def atan2d(y, x):
    return np.degrees(np.arctan2(y, x))
