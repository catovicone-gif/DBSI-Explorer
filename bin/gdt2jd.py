"""
Python port of MatlabCode/bin/gdt2jd.m (Gregorian date/time -> Julian Date).

NOTE -- a real inconsistency found while building this: gdt2jd.m and
jd2gdt.m (already ported as jd2gdt.py) disagree by exactly 12 hours on
which JD corresponds to midnight of a given calendar date.

  gdt2jd.m:  jd0 = 2451544.5  <-> dn0 = datenum(2000,1,1)   (midnight)
  jd2gdt.m:  jd1 = 2451545.0  <-> dn1 = datenum(2000,1,1)   (midnight)

By the standard astronomical convention (used throughout, e.g., the
DE441 files themselves -- see DE441_HeliCenter.txt's header, where JD
2459580.500000000 is printed against "2022-Jan-01 00:00:00.0000", i.e.
JD *.5 = midnight, JD *.0 = noon), gdt2jd.m's mapping is the correct one
and jd2gdt.m is off by 12 hours. Cross-check: JD 2451544.5 (gdt2jd's
midnight 2000-01-01) to JD 2459580.5 (DE441's midnight 2022-01-01) is
exactly 8036.0 days -- a whole number, as it must be between two
midnights (22 years incl. 6 leap days: 22*365+6=8036). Using jd2gdt.m's
2451545.0 instead gives 8035.5 days, not a whole number -- confirming
jd2gdt.m's epoch is off by half a day.

This port therefore follows gdt2jd.m's (correct) convention.

FIXED (2026-08-27): jd2gdt.py originally carried the 12-hour offset
described above. Per the user's confirmation, it's now fixed and is the
CANONICAL name for JD->Gregorian conversion throughout this project
(matching the MATLAB jd2gdt.m/gdt2jd.m naming pair) -- jd2gdt.py simply
re-exports this module's jd_to_datetime() under that name, so there is
only one implementation of the epoch arithmetic to keep correct. New
code should `from jd2gdt import jd2gdt` rather than importing
jd_to_datetime directly from here; jd_to_datetime() is kept for any
code still importing it that name (identical behavior).
"""
import numpy as np

_JD_MIDNIGHT_2000 = 2451544.5
_EPOCH_2000 = np.datetime64('2000-01-01T00:00:00', 'us')


def gdt2jd(dt):
    """
    dt: a numpy datetime64, Python datetime, or array-like of these
    (interpreted as UTC/TDB-agnostic calendar date/time, matching the
    original MATLAB function which does no timezone handling either).
    Returns the Julian Date(s), scalar or ndarray matching dt's shape.
    """
    dt64 = np.asarray(dt, dtype='datetime64[us]')
    delta_us = (dt64 - _EPOCH_2000).astype('int64')
    jd = _JD_MIDNIGHT_2000 + delta_us / 86400e6
    return jd


def jd_to_datetime(jd):
    """
    Inverse of gdt2jd() -- Julian Date(s) -> numpy datetime64[us]. Also
    available (and now preferred for new code) as jd2gdt() in
    bin/jd2gdt.py, which re-exports this same function under the
    MATLAB-matching name -- see this module's docstring.
    """
    jd = np.asarray(jd, dtype=float)
    delta_us = np.round((jd - _JD_MIDNIGHT_2000) * 86400e6).astype('int64')
    return _EPOCH_2000 + delta_us.astype('timedelta64[us]')
