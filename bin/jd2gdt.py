"""
Python port of MatlabCode/bin/jd2gdt.m (Julian Date -> Gregorian date/time).

BUG FIX (2026-08-27): this port originally used jd2gdt.m's own reference
epoch, JD 2451545.0 = 2000-01-01 00:00 -- i.e. treating JD *.0 as
midnight. That is off by 12 hours from the standard astronomical
convention (JD *.5 = midnight, JD *.0 = noon), which gdt2jd.py's
gdt2jd()/jd_to_datetime() already followed correctly, and which DE441's
own files use (see gdt2jd.py's docstring for the original cross-check
that found this: JD 2451544.5 <-> midnight 2000-01-01 <-> JD 2459580.5
<-> midnight 2022-01-01 is exactly 8036.0 days, a whole number as it
must be between two midnights; JD 2451545.0 gives a non-whole 8035.5).

Per the user's confirmation (2026-08-27) -- "Good catch! ... start using
it, I am more familiar with that naming convention" -- jd2gdt() is now
the canonical JD->Gregorian function name used throughout this project,
matching the MATLAB jd2gdt.m/gdt2jd.m naming pair. Rather than duplicate
the epoch arithmetic a second time (and risk it drifting out of sync
again), this now simply re-exports gdt2jd.py's already-correct, already
independently-verified implementation under the jd2gdt name.
"""
from gdt2jd import jd_to_datetime as jd2gdt
