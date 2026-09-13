# DBSI Explorer

A coarse, first-principles demonstration of **Direct Beam Solar Irradiance
(DBSI)** approximation, companion to the paper *"Parametric Expansion of
Clear-Sky Beam Solar Irradiance into Trigonometric Time-Series."*

Given a location (any latitude/longitude/altitude), a date/time range, and
a set of orbital elements, this tool computes closed-form DBSI0 and DBSI2
approximations -- no ephemeris file and no TSI dataset required -- and
optionally compares them against the Bird & Hulstrom (1981) broadband
transmittance model. Results are shown as a plot and can be exported to
CSV.

This repository is a focused demonstration tool, not the paper's full
research codebase (validation studies, sensitivity sweeps, and the paper
source itself live in a separate, larger repository). It's meant to let
anyone explore the approximation's behavior directly, with a small GUI and
no data downloads.

**Start here:** [`DBSI_GUI/README.md`](DBSI_GUI/README.md) for
requirements, usage, inputs, config file format, and outputs.

## Quick start

```bash
pip install -r requirements.txt
cd DBSI_GUI
python dbsi_gui_app.py
```

(Tkinter ships with the standard Windows/macOS Python installers; on Linux
it's typically a separate OS package, e.g. `sudo apt install python3-tk`.)

## Repository layout

```
DBSI-Explorer/
├── DBSI_GUI/        the GUI application and its computational engine
├── bin/             the underlying physics modules (vendored from the
│                    paper's own validated codebase; DBSI_GUI/ imports
│                    these directly, unmodified)
├── requirements.txt
└── LICENSE
```

## License

MIT -- see [`LICENSE`](LICENSE).
