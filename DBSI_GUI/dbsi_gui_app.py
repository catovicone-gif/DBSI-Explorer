"""
dbsi_gui_app.py -- small desktop GUI for exploring the paper's closed-form
DBSI0 (eq. 3.1) and DBSI2 (eq. 3.2) approximations, with every orbital
element, location, time range, atmosphere parameter, and (optionally) Bird
& Hulstrom (1981) turbidity input adjustable by the user -- no ephemeris
file, no TSI dataset, nothing to download: everything this tool computes
comes from dbsi_core.run_dbsi(), which calls straight into the same
PythonCode/bin/*.py modules the paper's own validation scripts use.

USAGE
    python dbsi_gui_app.py

Requires: numpy, matplotlib, and a Tk-enabled Python (tkinter is part of
the standard library on Windows/macOS installers; on Linux it may need a
separate OS package, e.g. `sudo apt install python3-tk`).

See README.md in this folder for a full walkthrough, and
example_config.ini for the config-file format ("Load Config..." button).
"""
import os
import sys
import time
import traceback
import webbrowser

import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dbsi_core
import dbsi_config
from dbsi_core import LOCATIONS, Location

MAX_SAMPLES = 2_000_000  # sanity cap so a fat-fingered "1-minute steps for 50 years" can't hang the GUI
MAX_HATM_KM = 300.0      # transmittance table cost grows as H_a^2: ~0.6 s at 100 km, ~1 minute at 1000 km
MAX_BIRD_ALT_M = 40_000.0  # isa_surface_pressure_mb() goes NaN above ~44 km

# The paper this tool accompanies, shipped next to this script (and bundled
# into the PyInstaller executable, which unpacks its data files to _MEIPASS).
PAPER_PDF = os.path.join(getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__))),
                         'ArxivPaperR04InsolationExpansion.pdf')


class DBSIApp:
    def __init__(self, root):
        self.root = root
        root.title("DBSI0 / DBSI2 Explorer -- companion to ArxivPaperR04InsolationExpansion")
        root.geometry("1280x900")

        self.results = None  # last run_dbsi() output dict, for CSV/PNG export

        self._build_menu()
        self._build_input_panel()
        self._build_status_bar()  # before the plot panel, so the plot (not the status bar) gives up space
        self._build_plot_panel()

        self._apply_params(dbsi_config.DEFAULT_PARAMS)

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------
    def _build_menu(self):
        bar = tk.Menu(self.root)
        filemenu = tk.Menu(bar, tearoff=0)
        filemenu.add_command(label="Load Config...", command=self.on_load_config)
        filemenu.add_command(label="Save Config...", command=self.on_save_config)
        filemenu.add_separator()
        filemenu.add_command(label="Export CSV...", command=self.on_export_csv)
        filemenu.add_command(label="Save Plot as PNG...", command=self.on_save_png)
        filemenu.add_separator()
        filemenu.add_command(label="Quit", command=self.root.destroy)
        bar.add_cascade(label="File", menu=filemenu)
        self.root.config(menu=bar)

    def _build_input_panel(self):
        outer = ttk.Frame(self.root, padding=8)
        outer.pack(side=tk.TOP, fill=tk.X)

        left = ttk.Frame(outer)
        left.grid(row=0, column=0, sticky="nw", padx=(0, 12))
        mid = ttk.Frame(outer)
        mid.grid(row=0, column=1, sticky="nw", padx=(0, 12))
        right = ttk.Frame(outer)
        right.grid(row=0, column=2, sticky="nw")

        # --- Orbital elements ---
        f = ttk.LabelFrame(left, text="Orbital elements", padding=8)
        f.pack(fill=tk.X, pady=4)
        self.v_ecc = tk.StringVar()
        self.v_obliquity = tk.StringVar()
        self.v_omega = tk.StringVar()
        self.v_tp = tk.StringVar()
        self.v_period = tk.StringVar()
        self._labeled_entry(f, "Eccentricity, e:", self.v_ecc, 0)
        self._labeled_entry(f, "Obliquity, ε (deg):", self.v_obliquity, 1)
        self._labeled_entry(f, "Longitude of perihelion, ω~ (deg):", self.v_omega, 2)
        self._labeled_entry(f, "Perihelion date/time (UTC, ISO):", self.v_tp, 3, width=28)
        self._labeled_entry(f, "Orbital period, Pa (days):", self.v_period, 4, width=28)

        # --- Location ---
        f = ttk.LabelFrame(left, text="Location", padding=8)
        f.pack(fill=tk.X, pady=4)
        ttk.Label(f, text="Preset:").grid(row=0, column=0, sticky="w")
        self.v_preset = tk.StringVar()
        preset_names = list(LOCATIONS.keys()) + ["Custom"]
        preset_menu = ttk.Combobox(f, textvariable=self.v_preset, values=preset_names,
                                    state="readonly", width=18)
        preset_menu.grid(row=0, column=1, sticky="w", padx=4, pady=2)
        preset_menu.bind("<<ComboboxSelected>>", self.on_preset_change)
        self.v_lat = tk.StringVar()
        self.v_lon = tk.StringVar()
        self.v_alt = tk.StringVar()
        self._labeled_entry(f, "Latitude (deg, N+):", self.v_lat, 1)
        self._labeled_entry(f, "Longitude (deg, E+):", self.v_lon, 2)
        self._labeled_entry(f, "Altitude (m):", self.v_alt, 3)

        # --- Time range ---
        f = ttk.LabelFrame(mid, text="Time range", padding=8)
        f.pack(fill=tk.X, pady=4)
        self.v_start = tk.StringVar()
        self.v_end = tk.StringVar()
        self.v_spd = tk.StringVar()
        self._labeled_entry(f, "Start (UTC, ISO):", self.v_start, 0, width=22)
        self._labeled_entry(f, "End (UTC, ISO):", self.v_end, 1, width=22)
        ttk.Label(f, text="Resolution:").grid(row=2, column=0, sticky="w")
        spd_menu = ttk.Combobox(f, textvariable=self.v_spd, state="readonly", width=18,
                                 values=["4  (6-hourly)", "24  (hourly)", "96  (15-min)",
                                         "288  (5-min)", "1440  (1-min)"])
        spd_menu.grid(row=2, column=1, sticky="w", padx=4, pady=2)
        ttk.Label(f, text="(24h range -> use 96 or 1440; 1-year range -> use 4 or 24)",
                  foreground="gray30", font=("TkDefaultFont", 8)).grid(
            row=3, column=0, columnspan=2, sticky="w")

        # --- Atmosphere ---
        f = ttk.LabelFrame(mid, text="Atmosphere (this paper's own model)", padding=8)
        f.pack(fill=tk.X, pady=4)
        self.v_tsi = tk.StringVar()
        self.v_sigext = tk.StringVar()
        self.v_hatm = tk.StringVar()
        self.v_rho0 = tk.StringVar()
        self._labeled_entry(f, "TSI / solar constant (W/m²):", self.v_tsi, 0)
        self._labeled_entry(f, "σ_ext, extinction cross section (m²/kg):", self.v_sigext, 1)
        self._labeled_entry(f, "H_a, atmosphere thickness (km):", self.v_hatm, 2)
        self._labeled_entry(f, "ρ₀, surface density (kg/m³):", self.v_rho0, 3)

        # --- Bird & Hulstrom ---
        f = ttk.LabelFrame(right, text="Bird & Hulstrom (1981) comparison", padding=8)
        f.pack(fill=tk.X, pady=4)
        self.v_bird_on = tk.BooleanVar()
        chk = ttk.Checkbutton(f, text="Also compute DBSI0/DBSI2 with Bird & Hulstrom transmittance",
                               variable=self.v_bird_on, command=self.on_bird_toggle)
        chk.grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 4))
        self.v_tau038 = tk.StringVar()
        self.v_tau05 = tk.StringVar()
        self.v_ozone = tk.StringVar()
        self.v_water = tk.StringVar()
        self._bird_rows = []
        self._bird_rows.append(self._labeled_entry(f, "Turbidity τ @ 0.38µm:", self.v_tau038, 1))
        self._bird_rows.append(self._labeled_entry(f, "Turbidity τ @ 0.5µm:", self.v_tau05, 2))
        self._bird_rows.append(self._labeled_entry(f, "Ozone column (atm-cm):", self.v_ozone, 3))
        self._bird_rows.append(self._labeled_entry(f, "Precipitable water (cm):", self.v_water, 4))

        # --- Run / export buttons ---
        f = ttk.Frame(right, padding=8)
        f.pack(fill=tk.X, pady=4)
        ttk.Button(f, text="Run", command=self.on_run).grid(row=0, column=0, padx=2)
        ttk.Button(f, text="Load Config...", command=self.on_load_config).grid(row=0, column=1, padx=2)
        ttk.Button(f, text="Save Config...", command=self.on_save_config).grid(row=0, column=2, padx=2)
        ttk.Button(f, text="Export CSV...", command=self.on_export_csv).grid(row=1, column=0, padx=2, pady=4)
        ttk.Button(f, text="Save Plot PNG...", command=self.on_save_png).grid(row=1, column=1, padx=2, pady=4)

        self.on_bird_toggle()  # set initial enabled/disabled state

    def _labeled_entry(self, parent, label, var, row, width=14):
        ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w")
        entry = ttk.Entry(parent, textvariable=var, width=width)
        entry.grid(row=row, column=1, sticky="w", padx=4, pady=2)
        return entry

    def _build_plot_panel(self):
        frame = ttk.Frame(self.root)
        frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.fig = Figure(figsize=(11, 5.5), dpi=100)
        self.ax0 = self.fig.add_subplot(211)
        self.ax1 = self.fig.add_subplot(212)
        self.fig.tight_layout(pad=3.0)
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        # Toolbar is packed before the canvas so that, when the window is
        # shorter than the figure's natural height, the canvas is what
        # shrinks instead of the toolbar being pushed off-screen.
        toolbar = NavigationToolbar2Tk(self.canvas, frame, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _build_status_bar(self):
        self.status = tk.StringVar(value="Ready.")
        bar = ttk.Frame(self.root)
        bar.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Button(bar, text="Exit", command=self.root.destroy).pack(side=tk.RIGHT, padx=(2, 4), pady=2)
        ttk.Button(bar, text="Open Paper (PDF)", command=self.on_open_paper).pack(side=tk.RIGHT, padx=2, pady=2)
        ttk.Label(bar, textvariable=self.status, relief=tk.SUNKEN, anchor="w", padding=4).pack(
            side=tk.LEFT, fill=tk.X, expand=True)

    # ------------------------------------------------------------------
    # Parameter <-> widget plumbing
    # ------------------------------------------------------------------
    def _apply_params(self, params):
        o, l, t, a, b = (params['orbital'], params['location'], params['time'],
                         params['atmosphere'], params['bird'])
        self.v_ecc.set(o['eccentricity'])
        self.v_obliquity.set(o['obliquity_deg'])
        self.v_omega.set(o['omega_tilde_deg'])
        self.v_tp.set(o['perihelion_datetime'])
        self.v_period.set(o['orbital_period_days'])

        name = l.get('name', 'Custom')
        self.v_preset.set(name if name in LOCATIONS else "Custom")
        self.v_lat.set(l['latitude_deg'])
        self.v_lon.set(l['longitude_deg'])
        self.v_alt.set(l['altitude_m'])

        self.v_start.set(t['start_datetime'])
        self.v_end.set(t['end_datetime'])
        spd = int(t['samples_per_day'])
        match = next((s for s in ["4  (6-hourly)", "24  (hourly)", "96  (15-min)",
                                   "288  (5-min)", "1440  (1-min)"] if s.startswith(str(spd))), None)
        self.v_spd.set(match or f"{spd}")

        self.v_tsi.set(a['tsi_wm2'])
        self.v_sigext.set(a['sigma_ext_m2kg'])
        self.v_hatm.set(a['atm_thickness_km'])
        self.v_rho0.set(a['atm_density_kgm3'])

        self.v_bird_on.set(bool(b['include_bird']))
        self.v_tau038.set(b['tau_a038'])
        self.v_tau05.set(b['tau_a05'])
        self.v_ozone.set(b['ozone_atmcm'])
        self.v_water.set(b['water_cm'])
        self.on_bird_toggle()

    def _collect_params(self):
        """Reads every widget back into a nested dict matching
        dbsi_config.DEFAULT_PARAMS's shape. Raises ValueError with a plain
        message (shown in a dialog by the caller) if any field doesn't
        parse -- never lets a bad float crash the app."""
        def f(var, label):
            try:
                return float(var.get())
            except ValueError:
                raise ValueError(f"'{label}' must be a number (got {var.get()!r}).")

        params = {
            'orbital': {
                'eccentricity': f(self.v_ecc, "Eccentricity"),
                'obliquity_deg': f(self.v_obliquity, "Obliquity"),
                'omega_tilde_deg': f(self.v_omega, "Longitude of perihelion"),
                'perihelion_datetime': self.v_tp.get().strip(),
                'orbital_period_days': f(self.v_period, "Orbital period"),
            },
            'location': {
                'name': self.v_preset.get(),
                'latitude_deg': f(self.v_lat, "Latitude"),
                'longitude_deg': f(self.v_lon, "Longitude"),
                'altitude_m': f(self.v_alt, "Altitude"),
            },
            'time': {
                'start_datetime': self.v_start.get().strip(),
                'end_datetime': self.v_end.get().strip(),
                'samples_per_day': int(self.v_spd.get().split()[0]),
            },
            'atmosphere': {
                'tsi_wm2': f(self.v_tsi, "TSI"),
                'sigma_ext_m2kg': f(self.v_sigext, "sigma_ext"),
                'atm_thickness_km': f(self.v_hatm, "H_a"),
                'atm_density_kgm3': f(self.v_rho0, "rho_0"),
            },
            'bird': {
                'include_bird': bool(self.v_bird_on.get()),
                'tau_a038': f(self.v_tau038, "Turbidity @0.38um"),
                'tau_a05': f(self.v_tau05, "Turbidity @0.5um"),
                'ozone_atmcm': f(self.v_ozone, "Ozone column"),
                'water_cm': f(self.v_water, "Precipitable water"),
            },
        }
        return params

    @staticmethod
    def _validate_params(params):
        """Range-checks a _collect_params() dict before it reaches
        dbsi_core.run_dbsi(), which does no input checking of its own:
        out-of-range values there either raise an opaque exception (Pa=0,
        e=1, H_a<0), freeze the GUI (large H_a -- the transmittance table
        costs O(H_a^2)), or silently return NaN/complex/unphysical numbers
        (H_a=0, negative sigma_ext or turbidity). Returns the parsed
        (start, end, perihelion) datetime64 values; raises ValueError with
        a plain message otherwise."""
        o, l, t, a, b = (params['orbital'], params['location'], params['time'],
                         params['atmosphere'], params['bird'])

        def check(ok, msg):
            if not ok:
                raise ValueError(msg)

        def when(text, label):
            try:
                val = np.datetime64(text)
            except ValueError:
                val = np.datetime64('NaT')
            check(not np.isnat(val), f"'{label}' must be an ISO date/time such as "
                                     f"2024-01-01T00:00:00 (got {text!r}).")
            return val

        check(0.0 <= o['eccentricity'] < 1.0, "Eccentricity must be in [0, 1).")
        check(o['orbital_period_days'] > 0.0, "Orbital period must be positive.")
        check(-90.0 <= l['latitude_deg'] <= 90.0, "Latitude must be between -90 and 90 deg.")
        check(a['tsi_wm2'] >= 0.0, "TSI must not be negative.")
        check(a['sigma_ext_m2kg'] >= 0.0, "sigma_ext must not be negative.")
        check(a['atm_density_kgm3'] >= 0.0, "rho_0 must not be negative.")
        check(0.0 < a['atm_thickness_km'] <= MAX_HATM_KM,
              f"H_a must be greater than 0 and at most {MAX_HATM_KM:g} km.")
        check(t['samples_per_day'] > 0, "Resolution must be a positive number of samples per day.")
        if b['include_bird']:
            check(min(b['tau_a038'], b['tau_a05'], b['ozone_atmcm'], b['water_cm']) >= 0.0,
                  "Bird & Hulstrom turbidity, ozone and water inputs must not be negative.")
            check(l['altitude_m'] <= MAX_BIRD_ALT_M,
                  f"The Bird & Hulstrom model's standard-atmosphere pressure is only defined "
                  f"up to {MAX_BIRD_ALT_M:,.0f} m altitude.")

        start = when(t['start_datetime'], "Start")
        end = when(t['end_datetime'], "End")
        t_p = when(o['perihelion_datetime'], "Perihelion date/time")
        check(end > start, "End must be after start.")
        days = (end - start) / np.timedelta64(1, 'D')
        n = days * t['samples_per_day']
        check(n <= MAX_SAMPLES,
              f"This time range/resolution would compute {n:,.0f} points "
              f"(limit {MAX_SAMPLES:,}). Shorten the range or reduce the resolution.")
        return start, end, t_p

    # ------------------------------------------------------------------
    # Callbacks
    # ------------------------------------------------------------------
    def on_preset_change(self, event=None):
        name = self.v_preset.get()
        if name in LOCATIONS:
            loc = LOCATIONS[name]
            self.v_lat.set(loc.lat_deg)
            self.v_lon.set(loc.lon_deg)
            self.v_alt.set(loc.alt_m)

    def on_bird_toggle(self):
        state = "normal" if self.v_bird_on.get() else "disabled"
        for entry in self._bird_rows:
            entry.configure(state=state)

    def on_load_config(self):
        path = filedialog.askopenfilename(title="Load DBSI config",
                                           filetypes=[("INI config", "*.ini"), ("All files", "*.*")])
        if not path:
            return
        try:
            params = dbsi_config.load_config(path)
        except Exception as exc:
            messagebox.showerror("Load Config failed", str(exc))
            return
        self._apply_params(params)
        self.status.set(f"Loaded config: {path}")

    def on_save_config(self):
        path = filedialog.asksaveasfilename(title="Save DBSI config", defaultextension=".ini",
                                             filetypes=[("INI config", "*.ini"), ("All files", "*.*")])
        if not path:
            return
        try:
            params = self._collect_params()
            dbsi_config.save_config(path, params)
        except Exception as exc:
            messagebox.showerror("Save Config failed", str(exc))
            return
        self.status.set(f"Saved config: {path}")

    def on_run(self):
        try:
            params = self._collect_params()
            start, end, t_p = self._validate_params(params)
        except ValueError as exc:
            messagebox.showerror("Invalid input", str(exc))
            return

        o, l, t, a, b = (params['orbital'], params['location'], params['time'],
                         params['atmosphere'], params['bird'])
        try:
            t0 = time.time()
            tg = dbsi_core.make_time_grid(start, end, samples_per_day=t['samples_per_day'])

            loc = Location(name=l['name'], lon_deg=l['longitude_deg'], lat_deg=l['latitude_deg'],
                            tz_hours=0.0, alt_m=l['altitude_m'])
            t_p_jd = dbsi_core.datetime_to_jd(t_p)

            res = dbsi_core.run_dbsi(
                loc, tg,
                eccentricity=o['eccentricity'], obliquity_deg=o['obliquity_deg'],
                omega_tilde_deg=o['omega_tilde_deg'], t_p_jd=t_p_jd,
                orbital_period_days=o['orbital_period_days'],
                tsi_wm2=a['tsi_wm2'], atm_thickness_km=a['atm_thickness_km'],
                atm_density_kgm3=a['atm_density_kgm3'], sigma_ext_m2kg=a['sigma_ext_m2kg'],
                include_bird=b['include_bird'], bird_tau_a038=b['tau_a038'],
                bird_tau_a05=b['tau_a05'], bird_ozone_atmcm=b['ozone_atmcm'],
                bird_water_cm=b['water_cm'])
            elapsed = time.time() - t0
        except Exception:
            messagebox.showerror("Run failed", traceback.format_exc())
            return

        self.results = res
        self.results['_loc'] = loc
        self._refresh_plot(res, loc)
        self.status.set(f"Computed {len(tg.JD):,} points ({t['start_datetime']} to "
                         f"{t['end_datetime']}, {t['samples_per_day']} samples/day) in {elapsed:.2f} s.")

    def _refresh_plot(self, res, loc):
        for ax in (self.ax0, self.ax1):
            ax.clear()

        gd = res['GD']
        self.ax0.plot(gd, res['DBSI0'], lw=0.7, label='DBSI0 (eq. 3.1)')
        self.ax0.plot(gd, res['DBSI2'], lw=0.7, label='DBSI2 (eq. 3.2)')
        if res.get('include_bird'):
            self.ax0.plot(gd, res['DBSI0_bird'], lw=0.7, ls='--', label='DBSI0 (Bird & Hulstrom)')
            self.ax0.plot(gd, res['DBSI2_bird'], lw=0.7, ls='--', label='DBSI2 (Bird & Hulstrom)')
        self.ax0.set_ylabel('DBSI (W/m$^2$)')
        span = np.datetime_as_string(gd[[0, -1]], unit='s')
        self.ax0.set_title(f"{loc.name}: DBSI0 / DBSI2, {span[0]} to {span[1]}")
        self.ax0.grid(True)
        self.ax0.legend(fontsize=8)

        self.ax1.plot(gd, res['DBSI2'] - res['DBSI0'], lw=0.7, color='tab:green',
                       label='DBSI2 - DBSI0')
        self.ax1.set_ylabel('W/m$^2$')
        self.ax1.set_title('Contribution of higher-order harmonics (S2,S3,T1,T2)')
        self.ax1.grid(True)
        self.ax1.legend(fontsize=8)

        self.fig.autofmt_xdate()
        self.fig.tight_layout(pad=3.0)
        self.canvas.draw()

    def on_export_csv(self):
        if self.results is None:
            messagebox.showwarning("Nothing to export", "Click Run first.")
            return
        path = filedialog.asksaveasfilename(title="Export DBSI time series", defaultextension=".csv",
                                             filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        res = self.results
        cols = ['JD', 'GD', 'l_sun_approx_deg', 'decl_approx_deg', 'EOT_min',
                'Tsun_true_local_deg', 'cosz_approx', 'transmissivity', 'DBSI0', 'DBSI2']
        if res.get('include_bird'):
            cols += ['transmissivity_bird', 'DBSI0_bird', 'DBSI2_bird']
        try:
            with open(path, 'w') as fh:
                fh.write(",".join(cols) + "\n")
                gd_str = np.datetime_as_string(res['GD'], unit='s')
                n = len(res['JD'])
                for i in range(n):
                    row = []
                    for c in cols:
                        if c == 'GD':
                            row.append(gd_str[i])
                        elif c == 'JD':
                            row.append(f"{res['JD'][i]:.6f}")
                        else:
                            row.append(f"{res[c][i]:.6f}")
                    fh.write(",".join(row) + "\n")
        except Exception as exc:
            messagebox.showerror("Export failed", str(exc))
            return
        self.status.set(f"Exported {n:,} rows to {path}")

    def on_save_png(self):
        if self.results is None:
            messagebox.showwarning("Nothing to save", "Click Run first.")
            return
        path = filedialog.asksaveasfilename(title="Save plot", defaultextension=".png",
                                             filetypes=[("PNG image", "*.png"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.fig.savefig(path, dpi=150)
        except Exception as exc:
            messagebox.showerror("Save failed", str(exc))
            return
        self.status.set(f"Saved plot to {path}")

    def on_open_paper(self):
        """Opens the companion paper in the system's default PDF viewer."""
        if not os.path.isfile(PAPER_PDF):
            messagebox.showerror("Paper not found", f"Could not find the paper at:\n{PAPER_PDF}")
            return
        try:
            if sys.platform == 'win32':
                os.startfile(PAPER_PDF)
            else:
                webbrowser.open('file://' + PAPER_PDF)
        except OSError as exc:
            messagebox.showerror("Could not open paper", str(exc))
            return
        self.status.set(f"Opened {os.path.basename(PAPER_PDF)} in the default PDF viewer.")


def main():
    root = tk.Tk()
    app = DBSIApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
