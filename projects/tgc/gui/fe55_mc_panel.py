"""55Fe Monte Carlo panel for the TGC GUI.

Runs the signal studies from tools/fe55_mc.py in a background thread (in-process,
reusing its compute API) and renders them as interactive Plotly figures embedded in a
QWebEngineView. Four studies: charge & pulse-height spectrum, collimator comparison,
position grid, and single-event inspection. Run parameters are focused (study + key
knobs); the physics parameters live in config/fe55_mc.json.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np
from PyQt5.QtCore import Qt, QThread, QUrl, pyqtSignal
from PyQt5.QtWidgets import (
    QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox, QFormLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox, QProgressBar, QPushButton,
    QSpinBox, QVBoxLayout, QWidget,
)
from PyQt5.QtWebEngineWidgets import QWebEngineView

import plotly.graph_objects as go
from plotly.subplots import make_subplots

HERE = Path(__file__).resolve().parent
TGC_DIR = HERE.parent
sys.path.insert(0, str(TGC_DIR / "tools"))
import fe55_mc as mc  # noqa: E402 — the refactored compute API

DEFAULT_SHAPE = TGC_DIR / "results" / "scan_waveform_halfpitch"
DEFAULT_GAIN = TGC_DIR / "results" / "gain_scan_halfpitch"
DEFAULT_CONFIG = TGC_DIR / "config" / "fe55_mc.json"
SCAN_GAIN_TEMPLATE = TGC_DIR / "config" / "scan_gain.json"      # ion off (gain maps)
SCAN_WF_TEMPLATE = TGC_DIR / "config" / "scan_waveform.json"    # ion on (shape/drift)
GAIN_SCAN_PY = TGC_DIR / "tools" / "gain_scan.py"
TGC_BINARY = TGC_DIR / "build" / "tgc_sim"
GARFIELD_INSTALL = (TGC_DIR / "../../local/garfield").resolve()

STUDIES = {  # label -> internal key
    "Charge & pulse-height spectrum": "spectrum",
    "Collimator comparison": "collimator",
    "Position grid": "grid",
    "Single events": "single",
    "Underlying data (maps)": "underlying",
}
VIEWS = {  # key -> list of view labels for the view selector
    "spectrum":   ["Charge spectrum", "Pulse-height spectrum"],
    "collimator": ["Charge comparison", "Pulse-height comparison"],
    "grid":       ["Observable maps", "FWHM map", "Shape sweeps",
                   "Per-point distribution", "Corner waveforms"],
    "single":     ["Event pulses"],
    "underlying": ["Gain @ point", "Drift-time map", "Single-electron shape"],
}
_MAPS_CACHE: dict = {}   # (shape, gain, tcut) -> bundle (uproot read is ~seconds)

# Config-signature fields that must match between the detector config and the maps.
_SIG_GEOM = ("wire_pitch_cm", "gap_cm", "wire_diameter_um", "n_wires", "wire_voltage_V")
_SIG_GAS = ("gas1", "gas1_fraction_pct", "gas2", "temperature_K", "pressure_Torr")


def _cfg_signature(cfg: dict) -> dict:
    """Detector-config signature (geometry + gas) used to match a config to its maps."""
    g, gas = cfg.get("geometry", {}), cfg.get("gas", {})
    sig = {k: g.get(k) for k in _SIG_GEOM}
    sig.update({k: gas.get(k) for k in _SIG_GAS})
    return sig


def _maps_signature(gain_dir):
    """Signature of the config the maps were generated with (any point's run_config.json)."""
    for rc in sorted(Path(gain_dir).glob("**/run_config.json")):
        try:
            return _cfg_signature(json.loads(rc.read_text()))
        except Exception:  # noqa: BLE001
            continue
    return None


def _sig_diffs(a: dict, b: dict):
    """Fields where signatures a (config) and b (maps) differ, with float tolerance."""
    out = []
    for k in list(_SIG_GEOM) + list(_SIG_GAS):
        va, vb = a.get(k), b.get(k)
        if isinstance(va, (int, float)) and isinstance(vb, (int, float)):
            if abs(float(va) - float(vb)) > 1e-6:
                out.append((k, va, vb))
        elif str(va).lower() != str(vb).lower():
            out.append((k, va, vb))
    return out


def _gain_stats(gain_dir):
    """Minimum measured events/cell in the gain maps, from gain_map.csv; None if absent
    (the pooled n_events column written by gain_scan.py). Used to gate MC runs on stats."""
    import pandas as pd
    f = Path(gain_dir) / "gain_map.csv"
    if not f.exists():
        return None
    try:
        return int(pd.read_csv(f)["n_events"].min())
    except Exception:  # noqa: BLE001
        return None


# ───────────────────────────── background runner ───────────────────────────

class MCRunner(QThread):
    """Runs one study's computation in a worker thread and emits the result data."""
    log_line = pyqtSignal(str)
    progress = pyqtSignal(int, int)      # done, total (0,0 = indeterminate)
    finished_ok = pyqtSignal(object)     # the result-data dict
    failed = pyqtSignal(str)

    def __init__(self, study: str, params: dict, parent=None):
        super().__init__(parent)
        self._study = study
        self._p = params
        self._cancel = False

    def cancel(self):
        self._cancel = True

    def run(self):
        try:
            p = self._p
            cfg_path = Path(p["config"])
            cfg = mc._merge(mc.DEFAULTS,
                            json.loads(cfg_path.read_text()) if cfg_path.exists() else {})
            tcut = cfg["run"]["tcut_ns"]
            key = (str(p["shape_dir"]), str(p["gain_dir"]), float(tcut))
            if key not in _MAPS_CACHE:
                self.log_line.emit(f"[MC] loading maps: {Path(p['shape_dir']).name} / "
                                   f"{Path(p['gain_dir']).name} …")
                _MAPS_CACHE[key] = mc.load_maps(p["shape_dir"], p["gain_dir"], cfg, tcut)
            bundle = _MAPS_CACHE[key]
            phys = mc.build_phys(cfg, bundle["wires"])
            rng = np.random.default_rng(p["seed"])
            gain = p["gain"]
            n = p["n_photons"]

            if self._study == "spectrum":
                self.progress.emit(0, 0)
                self.log_line.emit(f"[MC] realistic exposure: {n} photons …")
                incl = mc.run_exposure(bundle, phys, rng, n, gain, keep_samples=8)
                data = {"kind": "spectrum", "incl": incl}

            elif self._study == "collimator":
                self.progress.emit(0, 0)
                wire0 = bundle["wires"][len(bundle["wires"]) // 2]
                pitch = cfg["geometry"]["wire_pitch_cm"]
                r = p["collimator"] / 2.0 / 10.0
                exposures = [("full-cell uniform",
                              mc.run_exposure(bundle, phys, rng, n, gain, keep_samples=0))]
                centers = (["wire", "gap"] if p["collimator_center"] == "both"
                           else [p["collimator_center"]])
                for ctr in centers:
                    self.log_line.emit(f"[MC] collimated on-{ctr}: {n} photons …")
                    cx = wire0 + (pitch / 2 if ctr == "gap" else 0.0)
                    coll = {"r_cm": r, "center_cm": cx}
                    res = mc.run_exposure(bundle, phys, rng, n, gain, collimator=coll, keep_samples=0)
                    exposures.append((f"collimated {p['collimator']:g}mm on-{ctr}", res))
                data = {"kind": "collimator", "exposures": exposures}

            elif self._study == "grid":
                g = p["grid"]
                depths = np.unique(np.linspace(g["dmin"], g["dmax"], g["ndepth"]))
                dists = np.unique(np.linspace(g["xmin"], g["xmax"], g["ndist"]))

                def _prog(done, total):
                    self.progress.emit(done, total)
                self.log_line.emit(f"[MC] grid scan: {len(depths)}×{len(dists)}="
                                   f"{len(depths) * len(dists)} points × {n} photons/point …")
                scan_res, scan_meta = mc.run_grid(bundle, phys, cfg, rng, gain, n,
                                                  depths=depths, dists=dists,
                                                  progress=_prog, cancel=lambda: self._cancel)
                if self._cancel:
                    self.failed.emit("cancelled"); return
                self.progress.emit(0, 0)
                self.log_line.emit("[MC] realistic exposure for waveform samples …")
                incl = mc.run_exposure(bundle, phys, rng, min(n * 4, 8000), gain, keep_samples=8)
                data = {"kind": "grid", "scan_res": scan_res, "scan_meta": scan_meta,
                        "incl": incl, "T": bundle["T"]}

            elif self._study == "underlying":
                self.progress.emit(0, 0)
                self.log_line.emit("[MC] loaded underlying measured maps.")
                data = {"kind": "underlying", "bundle": bundle,
                        "sig": _maps_signature(p["gain_dir"])}

            else:  # single
                self.progress.emit(0, 0)
                self.log_line.emit(f"[MC] single-event inspection: {p['inspect_events']} events/pos …")
                recs = mc.run_inspect(bundle, phys, cfg, rng, gain, p["inspect_events"])
                data = {"kind": "single", "records": recs, "T": bundle["T"]}

            self.finished_ok.emit(data)
        except Exception as exc:  # noqa: BLE001
            import traceback
            self.failed.emit(f"{exc}\n{traceback.format_exc()}")


# ─────────────────────── on-demand Garfield map generation ─────────────────

class GarfieldScanRunner(QThread):
    """Runs the tgc_sim scans (gain + waveform) that produce the MC maps for a given
    detector config and position grid, by driving tools/gain_scan.py as a subprocess."""
    log_line = pyqtSignal(str)
    progress = pyqtSignal(int, int)
    finished_ok = pyqtSignal()
    failed = pyqtSignal(str)

    def __init__(self, cfg_dict, depths, xpos, gain_dir, shape_dir,
                 n_gain, n_wf, jobs, parent=None):
        super().__init__(parent)
        self._cfg = cfg_dict
        self._depths = [float(d) for d in depths]
        self._xpos = [float(x) for x in xpos]
        self._gain_dir, self._shape_dir = str(gain_dir), str(shape_dir)
        self._n_gain, self._n_wf, self._jobs = int(n_gain), int(n_wf), int(jobs)
        self._cancel = False
        self._proc = None

    def cancel(self):
        self._cancel = True
        if self._proc and self._proc.poll() is None:
            self._proc.terminate()

    def _mkcfg(self, template, n_events, label):
        base = json.loads(Path(template).read_text())
        for blk in ("geometry", "gas"):           # adopt the detector config
            base.setdefault(blk, {}).update(self._cfg.get(blk, {}))
        base["source"]["source_distances_mm"] = self._depths
        base["source"]["x_positions_cm"] = self._xpos
        base["simulation"]["n_events"] = int(n_events)
        tf = tempfile.NamedTemporaryFile(mode="w", suffix=f"_{label}.json", delete=False)
        json.dump(base, tf, indent=2); tf.close()
        return tf.name

    def _env(self):
        env = os.environ.copy()
        env.setdefault("GARFIELD_INSTALL", str(GARFIELD_INSTALL))
        env.setdefault("HEED_DATABASE", str(GARFIELD_INSTALL / "share" / "Heed" / "database"))
        return env

    def _stream(self, cmd):
        self.log_line.emit("[scan] $ " + " ".join(str(c) for c in cmd))
        self._proc = subprocess.Popen([str(c) for c in cmd], cwd=str(TGC_DIR),
                                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                      text=True, bufsize=1, env=self._env())
        for line in self._proc.stdout:
            self.log_line.emit(line.rstrip())
            if self._cancel:
                self._proc.terminate(); break
        self._proc.wait()
        return self._proc.returncode

    def run(self):
        import gain_scan as gs
        tmp = []
        try:
            if not TGC_BINARY.exists():
                self.failed.emit(f"tgc_sim not built at {TGC_BINARY}"); return
            self.progress.emit(0, 0)
            gain_cfg = self._mkcfg(SCAN_GAIN_TEMPLATE, self._n_gain, "gain"); tmp.append(gain_cfg)
            wf_cfg = self._mkcfg(SCAN_WF_TEMPLATE, self._n_wf, "wf"); tmp.append(wf_cfg)

            # Serial gas-table warmup if missing (so parallel jobs don't race to build it).
            gas_file = TGC_DIR / gs._derive_gas_filename(json.loads(Path(gain_cfg).read_text())["gas"])
            if not gas_file.exists():
                self.log_line.emit(f"[scan] generating Magboltz gas table {gas_file.name} "
                                   f"(one-time, slow) …")
                warm = self._mkcfg(SCAN_GAIN_TEMPLATE, 1, "warm"); tmp.append(warm)
                if self._stream([TGC_BINARY, "--config", warm, "--distance", f"{self._depths[0]:g}",
                                 "--run-name", "warmup", "--out", tempfile.mkdtemp()]) != 0:
                    self.failed.emit("gas-table warmup failed"); return
            if self._cancel:
                self.failed.emit("cancelled"); return

            label = datetime.now().strftime("gui_%Y%m%d_%H%M%S")
            self.log_line.emit(f"[scan] GAIN maps → {self._gain_dir} …")
            if self._stream(["python3", GAIN_SCAN_PY, "run", "--config", gain_cfg,
                             "--out", self._gain_dir, "--label", label,
                             "--jobs", str(self._jobs)]) != 0 or self._cancel:
                self.failed.emit("cancelled" if self._cancel else "gain scan failed"); return
            self.log_line.emit(f"[scan] WAVEFORM maps → {self._shape_dir} …")
            if self._stream(["python3", GAIN_SCAN_PY, "run", "--config", wf_cfg,
                             "--out", self._shape_dir, "--label", label,
                             "--jobs", str(self._jobs)]) != 0 or self._cancel:
                self.failed.emit("cancelled" if self._cancel else "waveform scan failed"); return
            if self._stream(["python3", GAIN_SCAN_PY, "waveform",
                             "--out", self._shape_dir]) != 0:
                self.failed.emit("waveform analysis failed"); return

            _MAPS_CACHE.clear()          # new maps on disk → force a reload next MC run
            self.finished_ok.emit()
        except Exception as exc:  # noqa: BLE001
            import traceback
            self.failed.emit(f"{exc}\n{traceback.format_exc()}")
        finally:
            for f in tmp:
                try:
                    os.unlink(f)
                except OSError:
                    pass


# ─────────────────────────────── plotly builders ───────────────────────────

_C_ALL, _C_MAIN, _C_ESC = "#9e9e9e", "#1f77b4", "#d62728"
_XLABEL = {  # Plotly-friendly (HTML sub/sup) axis labels
    "charge": "collected charge  Q = Σg<sub>i</sub>  [10<sup>3</sup> gain-electrons]",
    "amplitude": "peak amplitude  max<sub>t</sub> i(t)  [10<sup>6</sup> arb.]",
}


def _hist_spectrum(res, quantity):
    _, scale, _, titleword, unit, cfmt = mc._QUANTITY[quantity]
    key = mc._QUANTITY[quantity][0]
    V = res[key] / scale
    tr = res["truth"]
    main, esc = V[tr == "main"], V[tr == "escape"]
    hi = float(np.percentile(V, 99.5)) * 1.1
    xb = dict(start=0.0, end=hi, size=hi / 120.0)
    cen = float(np.mean(main)) if len(main) else float("nan")
    sig = float(np.std(main) / np.mean(main) * 100) if len(main) else float("nan")
    cesc = float(np.mean(esc)) if len(esc) else float("nan")
    fig = go.Figure()
    fig.add_histogram(x=V, xbins=xb, marker_color=_C_ALL, name=f"all (N={len(V)})", opacity=0.7)
    fig.add_histogram(x=main, xbins=xb, marker_color=_C_MAIN, name="main 5.9 keV", marker_line_width=0)
    if len(esc):
        fig.add_histogram(x=esc, xbins=xb, marker_color=_C_ESC, name="Ar escape")
    fig.add_vline(x=cen, line=dict(color=_C_MAIN, dash="dot"))
    if len(esc):
        fig.add_vline(x=cesc, line=dict(color=_C_ESC, dash="dot"))
    # secondary top axis: self-calibrated energy (main peak = 5.9 keV)
    kev = 5.895 / cen if cen > 0 else float("nan")
    fig.add_scatter(x=[0.0, hi * kev], y=[0, 0], mode="markers", marker=dict(opacity=0),
                    xaxis="x2", showlegend=False, hoverinfo="skip")
    fig.update_layout(barmode="overlay", template="plotly_white",
                      title=f"55Fe {titleword}: main centroid={cen:{cfmt}}{unit}, "
                            f"σ/E={sig:.1f}%, esc/main={cesc/cen:.2f}",
                      xaxis=dict(title=_XLABEL[quantity]), yaxis_title="events / bin", bargap=0,
                      xaxis2=dict(overlaying="x", side="top", range=[0, hi * kev],
                                  title="energy [keV]  (self-calibrated to 5.9 keV peak)"))
    fig.update_traces(opacity=0.75, selector=dict(type="histogram"))
    return fig


def _hist_compare(exposures, quantity):
    key, scale, _, titleword, unit, cfmt = mc._QUANTITY[quantity]
    Vs = [(lab, res[key] / scale) for lab, res in exposures]
    hi = max(float(np.percentile(V, 99.5)) for _, V in Vs) * 1.1
    xb = dict(start=0.0, end=hi, size=hi / 120.0)
    cen0 = float(np.mean(Vs[0][1][exposures[0][1]["truth"] == "main"]))
    colors = ["#616161", "#ff7f0e", "#1f77b4", "#2ca02c"]
    fig = go.Figure()
    for (lab, V), (_, res), c in zip(Vs, exposures, colors):
        main = V[res["truth"] == "main"]
        cen = float(np.mean(main)); sig = float(np.std(main) / np.mean(main) * 100)
        fig.add_histogram(x=V, xbins=xb, histnorm="probability density", name=
                          f"{lab}: {cen:{cfmt}}{unit}, σ/E={sig:.1f}%, ×{cen/cen0:.2f}",
                          marker_color=c, opacity=0.55)
        fig.add_vline(x=cen, line=dict(color=c, dash="dot"))
    fig.update_layout(barmode="overlay", template="plotly_white",
                      title=f"55Fe {titleword}: collimator-position comparison",
                      xaxis_title=_XLABEL[quantity], yaxis_title="probability density", bargap=0)
    return fig


def _grid_axes(meta):
    ys_u = sorted({m["depth_mm"] for m in meta.values()})
    ds_u = sorted({m["d_wire_mm"] for m in meta.values()})
    return ys_u, ds_u


def _grid_z(meta, scan_res, ys_u, ds_u, valfn):
    Z = np.full((len(ys_u), len(ds_u)), np.nan)
    for name, m in meta.items():
        r = scan_res[name]
        if not r["n"]:
            continue
        Z[ys_u.index(m["depth_mm"]), ds_u.index(m["d_wire_mm"])] = valfn(r)
    return Z


def _fig_grid_maps(scan_res, meta):
    ys_u, ds_u = _grid_axes(meta)
    panels = [("Peak time [ns]", lambda r: float(np.nanmean(r["tpeak"])), "Viridis"),
              ("Rise 10-90% [ns]", lambda r: float(np.nanmean(r["rise"])), "Magma"),
              ("Mean amplitude [arb.]", lambda r: float(r["amp"].mean()), "Cividis"),
              ("Mean charge [1e3]", lambda r: float(r["Q"].mean()) / 1e3, "Plasma")]
    fig = make_subplots(rows=2, cols=2, subplot_titles=[p[0] for p in panels],
                        horizontal_spacing=0.12, vertical_spacing=0.14)
    for i, (title, vf, cs) in enumerate(panels):
        Z = _grid_z(meta, scan_res, ys_u, ds_u, vf)
        r, c = i // 2 + 1, i % 2 + 1
        fig.add_heatmap(z=Z, x=ds_u, y=ys_u, colorscale=cs, row=r, col=c,
                        colorbar=dict(len=0.42, y=0.79 if r == 1 else 0.21,
                                      x=0.46 if c == 1 else 1.0))
    fig.update_xaxes(title_text="distance from wire [mm]")
    fig.update_yaxes(title_text="depth [mm]")
    fig.update_layout(template="plotly_white", title="55Fe pulse observables across the grid")
    return fig


def _fig_fwhm(scan_res, meta):
    ys_u, ds_u = _grid_axes(meta)
    Z = _grid_z(meta, scan_res, ys_u, ds_u,
                lambda r: mc._fwhm_ns(r["t_rel"], r["aligned_mean"]))
    fig = go.Figure(go.Heatmap(z=Z, x=ds_u, y=ys_u, colorscale="Magma",
                               colorbar=dict(title="FWHM [ns]")))
    fig.update_layout(template="plotly_white", title="55Fe pulse FWHM vs position",
                      xaxis_title="distance from wire [mm]  (0 = on wire)",
                      yaxis_title="conversion depth [mm]")
    return fig


def _sample_color(frac, scale):
    import plotly.colors as pc
    return pc.sample_colorscale(scale, max(0.0, min(1.0, frac)))[0]


def _fig_sweep(scan_res, meta):
    ys_u, ds_u = _grid_axes(meta)
    depth_fixed, d_fixed = ys_u[len(ys_u) // 2], ds_u[0]

    def _find(depth, d):
        for name, m in meta.items():
            if abs(m["depth_mm"] - depth) < 1e-6 and abs(m["d_wire_mm"] - d) < 1e-6:
                return scan_res[name]
        return None

    fig = make_subplots(rows=1, cols=2, subplot_titles=(
        f"Transverse sweep at depth = {depth_fixed:g} mm",
        f"Depth sweep at distance = {d_fixed:g} mm"), horizontal_spacing=0.1)
    for dv in ds_u:
        r = _find(depth_fixed, dv)
        if r and r["n"]:
            fig.add_scatter(x=r["t_rel"], y=r["aligned_mean"], mode="lines",
                            line=dict(color=_sample_color(dv / ds_u[-1], "Viridis"), width=1.3),
                            name=f"d={dv:g}mm", row=1, col=1, showlegend=False)
    for yv in ys_u:
        r = _find(yv, d_fixed)
        if r and r["n"]:
            fig.add_scatter(x=r["t_rel"], y=r["aligned_mean"], mode="lines",
                            line=dict(color=_sample_color((yv - ys_u[0]) / (ys_u[-1] - ys_u[0]),
                                                          "Plasma"), width=1.3),
                            name=f"y={yv:g}mm", row=1, col=2, showlegend=False)
    fig.update_xaxes(title_text="time since peak [ns]", range=[-5, 80])
    fig.update_yaxes(title_text="mean pulse / peak")
    fig.update_layout(template="plotly_white",
                      title="Aligned pulse shapes: transverse (Viridis: wire→gap) "
                            "and depth (Plasma: shallow→deep) sweeps")
    return fig


def _point_result(scan_res, meta, depth, d):
    for name, m in meta.items():
        if abs(m["depth_mm"] - depth) < 1e-6 and abs(m["d_wire_mm"] - d) < 1e-6:
            return scan_res[name]
    return None


def _fig_point_dist(scan_res, meta, depth, d):
    r = _point_result(scan_res, meta, depth, d)
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Charge", "Peak amplitude"),
                        horizontal_spacing=0.12)
    if r is None or not r["n"]:
        return fig
    for col, quantity in ((1, "charge"), (2, "amplitude")):
        key, scale, xlabel, _, unit, cfmt = mc._QUANTITY[quantity]
        V = r[key] / scale
        tr = r["truth"]
        hi = float(np.percentile(V, 99.5)) * 1.15
        xb = dict(start=0.0, end=hi, size=hi / 60.0)
        fig.add_histogram(x=V, xbins=xb, marker_color=_C_ALL, opacity=0.7,
                          name="all", row=1, col=col, showlegend=False)
        fig.add_histogram(x=V[tr == "main"], xbins=xb, marker_color=_C_MAIN,
                          name="main", row=1, col=col, showlegend=(col == 1))
        if (tr == "escape").any():
            fig.add_histogram(x=V[tr == "escape"], xbins=xb, marker_color=_C_ESC,
                              name="escape", row=1, col=col, showlegend=(col == 1))
        short = "charge [10³ e⁻]" if quantity == "charge" else "amplitude [10⁶ arb.]"
        fig.update_xaxes(title_text=short, row=1, col=col)
    fig.update_layout(barmode="overlay", template="plotly_white", bargap=0,
                      title=f"Per-point distributions at depth={depth:g} mm, "
                            f"distance-from-wire={d:g} mm")
    return fig


def _fig_corner_waveforms(scan_res, meta, incl, T):
    ys_u, ds_u = _grid_axes(meta)
    corners = [("on-wire shallow", ys_u[0], ds_u[0]), ("on-wire deep", ys_u[-1], ds_u[0]),
               ("mid-gap shallow", ys_u[0], ds_u[-1]), ("mid-gap deep", ys_u[-1], ds_u[-1])]
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.1,
                        subplot_titles=("Mean pulse vs position (peak-aligned)",
                                        "Realistic-exposure pulses (peak-aligned)"))
    for lab, yv, dv in corners:
        r = _point_result(scan_res, meta, yv, dv)
        if r and r["n"]:
            fig.add_scatter(x=r["t_rel"], y=r["aligned_mean"], mode="lines",
                            name=lab, row=1, col=1)
    if incl and incl["n"]:
        for W in incl["samples"][:8]:
            pk = int(np.argmax(W))
            if W[pk] > 0:
                fig.add_scatter(x=T - T[pk], y=W / W[pk], mode="lines", opacity=0.3,
                                line=dict(color=_C_MAIN, width=0.8), showlegend=False,
                                row=1, col=2)
        fig.add_scatter(x=incl["t_rel"], y=incl["aligned_mean"], mode="lines",
                        line=dict(color="black", width=2), name="aligned mean", row=1, col=2)
    fig.update_xaxes(title_text="time since peak [ns]", range=[-5, 80])
    fig.update_yaxes(title_text="pulse / peak")
    fig.update_layout(template="plotly_white", title="55Fe pulse shapes")
    return fig


def _fig_single(records, T):
    positions = [p for p in records if records[p]]
    ncol = max((len(records[p]) for p in positions), default=1)
    specs = [[{"secondary_y": True} for _ in range(ncol)] for _ in positions]
    fig = make_subplots(rows=len(positions), cols=ncol, specs=specs,
                        vertical_spacing=0.07, horizontal_spacing=0.09)
    for i, pos in enumerate(positions):
        for j, rec in enumerate(records[pos][:ncol]):
            W, taus, g = rec["W"], rec["taus"], rec["g"]
            conv, tr = rec["conv"], rec["truth"]
            tp = T[int(np.argmax(W))]
            fig.add_scatter(x=T, y=W, mode="lines", line=dict(color=_C_MAIN, width=1.2),
                            showlegend=False, row=i + 1, col=j + 1, secondary_y=False)
            fig.add_scatter(x=taus, y=g, mode="markers",
                            marker=dict(color="#ff7f0e", size=3, opacity=0.5),
                            showlegend=False, row=i + 1, col=j + 1, secondary_y=True)
            fig.add_annotation(row=i + 1, col=j + 1, x=0.96, y=0.95, xref="x domain",
                               yref="y domain", align="right", showarrow=False, font=dict(size=8),
                               text=(f"<b>{pos}</b><br>N={len(g)}, Q={g.sum()/1e3:.0f}k<br>"
                                     f"{tr}, pk={tp:.1f}ns"))
            fig.update_xaxes(range=[max(T[0], tp - 8), tp + 45], row=i + 1, col=j + 1)
            fig.update_yaxes(showticklabels=(j + 1 == ncol), secondary_y=True,
                             row=i + 1, col=j + 1)
    fig.update_layout(template="plotly_white", height=250 * len(positions),
                      title="55Fe single events: W(t) [blue] + per-primary gains [orange]")
    return fig


def _nearest_cell(bundle, depth, d):
    ys, ds = bundle["ys"], bundle["ds"]
    return int(np.argmin(np.abs(ys - depth))), int(np.argmin(np.abs(ds - d)))


def _fig_underlying_gain(bundle, depth, d):
    iy, idd = _nearest_cell(bundle, depth, d)
    arr = np.asarray(bundle["maps"]["cells"][iy, idd], float)
    yv, dv = bundle["ys"][iy], bundle["ds"][idd]
    fig = go.Figure(go.Histogram(x=arr, nbinsx=60, marker_color=_C_MAIN))
    fig.add_vline(x=float(arr.mean()), line=dict(color="black", dash="dot"))
    zero = float(np.mean(arr <= 0) * 100)
    fig.update_layout(template="plotly_white", bargap=0,
                      title=f"Measured single-electron gain (avalanche_size) at depth={yv:g} mm, "
                            f"d={dv:g} mm<br>N={len(arr)}, mean={arr.mean():.0f}, "
                            f"RMS={arr.std():.0f}, zero-gain={zero:.0f}%  "
                            f"— the Polya the MC samples at this cell",
                      xaxis_title="avalanche size (electrons)", yaxis_title="events / bin")
    return fig


def _fig_underlying_drift(bundle):
    ys, ds = bundle["ys"], bundle["ds"]
    tau, rms = bundle["maps"]["tau"], bundle["maps"]["rms"]
    fig = go.Figure(go.Heatmap(z=tau, x=ds, y=ys, colorscale="Viridis",
                               colorbar=dict(title="mean drift [ns]"),
                               customdata=rms,
                               hovertemplate="d=%{x} mm, depth=%{y} mm<br>"
                                             "mean=%{z:.2f} ns, RMS=%{customdata:.2f} ns<extra></extra>"))
    fig.update_layout(template="plotly_white", title="Measured electron drift time to the wire",
                      xaxis_title="distance from wire [mm]", yaxis_title="conversion depth [mm]")
    return fig


def _fig_underlying_shape(bundle):
    fig = go.Figure(go.Scatter(x=bundle["ts"], y=bundle["S"], mode="lines",
                               line=dict(color=_C_MAIN, width=1.6)))
    fig.add_vline(x=0.0, line=dict(color="gray", dash="dot"))
    fig.update_layout(template="plotly_white",
                      title="Measured single-electron anode pulse shape S(t′) (position-independent)",
                      xaxis_title="time since peak [ns]", yaxis_title="peak-normalized current")
    return fig


def build_figure(data, view, extra=None):
    """Dispatch a cached result + a view label to a Plotly figure."""
    extra = extra or {}
    kind = data["kind"]
    if kind == "underlying":
        b = data["bundle"]
        if view == "Drift-time map":
            return _fig_underlying_drift(b)
        if view == "Single-electron shape":
            return _fig_underlying_shape(b)
        return _fig_underlying_gain(b, extra.get("depth", b["ys"][0]), extra.get("d", b["ds"][0]))
    if kind == "spectrum":
        return _hist_spectrum(data["incl"],
                              "charge" if view.startswith("Charge") else "amplitude")
    if kind == "collimator":
        return _hist_compare(data["exposures"],
                             "charge" if view.startswith("Charge") else "amplitude")
    if kind == "grid":
        sr, sm, incl, T = data["scan_res"], data["scan_meta"], data["incl"], data["T"]
        if view == "Observable maps":
            return _fig_grid_maps(sr, sm)
        if view == "FWHM map":
            return _fig_fwhm(sr, sm)
        if view == "Shape sweeps":
            return _fig_sweep(sr, sm)
        if view == "Per-point distribution":
            return _fig_point_dist(sr, sm, extra.get("depth"), extra.get("d"))
        return _fig_corner_waveforms(sr, sm, incl, T)
    return _fig_single(data["records"], data["T"])


# ──────────────────────────────── the panel ────────────────────────────────

class PlotWindow(QMainWindow):
    """Separate top-level window holding the interactive Plotly figure."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("55Fe MC — plot")
        self.resize(1000, 700)
        self.web = QWebEngineView()
        self.setCentralWidget(self.web)

    def show_html(self, path):
        self.web.load(QUrl.fromLocalFile(str(path)))
        self.show(); self.raise_(); self.activateWindow()


class Fe55MCPanel(QWidget):
    """GUI panel: run controls; interactive Plotly figures open in a separate window."""

    def __init__(self, log_cb=None, parent=None):
        super().__init__(parent)
        self._log_cb = log_cb or (lambda s: None)
        self._runner: MCRunner | None = None
        self._gscan = None                 # GarfieldScanRunner while a scan is running
        self._data = None
        self._plotwin: PlotWindow | None = None
        self.config_panel = None           # set by MainWindow; the live detector config
        self._html = Path(tempfile.gettempdir()) / f"fe55_mc_{id(self)}.html"

        root = QVBoxLayout(self)

        # ── run controls ──────────────────────────────────────────────────
        ctl = QGroupBox("55Fe Monte Carlo — run")
        form = QFormLayout(ctl)
        self.study = QComboBox(); self.study.addItems(STUDIES.keys())
        self.n_photons = QSpinBox(); self.n_photons.setRange(100, 500000)
        self.n_photons.setValue(20000); self.n_photons.setSingleStep(1000)
        self.seed = QSpinBox(); self.seed.setRange(0, 1_000_000); self.seed.setValue(1)
        self.gain = QComboBox(); self.gain.addItems(["empirical", "polya"])
        self.g_stats = QSpinBox(); self.g_stats.setRange(10, 100000); self.g_stats.setValue(1000)
        self.g_stats.setSingleStep(100)
        self.g_stats.setToolTip("Required measured events per cell in the gain (avalanche_size) "
                                "distribution. A run prompts to produce more via Garfield when the "
                                "maps are below this (or missing).")
        self.coll_d = QDoubleSpinBox(); self.coll_d.setRange(0.1, 10.0)
        self.coll_d.setValue(1.8); self.coll_d.setSingleStep(0.1); self.coll_d.setSuffix(" mm")
        self.coll_center = QComboBox(); self.coll_center.addItems(["both", "wire", "gap"])
        self.inspect_events = QSpinBox(); self.inspect_events.setRange(1, 8); self.inspect_events.setValue(3)
        # position-grid parameters (resolution + extent); defaults reproduce the native 5x10 grid
        _gtip = ("Position-grid resolution/extent. Drift & timing are interpolated; the GAIN is sampled "
                 "from the nearest measured cell, so finer grids leave the charge/amplitude maps stepwise "
                 "at the native resolution.")
        self.g_ndepth = QSpinBox(); self.g_ndepth.setRange(2, 21); self.g_ndepth.setValue(5)
        self.g_ndist = QSpinBox(); self.g_ndist.setRange(2, 41); self.g_ndist.setValue(10)
        self.g_dmin = QDoubleSpinBox(); self.g_dmin.setRange(0.1, 1.3); self.g_dmin.setDecimals(2)
        self.g_dmin.setSingleStep(0.1); self.g_dmin.setValue(0.1); self.g_dmin.setSuffix(" mm")
        self.g_dmax = QDoubleSpinBox(); self.g_dmax.setRange(0.1, 1.3); self.g_dmax.setDecimals(2)
        self.g_dmax.setSingleStep(0.1); self.g_dmax.setValue(1.3); self.g_dmax.setSuffix(" mm")
        self.g_xmin = QDoubleSpinBox(); self.g_xmin.setRange(0.0, 0.9); self.g_xmin.setDecimals(2)
        self.g_xmin.setSingleStep(0.1); self.g_xmin.setValue(0.0); self.g_xmin.setSuffix(" mm")
        self.g_xmax = QDoubleSpinBox(); self.g_xmax.setRange(0.0, 0.9); self.g_xmax.setDecimals(2)
        self.g_xmax.setSingleStep(0.1); self.g_xmax.setValue(0.9); self.g_xmax.setSuffix(" mm")
        self._grid_fields = [self.g_ndepth, self.g_ndist, self.g_dmin, self.g_dmax,
                             self.g_xmin, self.g_xmax]
        for w in self._grid_fields:
            w.setToolTip(_gtip)
        self.shape_dir = QLineEdit(str(DEFAULT_SHAPE))
        self.gain_dir = QLineEdit(str(DEFAULT_GAIN))
        self.config = QLineEdit(str(DEFAULT_CONFIG))

        def _pair(a, b):
            w = QWidget(); h = QHBoxLayout(w); h.setContentsMargins(0, 0, 0, 0)
            h.addWidget(a); h.addWidget(QLabel("–")); h.addWidget(b)
            return w

        form.addRow("Study", self.study)
        form.addRow("Photons (per point for grid)", self.n_photons)
        form.addRow("Seed", self.seed)
        form.addRow("Gain model", self.gain)
        form.addRow("Gain stats (events/cell)", self.g_stats)
        form.addRow("Collimator diameter", self.coll_d)
        form.addRow("Collimator center", self.coll_center)
        form.addRow("Events/position (single)", self.inspect_events)
        form.addRow("Grid depth points", self.g_ndepth)
        form.addRow("Grid depth range", _pair(self.g_dmin, self.g_dmax))
        form.addRow("Grid distance points", self.g_ndist)
        form.addRow("Grid distance range", _pair(self.g_xmin, self.g_xmax))
        form.addRow("Shape dir", self.shape_dir)
        form.addRow("Gain dir", self.gain_dir)
        form.addRow("Config", self.config)
        root.addWidget(ctl)

        # ── run / view / status row ───────────────────────────────────────
        row = QHBoxLayout()
        self.run_btn = QPushButton("Run"); self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        self.show_plot_btn = QPushButton("Show plot window"); self.show_plot_btn.setEnabled(False)
        row.addWidget(self.run_btn); row.addWidget(self.stop_btn); row.addWidget(self.show_plot_btn)
        row.addWidget(QLabel("View:"))
        self.view = QComboBox(); self.view.setEnabled(False)
        row.addWidget(self.view)
        self.pt_depth = QComboBox(); self.pt_d = QComboBox()
        self.pt_depth.setVisible(False); self.pt_d.setVisible(False)
        row.addWidget(QLabel("depth:")); row.addWidget(self.pt_depth)
        row.addWidget(QLabel("d:")); row.addWidget(self.pt_d)
        row.addStretch()
        self.progress = QProgressBar(); self.progress.setMaximumWidth(220)
        row.addWidget(self.progress)
        root.addLayout(row)
        # hide the per-point depth/d labels together with the combos
        self._pt_labels = [row.itemAt(i).widget() for i in range(row.count())
                           if isinstance(row.itemAt(i).widget(), QLabel)
                           and row.itemAt(i).widget().text() in ("depth:", "d:")]
        for lb in self._pt_labels:
            lb.setVisible(False)

        self.status = QLabel("Choose a study and press Run.  Figures open in a separate plot "
                             "window.  (Uses the committed scan maps; no Garfield run needed.)")
        self.status.setWordWrap(True)
        root.addWidget(self.status)
        root.addStretch(1)

        # ── wiring ─────────────────────────────────────────────────────────
        self.run_btn.clicked.connect(self._on_run)
        self.stop_btn.clicked.connect(self._on_stop)
        self.show_plot_btn.clicked.connect(self._show_plot_window)
        self.study.currentTextChanged.connect(self._sync_controls)
        self.view.currentTextChanged.connect(self._on_view_changed)
        self.pt_depth.currentTextChanged.connect(self._on_view_changed)
        self.pt_d.currentTextChanged.connect(self._on_view_changed)
        self._sync_controls()

    # ── control enabling per study ─────────────────────────────────────────
    def _sync_controls(self):
        key = STUDIES[self.study.currentText()]
        self.coll_d.setEnabled(key == "collimator")
        self.coll_center.setEnabled(key == "collimator")
        self.inspect_events.setEnabled(key == "single")
        for w in self._grid_fields:
            w.setEnabled(key == "grid")
        self.n_photons.setEnabled(key in ("spectrum", "collimator", "grid"))
        if key in ("grid", "spectrum", "collimator"):
            self.n_photons.setValue(2000 if key == "grid" else 20000)

    # ── separate plot window ────────────────────────────────────────────────
    def _ensure_plotwin(self):
        if self._plotwin is None:
            self._plotwin = PlotWindow(self.window())
        return self._plotwin

    def _show_plot_window(self):
        if self._plotwin is not None:
            self._plotwin.show(); self._plotwin.raise_(); self._plotwin.activateWindow()

    # ── data validation + on-demand Garfield generation ─────────────────────
    def _measured_points(self, shape_dir, config_path):
        """Measured (depth_mm, d_wire_mm) set from the drift map; empty if absent."""
        import pandas as pd
        f = Path(shape_dir) / "waveform_drift.csv"
        if not f.exists():
            return set()
        cfg = mc._merge(mc.DEFAULTS, json.loads(Path(config_path).read_text())
                        if Path(config_path).exists() else {})
        wires = mc._wire_positions_cm(cfg["geometry"])
        df = pd.read_csv(f)
        return {(round(float(r.source_distance_mm), 3),
                 round(float(np.min(np.abs(r.x_position_cm - wires)) * 10), 3))
                for r in df.itertuples()}

    def _require_data(self, study, params):
        """True = ok to run the MC; otherwise shows an error box (maybe launching a scan)."""
        gain_dir, shape_dir = Path(params["gain_dir"]), Path(params["shape_dir"])
        stats = int(params.get("stats", 1000))
        miss = []
        if not (shape_dir / "waveform_drift.csv").exists() \
                or not (shape_dir / "waveform_shapes.csv").exists():
            miss.append(f"{shape_dir}  (waveform_shapes.csv / waveform_drift.csv)")
        if not list(gain_dir.glob("**/tgc_sim.root")):
            miss.append(f"{gain_dir}  (gain avalanche-size roots)")
        if miss:
            return self._prompt_generate("No Garfield map data found in:\n  " + "\n  ".join(miss),
                                         None, None, allow_proceed=False, n_gain_default=stats)
        # detector-config match
        if self.config_panel is not None:
            try:
                cur = _cfg_signature(self.config_panel.to_config_dict())
            except Exception:  # noqa: BLE001
                cur = None
            msig = _maps_signature(gain_dir)
            if cur and msig and _sig_diffs(cur, msig):
                txt = ("The maps were generated for a DIFFERENT detector configuration:\n"
                       + "\n".join(f"  {k}: config={a}  vs  maps={b}"
                                   for k, a, b in _sig_diffs(cur, msig)))
                return self._prompt_generate(txt, None, None, allow_proceed=False, n_gain_default=stats)
        # gain-distribution stats (every study samples the per-cell avalanche-size Polya)
        have = _gain_stats(gain_dir)
        if have is not None and have < stats:
            txt = (f"The measured gain distribution has only {have} events/cell, below the requested "
                   f"{stats}. Produce more (Garfield, pooled incrementally into the existing maps) or "
                   f"use what's available.")
            if not self._prompt_generate(txt, None, None, allow_proceed=True, n_gain_default=stats):
                return False
        # grid coverage (grid study with custom points)
        if study == "grid":
            g = params["grid"]
            depths = np.unique(np.round(np.linspace(g["dmin"], g["dmax"], g["ndepth"]), 4))
            dreq = np.unique(np.round(np.linspace(g["xmin"], g["xmax"], g["ndist"]), 4))
            meas = self._measured_points(shape_dir, params["config"])
            tol = 0.02
            missing = [(float(y), float(d)) for y in depths for d in dreq
                       if not any(abs(y - my) < tol and abs(d - md) < tol for my, md in meas)]
            if missing:
                cfg = mc._merge(mc.DEFAULTS, json.loads(Path(params["config"]).read_text())
                                if Path(params["config"]).exists() else {})
                wires = mc._wire_positions_cm(cfg["geometry"])
                wire0 = wires[len(wires) // 2]
                xpos = sorted({round(wire0 + d / 10.0, 4) for _, d in missing})
                ydep = sorted({y for y, _ in missing})
                txt = (f"{len(missing)} of {len(depths) * len(dreq)} requested grid points are NOT "
                       f"measured. They would be interpolated (drift/timing) and taken from the "
                       f"nearest measured cell (gain).")
                return self._prompt_generate(txt, ydep, xpos, allow_proceed=True, n_gain_default=stats)
        return True

    def _prompt_generate(self, message, depths, xpos, allow_proceed, n_gain_default=100):
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Warning)
        box.setWindowTitle("Missing MC data for this configuration")
        box.setText(message)
        gen = box.addButton("Run Garfield…", QMessageBox.AcceptRole)
        proc = (box.addButton("Use what's available", QMessageBox.DestructiveRole)
                if allow_proceed else None)
        box.addButton("Cancel", QMessageBox.RejectRole)
        box.exec_()
        clicked = box.clickedButton()
        if clicked is gen:
            self._open_generate_dialog(depths, xpos, n_gain_default)
            return False
        if proc is not None and clicked is proc:
            self._log_cb("[MC] proceeding with the available (interpolated / lower-stat) maps.")
            return True
        return False

    def _open_generate_dialog(self, depths, xpos, n_gain_default=100):
        if self.config_panel is None:
            QMessageBox.critical(self, "No detector config",
                                 "The detector configuration panel is not available to this tab.")
            return
        if depths is None:                      # whole native grid, from the scan template
            t = json.loads(SCAN_GAIN_TEMPLATE.read_text())["source"]
            depths, xpos = t["source_distances_mm"], t["x_positions_cm"]
        dlg = QDialog(self); dlg.setWindowTitle("Run Garfield simulation → MC maps")
        form = QFormLayout(dlg)
        warn = QLabel(f"Runs tgc_sim at {len(depths)} depths × {len(xpos)} x = "
                      f"{len(depths) * len(xpos)} points for the CURRENT detector config, then "
                      f"rebuilds the gain + waveform maps.\nThis can take minutes–hours; a new "
                      f"gas/voltage also regenerates the Magboltz table (one-time, slow). Results "
                      f"accumulate incrementally, so partial runs aren't wasted.")
        warn.setWordWrap(True); form.addRow(warn)
        n_gain = QSpinBox(); n_gain.setRange(1, 100000); n_gain.setValue(int(n_gain_default))
        n_wf = QSpinBox(); n_wf.setRange(1, 2000); n_wf.setValue(30)
        jobs = QSpinBox(); jobs.setRange(1, 64); jobs.setValue(min(8, max(1, len(depths))))
        form.addRow("Gain events / point", n_gain)
        form.addRow("Waveform events / point", n_wf)
        form.addRow("Parallel jobs", jobs)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        form.addRow(bb); bb.accepted.connect(dlg.accept); bb.rejected.connect(dlg.reject)
        if dlg.exec_() == QDialog.Accepted:
            self._launch_garfield(depths, xpos, n_gain.value(), n_wf.value(), jobs.value())

    def _launch_garfield(self, depths, xpos, n_gain, n_wf, jobs):
        cfg = self.config_panel.to_config_dict()
        self._gscan = GarfieldScanRunner(cfg, depths, xpos, self.gain_dir.text(),
                                         self.shape_dir.text(), n_gain, n_wf, jobs)
        self._gscan.log_line.connect(self._log_cb)
        self._gscan.progress.connect(self._on_progress)
        self._gscan.finished_ok.connect(self._on_garfield_done)
        self._gscan.failed.connect(self._on_garfield_failed)
        self.run_btn.setEnabled(False); self.stop_btn.setEnabled(True)
        self.progress.setRange(0, 0)
        self.status.setText("Running Garfield scan to produce MC maps … (see Log tab)")
        self._gscan.start()

    def _on_garfield_done(self):
        self.run_btn.setEnabled(True); self.stop_btn.setEnabled(False)
        self.progress.setRange(0, 1); self.progress.setValue(1)
        self.status.setText("Garfield maps generated.  Press Run to use them.")

    def _on_garfield_failed(self, msg):
        self.run_btn.setEnabled(True); self.stop_btn.setEnabled(False)
        self.progress.setRange(0, 1); self.progress.setValue(0)
        self.status.setText("Garfield scan cancelled." if msg == "cancelled"
                            else "Garfield scan failed — see Log tab.")
        if msg != "cancelled":
            self._log_cb(f"[scan] ERROR: {msg}")

    # ── run / stop ──────────────────────────────────────────────────────────
    def _params(self):
        dmin, dmax = sorted((self.g_dmin.value(), self.g_dmax.value()))
        xmin, xmax = sorted((self.g_xmin.value(), self.g_xmax.value()))
        return dict(n_photons=self.n_photons.value(), seed=self.seed.value(),
                    gain=self.gain.currentText(), collimator=self.coll_d.value(),
                    collimator_center=self.coll_center.currentText(),
                    inspect_events=self.inspect_events.value(), stats=self.g_stats.value(),
                    grid=dict(ndepth=self.g_ndepth.value(), ndist=self.g_ndist.value(),
                              dmin=dmin, dmax=dmax, xmin=xmin, xmax=xmax),
                    shape_dir=self.shape_dir.text(), gain_dir=self.gain_dir.text(),
                    config=self.config.text())

    def _on_run(self):
        if self._gscan is not None and self._gscan.isRunning():
            return
        key = STUDIES[self.study.currentText()]
        params = self._params()
        if not self._require_data(key, params):   # validation + (maybe) run-Garfield prompt
            return
        self.run_btn.setEnabled(False); self.stop_btn.setEnabled(True)
        self.view.setEnabled(False)
        self.progress.setRange(0, 0)
        self.status.setText(f"Running: {self.study.currentText()} …")
        self._runner = MCRunner(key, params)
        self._runner.log_line.connect(self._log_cb)
        self._runner.progress.connect(self._on_progress)
        self._runner.finished_ok.connect(self._on_done)
        self._runner.failed.connect(self._on_failed)
        self._runner.start()

    def _on_stop(self):
        if self._gscan is not None and self._gscan.isRunning():
            self._gscan.cancel()
        if self._runner is not None and self._runner.isRunning():
            self._runner.cancel()
        self.status.setText("Cancelling …")

    def _on_progress(self, done, total):
        if total > 0:
            self.progress.setRange(0, total); self.progress.setValue(done)
        else:
            self.progress.setRange(0, 0)

    def _on_failed(self, msg):
        self.run_btn.setEnabled(True); self.stop_btn.setEnabled(False)
        self.progress.setRange(0, 1); self.progress.setValue(0)
        self.status.setText("Cancelled." if msg == "cancelled" else "Error — see Log tab.")
        if msg != "cancelled":
            self._log_cb(f"[MC] ERROR: {msg}")

    def _on_done(self, data):
        self._data = data
        self.run_btn.setEnabled(True); self.stop_btn.setEnabled(False)
        self.show_plot_btn.setEnabled(True)
        self.progress.setRange(0, 1); self.progress.setValue(1)
        self.status.setText(f"Done: {self.study.currentText()}.  Pick a view below; "
                            f"the figure opens in the plot window.")
        key = data["kind"]
        self.view.blockSignals(True)
        self.view.clear(); self.view.addItems(VIEWS[key]); self.view.setEnabled(True)
        self.view.blockSignals(False)
        pt_src = None
        if key == "grid":
            pt_src = _grid_axes(data["scan_meta"])
        elif key == "underlying":
            b = data["bundle"]
            pt_src = (list(b["ys"]), list(b["ds"]))
            s = data.get("sig") or {}
            if s:
                f1 = float(s.get("gas1_fraction_pct", 0))
                self.status.setText(
                    f"Maps config:  V={s.get('wire_voltage_V')} V,  {s.get('gas1')}/"
                    f"{s.get('gas2')} {f1:g}:{100 - f1:g},  pitch={s.get('wire_pitch_cm')} cm,  "
                    f"gap={s.get('gap_cm')} cm.   Pick a view (figure opens in the plot window).")
        if pt_src is not None:
            for cb, vals in ((self.pt_depth, pt_src[0]), (self.pt_d, pt_src[1])):
                cb.blockSignals(True); cb.clear()
                cb.addItems([f"{v:g}" for v in vals]); cb.blockSignals(False)
        self._on_view_changed()

    # ── view rendering ──────────────────────────────────────────────────────
    def _on_view_changed(self, *_):
        if self._data is None or not self.view.isEnabled():
            return
        view = self.view.currentText()
        is_pt = ((self._data["kind"] == "grid" and view == "Per-point distribution")
                 or (self._data["kind"] == "underlying" and view == "Gain @ point"))
        for w in (self.pt_depth, self.pt_d, *self._pt_labels):
            w.setVisible(is_pt)
        extra = {}
        if is_pt and self.pt_depth.currentText() and self.pt_d.currentText():
            extra = {"depth": float(self.pt_depth.currentText()),
                     "d": float(self.pt_d.currentText())}
        try:
            fig = build_figure(self._data, view, extra)
        except Exception as exc:  # noqa: BLE001
            self._log_cb(f"[MC] figure error: {exc}")
            return
        fig.write_html(str(self._html), include_plotlyjs=True, full_html=True)
        self._ensure_plotwin().show_html(self._html)
