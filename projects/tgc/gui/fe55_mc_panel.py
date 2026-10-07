"""55Fe Monte Carlo panel for the TGC GUI.

Runs the signal studies from tools/fe55_mc.py in a background thread (in-process,
reusing its compute API) and renders them as interactive Plotly figures embedded in a
QWebEngineView. Four studies: charge & pulse-height spectrum, collimator comparison,
position grid, and single-event inspection. Run parameters are focused (study + key
knobs); the physics parameters live in config/fe55_mc.json.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np
from PyQt5.QtCore import Qt, QThread, QUrl, pyqtSignal
from PyQt5.QtWidgets import (
    QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
    QProgressBar, QPushButton, QSpinBox, QVBoxLayout, QWidget,
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

STUDIES = {  # label -> internal key
    "Charge & pulse-height spectrum": "spectrum",
    "Collimator comparison": "collimator",
    "Position grid": "grid",
    "Single events": "single",
}
VIEWS = {  # key -> list of view labels for the view selector
    "spectrum":   ["Charge spectrum", "Pulse-height spectrum"],
    "collimator": ["Charge comparison", "Pulse-height comparison"],
    "grid":       ["Observable maps", "FWHM map", "Shape sweeps",
                   "Per-point distribution", "Corner waveforms"],
    "single":     ["Event pulses"],
}
_MAPS_CACHE: dict = {}   # (shape, gain, tcut) -> bundle (uproot read is ~seconds)


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
                def _prog(done, total):
                    self.progress.emit(done, total)
                self.log_line.emit(f"[MC] grid scan: {n} photons/point …")
                scan_res, scan_meta = mc.run_grid(bundle, phys, cfg, rng, gain, n,
                                                  progress=_prog, cancel=lambda: self._cancel)
                if self._cancel:
                    self.failed.emit("cancelled"); return
                self.progress.emit(0, 0)
                self.log_line.emit("[MC] realistic exposure for waveform samples …")
                incl = mc.run_exposure(bundle, phys, rng, min(n * 4, 8000), gain, keep_samples=8)
                data = {"kind": "grid", "scan_res": scan_res, "scan_meta": scan_meta,
                        "incl": incl, "T": bundle["T"]}

            else:  # single
                self.progress.emit(0, 0)
                self.log_line.emit(f"[MC] single-event inspection: {p['inspect_events']} events/pos …")
                recs = mc.run_inspect(bundle, phys, cfg, rng, gain, p["inspect_events"])
                data = {"kind": "single", "records": recs, "T": bundle["T"]}

            self.finished_ok.emit(data)
        except Exception as exc:  # noqa: BLE001
            import traceback
            self.failed.emit(f"{exc}\n{traceback.format_exc()}")


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


def build_figure(data, view, extra=None):
    """Dispatch a cached result + a view label to a Plotly figure."""
    extra = extra or {}
    kind = data["kind"]
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

class Fe55MCPanel(QWidget):
    """GUI panel: run controls + a QWebEngineView showing interactive Plotly figures."""

    def __init__(self, log_cb=None, parent=None):
        super().__init__(parent)
        self._log_cb = log_cb or (lambda s: None)
        self._runner: MCRunner | None = None
        self._data = None
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
        self.coll_d = QDoubleSpinBox(); self.coll_d.setRange(0.1, 10.0)
        self.coll_d.setValue(1.8); self.coll_d.setSingleStep(0.1); self.coll_d.setSuffix(" mm")
        self.coll_center = QComboBox(); self.coll_center.addItems(["both", "wire", "gap"])
        self.inspect_events = QSpinBox(); self.inspect_events.setRange(1, 8); self.inspect_events.setValue(3)
        self.shape_dir = QLineEdit(str(DEFAULT_SHAPE))
        self.gain_dir = QLineEdit(str(DEFAULT_GAIN))
        self.config = QLineEdit(str(DEFAULT_CONFIG))
        form.addRow("Study", self.study)
        form.addRow("Photons (per point for grid)", self.n_photons)
        form.addRow("Seed", self.seed)
        form.addRow("Gain model", self.gain)
        form.addRow("Collimator diameter", self.coll_d)
        form.addRow("Collimator center", self.coll_center)
        form.addRow("Events/position (single)", self.inspect_events)
        form.addRow("Shape dir", self.shape_dir)
        form.addRow("Gain dir", self.gain_dir)
        form.addRow("Config", self.config)
        root.addWidget(ctl)

        # ── run / view / status row ───────────────────────────────────────
        row = QHBoxLayout()
        self.run_btn = QPushButton("Run"); self.stop_btn = QPushButton("Stop")
        self.stop_btn.setEnabled(False)
        row.addWidget(self.run_btn); row.addWidget(self.stop_btn)
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

        self.status = QLabel("Choose a study and press Run.  "
                             "(Uses the committed scan maps; no Garfield run needed.)")
        root.addWidget(self.status)

        self.web = QWebEngineView()
        root.addWidget(self.web, stretch=1)

        # ── wiring ─────────────────────────────────────────────────────────
        self.run_btn.clicked.connect(self._on_run)
        self.stop_btn.clicked.connect(self._on_stop)
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
        self.n_photons.setEnabled(key != "single")
        self.n_photons.setValue(2000 if key == "grid" else 20000)

    # ── run / stop ──────────────────────────────────────────────────────────
    def _on_run(self):
        for d, what in ((self.shape_dir.text(), "shape"), (self.gain_dir.text(), "gain")):
            if not Path(d).exists():
                self.status.setText(f"⚠ {what} dir not found: {d}  —  generate the position "
                                    f"scan products first (see docs/manual.md §14).")
                return
        key = STUDIES[self.study.currentText()]
        params = dict(n_photons=self.n_photons.value(), seed=self.seed.value(),
                      gain=self.gain.currentText(), collimator=self.coll_d.value(),
                      collimator_center=self.coll_center.currentText(),
                      inspect_events=self.inspect_events.value(),
                      shape_dir=self.shape_dir.text(), gain_dir=self.gain_dir.text(),
                      config=self.config.text())
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
        if self._runner:
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
        self.progress.setRange(0, 1); self.progress.setValue(1)
        self.status.setText(f"Done: {self.study.currentText()}.  Explore views below.")
        key = data["kind"]
        self.view.blockSignals(True)
        self.view.clear(); self.view.addItems(VIEWS[key]); self.view.setEnabled(True)
        self.view.blockSignals(False)
        if key == "grid":
            ys_u, ds_u = _grid_axes(data["scan_meta"])
            for cb, vals in ((self.pt_depth, ys_u), (self.pt_d, ds_u)):
                cb.blockSignals(True); cb.clear()
                cb.addItems([f"{v:g}" for v in vals]); cb.blockSignals(False)
        self._on_view_changed()

    # ── view rendering ──────────────────────────────────────────────────────
    def _on_view_changed(self, *_):
        if self._data is None or not self.view.isEnabled():
            return
        view = self.view.currentText()
        is_pt = (self._data["kind"] == "grid" and view == "Per-point distribution")
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
        self.web.load(QUrl.fromLocalFile(str(self._html)))
