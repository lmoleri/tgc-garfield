#!/usr/bin/env python3
"""
TGC Simulation GUI
A PyQt5 desktop application for configuring, running, and displaying results
from the tgc_sim Garfield++ binary.

Launch from anywhere:
    python3 projects/tgc/gui/app.py
"""

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PyQt5.QtCore import Qt, QThread, QTimer, pyqtSignal
from PyQt5.QtGui import QFont, QFontDatabase
from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSlider,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR       = Path(__file__).parent.resolve()              # …/projects/tgc/gui/
TGC_DIR          = (SCRIPT_DIR / "..").resolve()                # …/projects/tgc/
BINARY           = TGC_DIR / "build" / "tgc_sim"
GARFIELD_INSTALL = (TGC_DIR / "../../local/garfield").resolve() # …/local/garfield/


# ---------------------------------------------------------------------------
# Gas filename derivation
# ---------------------------------------------------------------------------

def derive_gas_filename(gas: dict) -> str:
    """Return a deterministic .gas filename from gas config parameters."""
    T   = round(gas.get("temperature_K", 293.15))
    P   = round(gas.get("pressure_Torr", 760.0))
    Ee  = round(gas.get("max_electron_energy_eV", 2000.0))
    Ef  = round(gas.get("e_field_max_vcm", 300000.0) / 1000)
    n   = gas.get("n_field_points", 20)
    c   = gas.get("n_magboltz_collisions", 10)
    pen = "pen" if gas.get("enable_penning", True) else "nopen"
    return f"ar70_co2_30_T{T}_P{P}_Ee{Ee}_Ef{Ef}k_n{n}_c{c}_{pen}.gas"


# ---------------------------------------------------------------------------
# Background simulation runner
# ---------------------------------------------------------------------------

class SimRunner(QThread):
    """Runs tgc_sim in a background thread and emits stdout line-by-line."""

    log_line = pyqtSignal(str)   # one stdout line
    finished = pyqtSignal(str)   # emits the run output directory on success
    failed   = pyqtSignal(str)   # emits an error message on failure

    def __init__(self, config_dict: dict, out_dir: str, parent=None):
        super().__init__(parent)
        self._config   = config_dict
        self._out_dir  = out_dir
        self._proc: subprocess.Popen | None = None

    # ── public ──────────────────────────────────────────────────────────

    def stop(self):
        """Ask the subprocess to terminate (called from the main thread)."""
        if self._proc and self._proc.poll() is None:
            self._proc.terminate()

    # ── private ─────────────────────────────────────────────────────────

    def run(self):  # noqa: PLR0912 — runs in the worker thread
        # Write a temporary config file
        try:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".json", delete=False
            ) as tf:
                json.dump(self._config, tf, indent=2)
                tmp_cfg = tf.name
        except OSError as exc:
            self.failed.emit(f"Could not write temporary config: {exc}")
            return

        cmd = [str(BINARY), "--config", tmp_cfg, "--out", self._out_dir]

        env = os.environ.copy()
        env["GARFIELD_INSTALL"] = str(GARFIELD_INSTALL)
        env["HEED_DATABASE"] = str(GARFIELD_INSTALL / "share" / "Heed" / "database")

        try:
            self._proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                cwd=str(TGC_DIR),   # relative gas-file paths resolve from here
                env=env,
            )
            for line in self._proc.stdout:
                self.log_line.emit(line.rstrip())
            self._proc.wait()
            ret = self._proc.returncode
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))
            return
        finally:
            try:
                os.unlink(tmp_cfg)
            except OSError:
                pass

        if ret != 0:
            self.failed.emit(f"Binary exited with code {ret}")
            return

        # Locate the sub-directory the binary created: <out_dir>/V<V>V__n<n>/
        out_path = Path(self._out_dir)
        subdirs  = sorted(out_path.glob("V*__n*"),
                          key=lambda p: p.stat().st_mtime)
        run_dir  = str(subdirs[-1]) if subdirs else self._out_dir
        self.finished.emit(run_dir)


# ---------------------------------------------------------------------------
# Config panel (left side)
# ---------------------------------------------------------------------------

class ConfigPanel(QScrollArea):
    """Scrollable panel with one QGroupBox per config section."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWidgetResizable(True)
        self.setMinimumWidth(330)

        container = QWidget()
        root_layout = QVBoxLayout(container)
        root_layout.setSpacing(8)

        # ── Geometry ──────────────────────────────────────────────────────
        geo_box  = QGroupBox("Geometry")
        geo_form = QFormLayout(geo_box)

        self.wire_pitch = self._dspin(0.05, 5.0, 0.01, 3, 0.18)
        self.wire_diam  = self._dspin(10.0, 500.0, 1.0, 1, 50.0)
        self.gap_cm     = self._dspin(0.05, 2.0, 0.01, 3, 0.14)
        self.n_wires    = self._spin(2, 100, 10)
        self.wire_volts = self._dspin(100.0, 5000.0, 50.0, 1, 1900.0)

        geo_form.addRow("Wire pitch [cm]",    self.wire_pitch)
        geo_form.addRow("Wire diameter [μm]", self.wire_diam)
        geo_form.addRow("Gap [cm]",           self.gap_cm)
        geo_form.addRow("N wires",            self.n_wires)
        geo_form.addRow("Wire voltage [V]",   self.wire_volts)
        root_layout.addWidget(geo_box)

        # ── Readout ───────────────────────────────────────────────────────
        ro_box  = QGroupBox("Readout")
        ro_form = QFormLayout(ro_box)

        self.readout_type = QComboBox()
        self.readout_type.addItems(["Conductive", "Resistive"])

        self.insulator_material = QComboBox()
        self.insulator_material.addItems(["Kapton", "FR4"])

        self.insulator_thickness = self._dspin(1.0, 10000.0, 10.0, 1, 100.0)
        self.insulator_thickness.setToolTip("Insulating substrate thickness [μm]")

        self.surface_resistivity = self._dspin(1.0, 1e7, 100.0, 1, 500.0)
        self.surface_resistivity.setToolTip("Resistive layer surface resistivity [kΩ/sq]")

        ro_form.addRow("Type",                   self.readout_type)
        ro_form.addRow("Insulator material",     self.insulator_material)
        ro_form.addRow("Thickness [μm]",         self.insulator_thickness)
        ro_form.addRow("Resistivity [kΩ/sq]",    self.surface_resistivity)
        root_layout.addWidget(ro_box)

        self._update_readout_widgets()
        self.readout_type.currentIndexChanged.connect(self._update_readout_widgets)

        # ── Source ────────────────────────────────────────────────────────
        src_box  = QGroupBox("Source")
        src_form = QFormLayout(src_box)

        self.energy_kev = self._dspin(0.1, 100.0, 0.1, 2, 5.9)
        self.distances  = QLineEdit("0.2,0.5,0.9,1.2")
        self.distances.setToolTip("Comma-separated source y-distances from wire plane [mm]")

        self.x_random = QCheckBox("Random (uniform over wire span)")
        self.x_random.setChecked(True)
        self.x_pos = self._dspin(-10.0, 10.0, 0.01, 3, 0.0)
        self.x_pos.setEnabled(False)
        self.x_pos.setToolTip("Fixed photon x-position [cm]")
        self.x_random.toggled.connect(lambda on: self.x_pos.setEnabled(not on))

        src_form.addRow("Energy [keV]",   self.energy_kev)
        src_form.addRow("Distances [mm]", self.distances)
        src_form.addRow("X position",     self.x_random)
        src_form.addRow("  fixed x [cm]", self.x_pos)
        root_layout.addWidget(src_box)

        # ── Gas ───────────────────────────────────────────────────────────
        gas_box  = QGroupBox("Gas")
        gas_form = QFormLayout(gas_box)

        self.temperature = self._dspin(200.0, 500.0, 1.0, 2, 293.15)
        self.pressure    = self._dspin(100.0, 3000.0, 10.0, 1, 760.0)

        self.penning = QCheckBox()
        self.penning.setChecked(True)
        self.ncoll = self._spin(1, 100, 10)
        self.ncoll.setToolTip("Magboltz collision cycles per field point (higher = more accurate)")
        self.w_value = self._dspin(10.0, 100.0, 0.5, 1, 26.0)
        self.w_value.setToolTip("Effective ionisation energy W [eV per ion pair] for primary electron count")

        self.max_electron_energy = self._dspin(100.0, 100_000.0, 100.0, 0, 2000.0)
        self.max_electron_energy.setToolTip(
            "Upper electron energy for Magboltz cross-section table [eV].\n"
            "Must exceed peak electron energy near the wire (~500–1000 eV)."
        )
        self.n_field_pts = self._spin(5, 500, 20)
        self.n_field_pts.setToolTip(
            "Number of log-spaced E-field points for the Magboltz transport table.\n"
            "More points → smoother interpolation; fewer → faster gas generation."
        )
        self.e_field_max = self._dspin(10_000.0, 1_000_000.0, 10_000.0, 0, 300_000.0)
        self.e_field_max.setToolTip(
            "Maximum E-field in the Magboltz table [V/cm].\n"
            "Must exceed the peak near-wire field (~200–400 kV/cm at 1900 V)."
        )

        self.gas_file_label = QLabel()
        self.gas_file_label.setWordWrap(True)
        self.gas_file_label.setStyleSheet("font-size: 10px;")

        gas_form.addRow("Temperature [K]",     self.temperature)
        gas_form.addRow("Pressure [Torr]",     self.pressure)
        gas_form.addRow("Penning transfer",    self.penning)
        gas_form.addRow("Magboltz ncoll",      self.ncoll)
        gas_form.addRow("W-value [eV]",        self.w_value)
        gas_form.addRow("Max e⁻ energy [eV]",  self.max_electron_energy)
        gas_form.addRow("Field points",        self.n_field_pts)
        gas_form.addRow("E-field max [V/cm]",  self.e_field_max)
        gas_form.addRow("Gas file (auto)",     self.gas_file_label)
        root_layout.addWidget(gas_box)

        # Update gas file label whenever a relevant parameter changes
        self.temperature.valueChanged.connect(self._update_gas_file_label)
        self.pressure.valueChanged.connect(self._update_gas_file_label)
        self.penning.toggled.connect(self._update_gas_file_label)
        self.ncoll.valueChanged.connect(self._update_gas_file_label)
        self.max_electron_energy.valueChanged.connect(self._update_gas_file_label)
        self.n_field_pts.valueChanged.connect(self._update_gas_file_label)
        self.e_field_max.valueChanged.connect(self._update_gas_file_label)

        # ── Simulation ────────────────────────────────────────────────────
        sim_box  = QGroupBox("Simulation")
        sim_form = QFormLayout(sim_box)

        self.n_events    = self._spin(1, 100000, 1000)
        self.max_aval    = self._spin(1000, 10000000, 500000)
        self.time_window = self._dspin(10.0, 20000.0, 10.0, 1, 300.0)
        self.time_step   = self._dspin(0.1, 10.0, 0.1, 2, 0.5)
        self.enable_ion_drift = QCheckBox()
        self.enable_ion_drift.setChecked(True)
        self.enable_ion_drift.setToolTip(
            "Drift positive ions after each avalanche (DriftLineRKF).\n"
            "Disabling skips ion signal computation and greatly speeds up runs."
        )
        self.store_drift_lines = QCheckBox()
        self.store_drift_lines.setChecked(False)
        self.store_drift_lines.setToolTip(
            "Store intermediate drift-line steps for the 3D track view.\n"
            "Off: primary electron shown as straight start→end line.\n"
            "On: full curved trajectory, but adds ~15 % CPU time per event."
        )

        sim_form.addRow("Events",             self.n_events)
        sim_form.addRow("Max avalanche size", self.max_aval)
        sim_form.addRow("Time window [ns]",   self.time_window)
        sim_form.addRow("Time step [ns]",     self.time_step)
        sim_form.addRow("Ion transport",      self.enable_ion_drift)
        sim_form.addRow("Store drift lines",  self.store_drift_lines)
        root_layout.addWidget(sim_box)

        # ── Output ────────────────────────────────────────────────────────
        out_box  = QGroupBox("Output")
        out_form = QFormLayout(out_box)

        out_row = QWidget()
        out_h   = QHBoxLayout(out_row)
        out_h.setContentsMargins(0, 0, 0, 0)
        self.out_dir = QLineEdit("results")
        self.out_dir.setToolTip("Output base directory (relative paths resolve from projects/tgc/)")
        btn_out = QPushButton("…")
        btn_out.setFixedWidth(28)
        btn_out.clicked.connect(self._browse_out_dir)
        out_h.addWidget(self.out_dir)
        out_h.addWidget(btn_out)

        out_form.addRow("Directory", out_row)
        root_layout.addWidget(out_box)

        root_layout.addStretch()
        self.setWidget(container)

        self._update_gas_file_label()

    # ── widget factories ─────────────────────────────────────────────────

    @staticmethod
    def _dspin(lo: float, hi: float, step: float, dec: int, val: float) -> QDoubleSpinBox:
        w = QDoubleSpinBox()
        w.setRange(lo, hi)
        w.setSingleStep(step)
        w.setDecimals(dec)
        w.setValue(val)
        return w

    @staticmethod
    def _spin(lo: int, hi: int, val: int) -> QSpinBox:
        w = QSpinBox()
        w.setRange(lo, hi)
        w.setValue(val)
        return w

    # ── gas file label ───────────────────────────────────────────────────

    def _update_gas_file_label(self):
        gas = {
            "temperature_K":          self.temperature.value(),
            "pressure_Torr":          self.pressure.value(),
            "enable_penning":         self.penning.isChecked(),
            "n_magboltz_collisions":  self.ncoll.value(),
            "max_electron_energy_eV": self.max_electron_energy.value(),
            "n_field_points":         self.n_field_pts.value(),
            "e_field_max_vcm":        self.e_field_max.value(),
        }
        name = derive_gas_filename(gas)
        exists = (TGC_DIR / name).exists()
        status = "✓ exists" if exists else "will be generated"
        self.gas_file_label.setText(f"{name}\n[{status}]")

    def _update_readout_widgets(self):
        resistive = self.readout_type.currentText() == "Resistive"
        self.insulator_material.setEnabled(resistive)
        self.insulator_thickness.setEnabled(resistive)
        self.surface_resistivity.setEnabled(resistive)

    # ── file dialogs ─────────────────────────────────────────────────────

    def _browse_out_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Select output directory")
        if path:
            self.out_dir.setText(path)

    # ── serialisation ─────────────────────────────────────────────────────

    def to_config_dict(self) -> dict:
        """Assemble widget values into a config dict suitable for JSON dump."""
        raw = self.distances.text().strip()
        try:
            dists = [float(x.strip()) for x in raw.split(",") if x.strip()]
        except ValueError:
            dists = [0.7]

        x_pos = None if self.x_random.isChecked() else self.x_pos.value()

        ro_type = self.readout_type.currentText().lower()
        ins_mat = self.insulator_material.currentText().lower()

        return {
            "geometry": {
                "wire_pitch_cm":    self.wire_pitch.value(),
                "wire_diameter_um": self.wire_diam.value(),
                "gap_cm":           self.gap_cm.value(),
                "n_wires":          self.n_wires.value(),
                "wire_voltage_V":   self.wire_volts.value(),
            },
            "readout": {
                "type":                       ro_type,
                "insulator_material":         ins_mat,
                "insulator_thickness_um":     self.insulator_thickness.value(),
                "surface_resistivity_ohm_sq": self.surface_resistivity.value() * 1000.0,
            },
            "source": {
                "energy_keV":          self.energy_kev.value(),
                "source_distances_mm": dists,
                "x_position_cm":       x_pos,
            },
            "gas": {
                "temperature_K":          self.temperature.value(),
                "pressure_Torr":          self.pressure.value(),
                "enable_penning":         self.penning.isChecked(),
                "n_magboltz_collisions":  self.ncoll.value(),
                "w_value_eV":             self.w_value.value(),
                "max_electron_energy_eV": self.max_electron_energy.value(),
                "n_field_points":         self.n_field_pts.value(),
                "e_field_max_vcm":        self.e_field_max.value(),
            },
            "simulation": {
                "n_events":           self.n_events.value(),
                "max_avalanche_size": self.max_aval.value(),
                "time_window_ns":     self.time_window.value(),
                "time_step_ns":       self.time_step.value(),
                "enable_ion_drift":   self.enable_ion_drift.isChecked(),
                "store_drift_lines":  self.store_drift_lines.isChecked(),
            },
        }

    def load_from_dict(self, d: dict):
        """Populate all widgets from a config dict (e.g. loaded from JSON)."""
        g = d.get("geometry", {})
        self.wire_pitch.setValue(g.get("wire_pitch_cm", 0.18))
        self.wire_diam.setValue( g.get("wire_diameter_um", 50.0))
        self.gap_cm.setValue(    g.get("gap_cm", 0.14))
        self.n_wires.setValue(   g.get("n_wires", 10))
        self.wire_volts.setValue(g.get("wire_voltage_V", 1900.0))

        ro = d.get("readout", {})
        ro_type = ro.get("type", "conductive")
        self.readout_type.setCurrentIndex(0 if ro_type == "conductive" else 1)
        ins_mat = ro.get("insulator_material", "kapton")
        self.insulator_material.setCurrentIndex(0 if ins_mat == "kapton" else 1)
        self.insulator_thickness.setValue(ro.get("insulator_thickness_um", 100.0))
        self.surface_resistivity.setValue(ro.get("surface_resistivity_ohm_sq", 500000.0) / 1000.0)

        s = d.get("source", {})
        self.energy_kev.setValue(s.get("energy_keV", 5.9))
        dists = s.get("source_distances_mm", [0.2, 0.5, 0.9, 1.2])
        self.distances.setText(",".join(str(v) for v in dists))
        x_pos = s.get("x_position_cm", None)
        if x_pos is None:
            self.x_random.setChecked(True)
        else:
            self.x_random.setChecked(False)
            self.x_pos.setValue(float(x_pos))

        gas = d.get("gas", {})
        self.temperature.setValue(gas.get("temperature_K", 293.15))
        self.pressure.setValue(   gas.get("pressure_Torr", 760.0))
        self.penning.setChecked(  gas.get("enable_penning", True))
        self.ncoll.setValue(      gas.get("n_magboltz_collisions", 10))
        self.w_value.setValue(    gas.get("w_value_eV", 26.0))
        self.max_electron_energy.setValue(gas.get("max_electron_energy_eV", 2000.0))
        self.n_field_pts.setValue(        gas.get("n_field_points", 20))
        self.e_field_max.setValue(        gas.get("e_field_max_vcm", 300_000.0))

        sim = d.get("simulation", {})
        self.n_events.setValue(        sim.get("n_events", 1000))
        self.max_aval.setValue(        sim.get("max_avalanche_size", 500000))
        self.time_window.setValue(     sim.get("time_window_ns", 300.0))
        self.time_step.setValue(       sim.get("time_step_ns", 0.5))
        self.enable_ion_drift.setChecked(sim.get("enable_ion_drift", True))
        self.store_drift_lines.setChecked(sim.get("store_drift_lines", False))


# ---------------------------------------------------------------------------
# Matplotlib canvas wrapper
# ---------------------------------------------------------------------------

class MplCanvas(FigureCanvasQTAgg):
    """A matplotlib Figure embedded in a Qt widget."""

    def __init__(self, nrows: int = 1, ncols: int = 1,
                 figsize: tuple | None = None, parent=None):
        self.figure = Figure(figsize=figsize or (6, 4))
        self.axes: list = [
            self.figure.add_subplot(nrows, ncols, i + 1)
            for i in range(nrows * ncols)
        ]
        super().__init__(self.figure)


class Mpl3DCanvas(FigureCanvasQTAgg):
    """A single matplotlib 3D axes embedded in a Qt widget."""

    def __init__(self, figsize: tuple = (8, 5), parent=None):
        self.figure = Figure(figsize=figsize)
        self.ax = self.figure.add_subplot(111, projection="3d")
        super().__init__(self.figure)
        self._user_dist = None  # persists user zoom across redraws
        self.mpl_connect("scroll_event", self._on_scroll)

    def _on_scroll(self, event):
        if event.inaxes is not self.ax:
            return
        self._user_dist = self.ax.dist * (0.9 if event.step > 0 else 1.1)
        self.ax.dist = self._user_dist
        self.draw_idle()


# ---------------------------------------------------------------------------
# Results panel (right side)
# ---------------------------------------------------------------------------

class ResultsPanel(QTabWidget):
    """Four-tab panel: Log | Summary | Plots | Waveforms."""

    def __init__(self, parent=None):
        super().__init__(parent)

        # ── Log tab ───────────────────────────────────────────────────────
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        mono = QFontDatabase.systemFont(QFontDatabase.FixedFont)
        mono.setPointSize(11)
        self.log.setFont(mono)
        self.addTab(self.log, "Log")

        # ── Summary tab ───────────────────────────────────────────────────
        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.addTab(self.table, "Summary")

        # ── Plots tab: 2×3 matplotlib figure ─────────────────────────────
        self.plots_canvas = MplCanvas(nrows=2, ncols=3, figsize=(11, 5))
        self.addTab(self.plots_canvas, "Plots")

        # ── Waveforms tab: ROOT TCanvas browser ──────────────────────────
        self._waveform_data: dict = {}
        self._root_canvas  = None    # ROOT TCanvas (kept alive between events)
        self._root_objects: list = []  # TGraph/TLegend objects (prevent Python GC)
        self._charge_canvas  = None   # ROOT TCanvas for charge integrals
        self._charge_objects: list = []
        self._track_data:   dict = {}  # label → dict of numpy object arrays
        self._track_geom:   dict | None = None

        wave_widget = QWidget()
        wave_layout = QVBoxLayout(wave_widget)
        wave_layout.setContentsMargins(8, 6, 8, 6)
        wave_layout.setSpacing(6)

        # — selector row —
        sel_row = QWidget()
        sel_h   = QHBoxLayout(sel_row)
        sel_h.setContentsMargins(0, 0, 0, 0)
        sel_h.addWidget(QLabel("Distance:"))
        self.wave_dist_combo = QComboBox()
        sel_h.addWidget(self.wave_dist_combo)
        sel_h.addSpacing(16)
        sel_h.addWidget(QLabel("Event:"))
        self.wave_event_slider = QSlider(Qt.Horizontal)
        self.wave_event_slider.setMinimum(0)
        self.wave_event_slider.setMaximum(0)
        self.wave_event_slider.setSingleStep(1)
        sel_h.addWidget(self.wave_event_slider)
        self.wave_event_label = QLabel("— / —")
        self.wave_event_label.setMinimumWidth(55)
        sel_h.addWidget(self.wave_event_label)
        wave_layout.addWidget(sel_row)

        # — per-event info —
        info_box  = QGroupBox("Current event")
        info_form = QFormLayout(info_box)
        info_form.setVerticalSpacing(2)
        self.wave_qa_lbl    = QLabel("—")
        self.wave_qc_lbl    = QLabel("—")
        self.wave_ratio_lbl = QLabel("—")
        info_form.addRow("Q anode [fC]:",   self.wave_qa_lbl)
        info_form.addRow("Q cathode [fC]:", self.wave_qc_lbl)
        info_form.addRow("Ratio:",          self.wave_ratio_lbl)
        wave_layout.addWidget(info_box)

        # — hint —
        hint = QLabel(
            "ROOT canvas opens automatically when results are loaded.\n"
            "Right-click inside the ROOT window to zoom, change axes, or save."
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color: grey; font-size: 11px;")
        wave_layout.addWidget(hint)
        wave_layout.addStretch()

        self.wave_dist_combo.currentIndexChanged.connect(self._on_wave_dist_changed)
        self.wave_event_slider.valueChanged.connect(self._update_waveform_plot)

        self.addTab(wave_widget, "Waveforms")

        # ── Charges tab: cumulative integral of waveforms ────────────────────
        charge_widget  = QWidget()
        charge_layout  = QVBoxLayout(charge_widget)
        charge_layout.setContentsMargins(8, 6, 8, 6)
        charge_layout.setSpacing(6)

        # — selector row —
        ch_sel_row = QWidget()
        ch_sel_h   = QHBoxLayout(ch_sel_row)
        ch_sel_h.setContentsMargins(0, 0, 0, 0)
        ch_sel_h.addWidget(QLabel("Distance:"))
        self.charge_dist_combo = QComboBox()
        ch_sel_h.addWidget(self.charge_dist_combo)
        ch_sel_h.addSpacing(16)
        ch_sel_h.addWidget(QLabel("Event:"))
        self.charge_event_slider = QSlider(Qt.Horizontal)
        self.charge_event_slider.setMinimum(0)
        self.charge_event_slider.setMaximum(0)
        self.charge_event_slider.setSingleStep(1)
        ch_sel_h.addWidget(self.charge_event_slider)
        self.charge_event_label = QLabel("— / —")
        self.charge_event_label.setMinimumWidth(55)
        ch_sel_h.addWidget(self.charge_event_label)
        charge_layout.addWidget(ch_sel_row)

        # — hint —
        ch_hint = QLabel(
            "ROOT canvas opens automatically when results are loaded.\n"
            "Right-click inside the ROOT window to zoom, change axes, or save."
        )
        ch_hint.setWordWrap(True)
        ch_hint.setStyleSheet("color: grey; font-size: 11px;")
        charge_layout.addWidget(ch_hint)
        charge_layout.addStretch()

        self.charge_dist_combo.currentIndexChanged.connect(self._on_charge_dist_changed)
        self.charge_event_slider.valueChanged.connect(self._update_charge_plot)

        self.addTab(charge_widget, "Charges")

        # ── 3D Tracks tab ─────────────────────────────────────────────────────
        tracks_widget = QWidget()
        tracks_layout = QVBoxLayout(tracks_widget)
        tracks_layout.setContentsMargins(8, 6, 8, 6)
        tracks_layout.setSpacing(6)

        # — selector row —
        trk_sel_row = QWidget()
        trk_sel_h   = QHBoxLayout(trk_sel_row)
        trk_sel_h.setContentsMargins(0, 0, 0, 0)
        trk_sel_h.addWidget(QLabel("Distance:"))
        self.trk_dist_combo = QComboBox()
        trk_sel_h.addWidget(self.trk_dist_combo)
        trk_sel_h.addSpacing(16)
        trk_sel_h.addWidget(QLabel("Event:"))
        self.trk_event_slider = QSlider(Qt.Horizontal)
        self.trk_event_slider.setMinimum(0)
        self.trk_event_slider.setMaximum(0)
        self.trk_event_slider.setSingleStep(1)
        trk_sel_h.addWidget(self.trk_event_slider)
        self.trk_event_label = QLabel("— / —")
        self.trk_event_label.setMinimumWidth(55)
        trk_sel_h.addWidget(self.trk_event_label)
        rel_hint = QLabel("(release slider to update)")
        rel_hint.setStyleSheet("color: grey; font-size: 10px;")
        trk_sel_h.addWidget(rel_hint)
        trk_sel_h.addStretch()
        tracks_layout.addWidget(trk_sel_row)

        self.tracks_canvas = Mpl3DCanvas(figsize=(8, 5))
        tracks_layout.addWidget(self.tracks_canvas, stretch=1)

        self.trk_dist_combo.currentIndexChanged.connect(self._on_trk_dist_changed)
        self.trk_event_slider.sliderReleased.connect(self._update_track_plot)

        self.addTab(tracks_widget, "3D Tracks")

        # — timer to keep ROOT canvas responsive —
        self._root_timer = QTimer(self)
        self._root_timer.setInterval(100)
        self._root_timer.timeout.connect(self._process_root_events)

    # ── Log ───────────────────────────────────────────────────────────────

    def append_log(self, line: str):
        self.log.appendPlainText(line)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())

    def clear_log(self):
        self.log.clear()

    # ── Summary table ─────────────────────────────────────────────────────

    def populate_table(self, csv_path: str):
        try:
            df = pd.read_csv(csv_path)
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Could not read CSV: {exc}")
            return

        df.columns = [c.strip() for c in df.columns]
        self.table.clear()
        self.table.setRowCount(len(df))
        self.table.setColumnCount(len(df.columns))
        self.table.setHorizontalHeaderLabels(list(df.columns))

        for r_idx, row in df.iterrows():
            for c_idx, val in enumerate(row):
                text = f"{val:.4g}" if isinstance(val, float) else str(val)
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(r_idx, c_idx, item)

        self.table.resizeColumnsToContents()

    # ── Plots ─────────────────────────────────────────────────────────────

    def draw_plots(self, csv_path: str):
        try:
            df = pd.read_csv(csv_path)
        except Exception:  # noqa: BLE001
            return

        df.columns = [c.strip() for c in df.columns]
        axes = self.plots_canvas.axes
        for ax in axes:
            ax.cla()

        x = df["source_distance_mm"].to_numpy()

        def _errplot(ax, y_col: str, yerr_col: str | None,
                     ylabel: str, title: str):
            if y_col not in df.columns:
                ax.set_title(title + " (no data)")
                return
            y  = df[y_col].to_numpy()
            ye = df[yerr_col].to_numpy() if (yerr_col and yerr_col in df.columns) else None
            ax.errorbar(x, y, yerr=ye, fmt="o-", capsize=4, lw=1.5, ms=5)
            ax.set_xlabel("Source distance [mm]")
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.grid(True, alpha=0.3)

        _errplot(axes[0], "mean_anode_charge_fC",       "sem_anode_charge_fC",
                 "⟨Q_anode⟩ [fC]",            "Anode charge vs distance")
        _errplot(axes[1], "mean_cathode_charge_fC",     "sem_cathode_charge_fC",
                 "⟨Q_cathode⟩ [fC]",          "Cathode (readout) vs distance")
        _errplot(axes[2], "mean_cathode_top_charge_fC", "sem_cathode_top_charge_fC",
                 "⟨Q_cathode_top⟩ [fC]",      "Cathode-top vs distance")
        _errplot(axes[3], "mean_charge_ratio",          "sem_charge_ratio",
                 "Q_cathode / Q_anode",        "Charge ratio vs distance")
        _errplot(axes[4], "mean_avalanche_size",        None,
                 "⟨Avalanche size⟩",           "Avalanche size vs distance")
        if len(axes) > 5:
            axes[5].axis("off")   # unused 6th cell

        self.plots_canvas.figure.tight_layout()
        self.plots_canvas.draw()

    # ── Waveforms ─────────────────────────────────────────────────────────

    def _process_root_events(self):
        """Keep the ROOT TCanvas window responsive while Qt runs."""
        try:
            import ROOT  # noqa: PLC0415
            ROOT.gSystem.ProcessEvents()
        except Exception:  # noqa: BLE001
            pass

    def load_waveform_data(self, root_path: str):
        """Load per-event TTree waveforms (and mean TProfiles) from the ROOT file."""
        try:
            import uproot  # noqa: PLC0415
        except ImportError:
            self.append_log("[GUI] uproot not available — waveform tab will be empty")
            return

        if not Path(root_path).exists():
            self.append_log(f"[GUI] ROOT file not found: {root_path}")
            return

        self._waveform_data.clear()
        self.wave_dist_combo.blockSignals(True)
        self.wave_dist_combo.clear()
        self.charge_dist_combo.blockSignals(True)
        self.charge_dist_combo.clear()

        try:
            with uproot.open(root_path) as f:
                dist_keys = sorted({
                    k.split("/")[0] for k in f.keys(cycle=False)
                    if k.split("/")[0].startswith("dist_")
                })
                for key in dist_keys:
                    label = key.removeprefix("dist_").replace("p", ".").replace("mm", " mm")
                    try:
                        pa     = f[f"{key}/p_anode_signal"]
                        pc     = f[f"{key}/p_cathode_signal"]
                        times  = pa.axis(0).centers()
                        mean_a = pa.values()
                        mean_c = pc.values()

                        # std::vector<float> branches → object array of 1D arrays;
                        # np.stack() converts to a proper (n_evt, nBins) 2D array.
                        tree    = f[f"{key}/t_signals"]
                        anode   = np.stack(tree["anode"].array(library="np"))
                        cathode = np.stack(tree["cathode"].array(library="np"))

                        self._waveform_data[label] = {
                            "times":   times,
                            "anode":   anode,
                            "cathode": cathode,
                            "mean_a":  mean_a,
                            "mean_c":  mean_c,
                        }
                        self.wave_dist_combo.addItem(label)
                        self.charge_dist_combo.addItem(label)
                    except Exception as exc:  # noqa: BLE001
                        self.append_log(f"[GUI] Waveforms: could not read {key}: {exc}")
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Could not open ROOT file: {exc}")

        self.wave_dist_combo.blockSignals(False)
        self.charge_dist_combo.blockSignals(False)
        if self.wave_dist_combo.count():
            self._on_wave_dist_changed(0)
            self._root_timer.start()
        if self.charge_dist_combo.count():
            self._on_charge_dist_changed(0)

    def _on_wave_dist_changed(self, index: int):
        label = self.wave_dist_combo.currentText()
        data  = self._waveform_data.get(label)
        if data is None:
            return
        n = len(data["anode"])
        self.wave_event_slider.blockSignals(True)
        self.wave_event_slider.setMaximum(max(0, n - 1))
        self.wave_event_slider.setValue(0)
        self.wave_event_slider.blockSignals(False)
        self.wave_event_label.setText(f"1 / {n}")
        self._update_waveform_plot()

    def _update_waveform_plot(self):
        """Draw the selected event in a ROOT TCanvas (anode top, cathode bottom)."""
        label = self.wave_dist_combo.currentText()
        data  = self._waveform_data.get(label)
        if data is None:
            return

        evt_idx = self.wave_event_slider.value()
        n       = len(data["anode"])
        self.wave_event_label.setText(f"{evt_idx + 1} / {n}")

        times   = data["times"].astype("f8")
        anode   = data["anode"][evt_idx].astype("f8")
        cathode = data["cathode"][evt_idx].astype("f8")
        mean_a  = data["mean_a"].astype("f8")
        mean_c  = data["mean_c"].astype("f8")
        nbins   = len(times)
        dt      = float(times[1] - times[0]) if nbins > 1 else 1.0

        qa    = float(-np.sum(anode)   * dt)
        qc    = float( np.sum(cathode) * dt)
        ratio = qc / qa if qa != 0 else float("nan")

        self.wave_qa_lbl.setText(f"{qa:.4g}")
        self.wave_qc_lbl.setText(f"{qc:.4g}")
        self.wave_ratio_lbl.setText(f"{ratio:.4g}")

        try:
            import ROOT  # noqa: PLC0415
            ROOT.gROOT.SetBatch(False)

            if self._root_canvas is None:
                self._root_canvas = ROOT.TCanvas(
                    "tgc_waveforms", "TGC Waveforms", 950, 700
                )
                self._root_canvas.SetWindowSize(950, 700)

            self._root_canvas.Clear()
            self._root_canvas.Divide(1, 2)
            self._root_objects.clear()   # release previous objects

            def _draw_pad(pad_idx: int, y_evt: np.ndarray, y_mean: np.ndarray,
                          signal_name: str, line_color: int):
                self._root_canvas.cd(pad_idx)
                ROOT.gPad.SetGrid()
                ga = ROOT.TGraph(nbins, times, y_evt)
                ga.SetTitle(
                    f"{signal_name} - {label}, event {evt_idx + 1};"
                    f"Time [ns];i [fC/ns]"
                )
                ga.SetLineColor(line_color)
                ga.SetLineWidth(2)
                gm = ROOT.TGraph(nbins, times, y_mean)
                gm.SetLineColor(ROOT.kGray + 1)
                gm.SetLineWidth(1)
                gm.SetLineStyle(2)    # dashed mean
                ga.Draw("AL")
                gm.Draw("L same")
                leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
                leg.AddEntry(ga, "this event", "L")
                leg.AddEntry(gm, "mean",       "L")
                leg.SetBorderSize(0)
                leg.Draw()
                self._root_objects.extend([ga, gm, leg])

            _draw_pad(1, anode,   mean_a, "Anode",   ROOT.kBlue + 1)
            _draw_pad(2, cathode, mean_c, "Cathode", ROOT.kRed  + 1)

            self._root_canvas.Update()

        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] ROOT canvas error: {exc}")

    # ── Charges (cumulative integral) ─────────────────────────────────────────

    def _on_charge_dist_changed(self, index: int):
        label = self.charge_dist_combo.currentText()
        data  = self._waveform_data.get(label)
        if data is None:
            return
        n = len(data["anode"])
        self.charge_event_slider.blockSignals(True)
        self.charge_event_slider.setMaximum(max(0, n - 1))
        self.charge_event_slider.setValue(0)
        self.charge_event_slider.blockSignals(False)
        self.charge_event_label.setText(f"1 / {n}")
        self._update_charge_plot()

    def _update_charge_plot(self):
        """Draw cumulative charge integrals Q(t) in a ROOT TCanvas (anode top, cathode bottom)."""
        label = self.charge_dist_combo.currentText()
        data  = self._waveform_data.get(label)
        if data is None:
            return

        evt_idx = self.charge_event_slider.value()
        n       = len(data["anode"])
        self.charge_event_label.setText(f"{evt_idx + 1} / {n}")

        times   = data["times"].astype("f8")
        anode   = data["anode"][evt_idx].astype("f8")
        cathode = data["cathode"][evt_idx].astype("f8")
        mean_a  = data["mean_a"].astype("f8")
        mean_c  = data["mean_c"].astype("f8")
        nbins   = len(times)
        dt      = float(times[1] - times[0]) if nbins > 1 else 1.0

        # Cumulative integrals [fC]; anode signal is negative by convention → negate.
        anode_int   = -np.cumsum(anode)  * dt
        cathode_int =  np.cumsum(cathode) * dt
        mean_a_int  = -np.cumsum(mean_a) * dt
        mean_c_int  =  np.cumsum(mean_c) * dt

        try:
            import ROOT  # noqa: PLC0415
            ROOT.gROOT.SetBatch(False)

            if self._charge_canvas is None:
                self._charge_canvas = ROOT.TCanvas(
                    "tgc_charges", "TGC Charges", 950, 700
                )
                self._charge_canvas.SetWindowSize(950, 700)

            self._charge_canvas.Clear()
            self._charge_canvas.Divide(1, 2)
            self._charge_objects.clear()

            def _draw_pad(pad_idx: int, y_evt: np.ndarray, y_mean: np.ndarray,
                          signal_name: str, line_color: int):
                self._charge_canvas.cd(pad_idx)
                ROOT.gPad.SetGrid()
                ga = ROOT.TGraph(nbins, times, y_evt)
                ga.SetTitle(
                    f"{signal_name} charge - {label}, event {evt_idx + 1};"
                    f"Time [ns];Q [fC]"
                )
                ga.SetLineColor(line_color)
                ga.SetLineWidth(2)
                gm = ROOT.TGraph(nbins, times, y_mean)
                gm.SetLineColor(ROOT.kGray + 1)
                gm.SetLineWidth(1)
                gm.SetLineStyle(2)
                ga.Draw("AL")
                gm.Draw("L same")
                leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
                leg.AddEntry(ga, "this event", "L")
                leg.AddEntry(gm, "mean",       "L")
                leg.SetBorderSize(0)
                leg.Draw()
                self._charge_objects.extend([ga, gm, leg])

            _draw_pad(1, anode_int,   mean_a_int, "Anode",   ROOT.kBlue + 1)
            _draw_pad(2, cathode_int, mean_c_int, "Cathode", ROOT.kRed  + 1)

            self._charge_canvas.Update()

        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] ROOT charge canvas error: {exc}")

    # ── 3D Tracks ──────────────────────────────────────────────────────────────

    def load_track_data(self, root_path: str, run_dir: str):
        """Load per-event track branches and geometry from the ROOT file."""
        try:
            import uproot  # noqa: PLC0415
        except ImportError:
            self.append_log("[GUI] uproot not available — 3D Tracks tab disabled")
            return

        self._track_data.clear()
        self._track_geom = None
        self.trk_dist_combo.blockSignals(True)
        self.trk_dist_combo.clear()

        # Geometry comes from run_config.json written by the binary
        cfg_path = Path(run_dir) / "run_config.json"
        if cfg_path.exists():
            try:
                with open(cfg_path) as f:
                    rc = json.load(f)
                g = rc.get("geometry", {})
                self._track_geom = {
                    "wire_pitch_cm": g.get("wire_pitch_cm", 0.18),
                    "n_wires":       int(g.get("n_wires", 10)),
                    "gap_cm":        g.get("gap_cm", 0.14),
                }
            except Exception as exc:  # noqa: BLE001
                self.append_log(f"[GUI] run_config.json read error: {exc}")

        if not Path(root_path).exists():
            self.trk_dist_combo.blockSignals(False)
            return

        try:
            with uproot.open(root_path) as f:
                dist_keys = sorted({
                    k.split("/")[0] for k in f.keys(cycle=False)
                    if k.split("/")[0].startswith("dist_")
                })
                for key in dist_keys:
                    label = key.removeprefix("dist_").replace("p", ".").replace("mm", " mm")
                    try:
                        tree = f[f"{key}/t_signals"]
                        if "primary_x" not in tree.keys():
                            continue   # pre-feature ROOT file — skip silently
                        self._track_data[label] = {
                            "primary_x": tree["primary_x"].array(library="np"),
                            "primary_y": tree["primary_y"].array(library="np"),
                            "primary_z": tree["primary_z"].array(library="np"),
                            "cloud_x":   tree["cloud_x"].array(library="np"),
                            "cloud_y":   tree["cloud_y"].array(library="np"),
                            "cloud_z":   tree["cloud_z"].array(library="np"),
                            "ion_x":     tree["ion_x"].array(library="np"),
                            "ion_y":     tree["ion_y"].array(library="np"),
                            "ion_z":     tree["ion_z"].array(library="np"),
                            "ion_npts":  tree["ion_npts"].array(library="np"),
                        }
                        self.trk_dist_combo.addItem(label)
                    except Exception as exc:  # noqa: BLE001
                        self.append_log(f"[GUI] Track load error for {key}: {exc}")
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Could not open ROOT file for tracks: {exc}")

        self.trk_dist_combo.blockSignals(False)
        if self.trk_dist_combo.count():
            self._on_trk_dist_changed(0)

    def _on_trk_dist_changed(self, index: int):
        label = self.trk_dist_combo.currentText()
        data  = self._track_data.get(label)
        if data is None:
            return
        n = len(data["primary_x"])
        self.trk_event_slider.blockSignals(True)
        self.trk_event_slider.setMaximum(max(0, n - 1))
        self.trk_event_slider.setValue(0)
        self.trk_event_slider.blockSignals(False)
        self.trk_event_label.setText(f"1 / {n}")
        self._update_track_plot()

    def _update_track_plot(self):
        """Render detector geometry and per-event tracks in the 3D canvas."""
        label = self.trk_dist_combo.currentText()
        data  = self._track_data.get(label)
        if data is None:
            return

        ev = self.trk_event_slider.value()
        n  = len(data["primary_x"])
        self.trk_event_label.setText(f"{ev + 1} / {n}")

        ax = self.tracks_canvas.ax
        ax.cla()

        # ── Geometry ──────────────────────────────────────────────────────────
        geom    = self._track_geom or {}
        pitch   = geom.get("wire_pitch_cm", 0.18)
        n_wires = int(geom.get("n_wires", 10))
        gap     = geom.get("gap_cm", 0.14)
        z_half  = 0.5   # sensor z extent [cm]
        x_half  = (n_wires - 1) / 2.0 * pitch + pitch

        # Cathode planes as translucent rectangles
        try:
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: PLC0415
            for y_plane, fc in [(-gap, "cyan"), (+gap, "lightyellow")]:
                verts = [[
                    [-x_half, y_plane, -z_half], [ x_half, y_plane, -z_half],
                    [ x_half, y_plane,  z_half], [-x_half, y_plane,  z_half],
                ]]
                poly = Poly3DCollection(verts, alpha=0.13,
                                        facecolor=fc, edgecolor="grey",
                                        linewidths=0.5)
                ax.add_collection3d(poly)
        except ImportError:
            pass   # very unlikely — mpl_toolkits ships with matplotlib

        # Anode wires (gold lines along z)
        for i in range(n_wires):
            xw = (i - (n_wires - 1) / 2.0) * pitch
            ax.plot([xw, xw], [0.0, 0.0], [-z_half, z_half],
                    color="gold", lw=1.2, zorder=3)

        # ── Per-event tracks ──────────────────────────────────────────────────
        # Primary electron drift line (blue)
        ax.plot(data["primary_x"][ev], data["primary_y"][ev], data["primary_z"][ev],
                color="royalblue", lw=1.8, label="Primary e⁻", zorder=5)

        # Avalanche cloud: start positions of secondary electron tracks (orange)
        cx = data["cloud_x"][ev]
        cy = data["cloud_y"][ev]
        cz = data["cloud_z"][ev]
        if len(cx) > 0:
            ax.scatter(cx, cy, cz, c="darkorange", s=4, alpha=0.5,
                       label="Avalanche cloud", zorder=4, depthshade=False)

        # Ion drift paths, colour-coded by destination cathode
        ix    = data["ion_x"][ev]
        iy    = data["ion_y"][ev]
        iz    = data["ion_z"][ev]
        inpts = data["ion_npts"][ev]
        if len(inpts) > 0:
            splits = np.concatenate([[0], np.cumsum(inpts)])
            tol    = 0.01 * gap
            for k in range(len(inpts)):
                xs = ix[splits[k]:splits[k + 1]]
                ys = iy[splits[k]:splits[k + 1]]
                zs = iz[splits[k]:splits[k + 1]]
                if len(ys) == 0:
                    continue
                y_end = float(ys[-1])
                if abs(y_end - (-gap)) <= tol:
                    color = "limegreen"   # → readout cathode
                elif abs(y_end - gap) <= tol:
                    color = "magenta"     # → non-readout cathode
                else:
                    color = "grey"        # absorbed / out of window
                ax.plot(xs, ys, zs, color=color, lw=0.8, alpha=0.7)

        # ── Axes / labels / legend ────────────────────────────────────────────
        ax.set_xlabel("x [cm]", fontsize=8, labelpad=1)
        ax.set_ylabel("y [cm]", fontsize=8, labelpad=1)
        ax.set_zlabel("z [cm]", fontsize=8, labelpad=1)
        ax.set_title(f"{label},  event {ev + 1}", fontsize=9)
        ax.set_xlim(-x_half,    x_half)
        ax.set_ylim(-gap * 1.15, gap * 1.15)
        ax.set_zlim(-z_half,    z_half)
        ax.tick_params(labelsize=7)

        # Build a compact legend (deduplicate auto-entries, add ion colour proxies)
        from matplotlib.lines import Line2D  # noqa: PLC0415
        handles, labs = ax.get_legend_handles_labels()
        by_label = dict(zip(labs, handles))
        by_label["Ion → readout"]     = Line2D([0], [0], color="limegreen", lw=1.5)
        by_label["Ion → non-readout"] = Line2D([0], [0], color="magenta",   lw=1.5)
        ax.legend(by_label.values(), by_label.keys(),
                  loc="upper right", fontsize=7, framealpha=0.6)

        if self.tracks_canvas._user_dist is not None:
            self.tracks_canvas.ax.dist = self.tracks_canvas._user_dist
        self.tracks_canvas.figure.tight_layout()
        self.tracks_canvas.draw_idle()


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._runner: SimRunner | None = None
        self._last_loaded_config_path: str | None = None

        binary_ok = BINARY.exists()
        suffix = "" if binary_ok else "  ⚠ Binary not found — build first"
        self.setWindowTitle(f"TGC Simulation{suffix}")
        self.resize(1200, 740)

        # ── Toolbar ───────────────────────────────────────────────────────
        tb = QToolBar("Main")
        tb.setMovable(False)
        self.addToolBar(tb)

        self.act_run  = tb.addAction("▶  Run")
        self.act_stop = tb.addAction("■  Stop")
        tb.addSeparator()
        self.act_load = tb.addAction("Load Config")
        self.act_save = tb.addAction("Save Config")

        self.act_run.setEnabled(binary_ok)
        self.act_stop.setEnabled(False)

        self.act_run.triggered.connect(self._on_run)
        self.act_stop.triggered.connect(self._on_stop)
        self.act_load.triggered.connect(self._on_load_config)
        self.act_save.triggered.connect(self._on_save_config)

        # ── Central splitter ─────────────────────────────────────────────
        splitter = QSplitter(Qt.Horizontal)

        self.config_panel  = ConfigPanel()
        self.results_panel = ResultsPanel()

        splitter.addWidget(self.config_panel)
        splitter.addWidget(self.results_panel)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([350, 850])

        self.setCentralWidget(splitter)

        # ── Status bar ────────────────────────────────────────────────────
        if binary_ok:
            self.statusBar().showMessage(f"Ready — binary: {BINARY}")
        else:
            self.statusBar().showMessage(
                f"Binary not found at {BINARY} — run cmake to build"
            )

    # ── Toolbar actions ───────────────────────────────────────────────────

    def _on_run(self):
        if not BINARY.exists():
            QMessageBox.critical(
                self, "Binary not found",
                f"tgc_sim binary not found at:\n{BINARY}\n\n"
                "Build the project first:\n"
                "  cmake -S projects/tgc -B projects/tgc/build ...\n"
                "  cmake --build projects/tgc/build -j4"
            )
            return

        cfg     = self.config_panel.to_config_dict()
        out_str = self.config_panel.out_dir.text().strip() or "results"

        # Resolve relative paths from the tgc project directory
        out_path = Path(out_str)
        if not out_path.is_absolute():
            out_path = (TGC_DIR / out_str).resolve()
        out_path.mkdir(parents=True, exist_ok=True)

        self.results_panel.clear_log()
        self.results_panel.setCurrentIndex(0)   # show Log tab while running
        self.statusBar().showMessage("Running…")

        self._runner = SimRunner(cfg, str(out_path))
        self._runner.log_line.connect(self.results_panel.append_log)
        self._runner.finished.connect(self._on_run_finished)
        self._runner.failed.connect(self._on_run_failed)

        self.act_run.setEnabled(False)
        self.act_stop.setEnabled(True)

        if self._last_loaded_config_path:
            self.results_panel.append_log(
                f"[GUI] Config based on: {self._last_loaded_config_path}"
            )
        else:
            self.results_panel.append_log("[GUI] Config from widget defaults")

        self._runner.start()

    def _on_stop(self):
        if self._runner:
            self._runner.stop()
        self.act_run.setEnabled(True)
        self.act_stop.setEnabled(False)
        self.statusBar().showMessage("Stopped by user")

    def _on_run_finished(self, run_dir: str):
        self.act_run.setEnabled(True)
        self.act_stop.setEnabled(False)
        self.statusBar().showMessage(f"Done — output in {run_dir}")

        self.results_panel.append_log(f"\n[GUI] Simulation complete.  Output: {run_dir}")

        csv_path  = str(Path(run_dir) / "summary.csv")
        root_path = str(Path(run_dir) / "tgc_sim.root")

        self.results_panel.populate_table(csv_path)
        self.results_panel.draw_plots(csv_path)
        self.results_panel.load_waveform_data(root_path)
        self.results_panel.load_track_data(root_path, run_dir)
        self.results_panel.setCurrentIndex(1)   # switch to Summary tab

    def _on_run_failed(self, msg: str):
        self.act_run.setEnabled(True)
        self.act_stop.setEnabled(False)
        self.statusBar().showMessage(f"Failed: {msg}")
        self.results_panel.append_log(f"\n[GUI] ERROR: {msg}")
        QMessageBox.warning(self, "Simulation failed", msg)

    # ── Config load/save ──────────────────────────────────────────────────

    def _on_load_config(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Load config", str(TGC_DIR / "config"),
            "JSON files (*.json);;All files (*)"
        )
        if not path:
            return
        try:
            with open(path) as f:
                d = json.load(f)
            self.config_panel.load_from_dict(d)
            self._last_loaded_config_path = path
            self.statusBar().showMessage(f"Config loaded from {path}")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Load failed", str(exc))

    def _on_save_config(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save config", str(TGC_DIR / "config"),
            "JSON files (*.json);;All files (*)"
        )
        if not path:
            return
        try:
            cfg = self.config_panel.to_config_dict()
            with open(path, "w") as f:
                json.dump(cfg, f, indent=2)
            self.statusBar().showMessage(f"Config saved to {path}")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Save failed", str(exc))

    # ── Window lifecycle ──────────────────────────────────────────────────

    def closeEvent(self, event):
        if self._runner and self._runner.isRunning():
            self._runner.stop()
            self._runner.wait(3000)
        super().closeEvent(event)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("TGC Simulation")
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
