#!/usr/bin/env python3
"""
TGC Simulation GUI
A PyQt5 desktop application for configuring, running, and displaying results
from the tgc_sim Garfield++ binary.

Launch from anywhere:
    python3 projects/tgc/gui/app.py
"""

import json
import math
import os
import subprocess
import sys
import tempfile
from datetime import datetime
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
    g1  = gas.get("gas1", "ar")
    f1  = round(gas.get("gas1_fraction_pct", 70.0))
    g2  = gas.get("gas2", "co2")
    f2  = 100 - f1
    T   = round(gas.get("temperature_K", 293.15))
    P   = round(gas.get("pressure_Torr", 760.0))
    Ee  = round(gas.get("max_electron_energy_eV", 2000.0))
    Ef    = round(gas.get("e_field_max_vcm", 300000.0) / 1000)
    EfMin = round(gas.get("e_field_min_vcm", 100.0))
    n   = gas.get("n_field_points", 20)
    c   = gas.get("n_magboltz_collisions", 10)
    pen = "pen" if gas.get("enable_penning", True) else "nopen"
    return f"{g1}{f1}_{g2}_{f2}_T{T}_P{P}_Ee{Ee}_Ef{EfMin}v-{Ef}k_n{n}_c{c}_{pen}.gas"


def derive_gas_props_filename(gas: dict) -> str:
    """Return the sidecar CSV filename for Magboltz transport properties."""
    return derive_gas_filename(gas).replace(".gas", "_props.csv")


# ---------------------------------------------------------------------------
# Background simulation runner
# ---------------------------------------------------------------------------

class SimRunner(QThread):
    """Runs tgc_sim in a background thread and emits stdout line-by-line."""

    log_line = pyqtSignal(str)   # one stdout line
    finished = pyqtSignal(str)   # emits the run output directory on success
    failed   = pyqtSignal(str)   # emits an error message on failure

    def __init__(self, config_dict: dict, out_dir: str,
                 run_name: str = "", parent=None):
        super().__init__(parent)
        self._config    = config_dict
        self._out_dir   = out_dir
        self._run_name  = run_name          # passed as --run-name to binary
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
        if self._run_name:
            cmd += ["--run-name", self._run_name]

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

        # Locate the sub-directory the binary created.
        # If we passed --run-name we already know the exact path; otherwise
        # glob for any folder matching the naming pattern (supports both the
        # new date-prefixed names and old-style V<V>V__n<N> directories).
        out_path = Path(self._out_dir)
        if self._run_name:
            run_dir = str(out_path / self._run_name)
        else:
            subdirs = sorted(out_path.glob("*__n*"),
                             key=lambda p: p.stat().st_mtime)
            run_dir = str(subdirs[-1]) if subdirs else self._out_dir
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

        self.all_sense_wires = QCheckBox("All wires read out")
        self.all_sense_wires.setChecked(True)
        self.sense_wires = QLineEdit("")
        self.sense_wires.setEnabled(False)
        self.sense_wires.setPlaceholderText("e.g. 4,5  (0-based; 0 = leftmost)")
        self.sense_wires.setToolTip(
            "Comma-separated 0-based indices of the wires summed into the anode "
            "readout. Other wires stay at HV and shape the field but are not read out.")
        self.all_sense_wires.toggled.connect(
            lambda on: self.sense_wires.setEnabled(not on))

        geo_form.addRow("Wire pitch [cm]",    self.wire_pitch)
        geo_form.addRow("Wire diameter [μm]", self.wire_diam)
        geo_form.addRow("Gap [cm]",           self.gap_cm)
        geo_form.addRow("N wires",            self.n_wires)
        geo_form.addRow("Wire voltage [V]",   self.wire_volts)
        geo_form.addRow("Sense wires",        self.all_sense_wires)
        geo_form.addRow("",                   self.sense_wires)
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

        self.delayed_signal_cb = QCheckBox()
        self.delayed_signal_cb.setChecked(True)
        self.delayed_signal_cb.setToolTip(
            "When unchecked, skips the time-varying delayed weighting potential.\n"
            "The static α-corrected weighting potential is still applied.\n"
            "Disabling removes the ~200× per-step overhead and gives conductive-like speed."
        )

        ro_form.addRow("Type",                   self.readout_type)
        ro_form.addRow("Insulator material",     self.insulator_material)
        ro_form.addRow("Thickness [μm]",         self.insulator_thickness)
        ro_form.addRow("Resistivity [kΩ/sq]",    self.surface_resistivity)
        ro_form.addRow("Delayed signal",         self.delayed_signal_cb)
        root_layout.addWidget(ro_box)

        self._update_readout_widgets()
        self.readout_type.currentIndexChanged.connect(self._update_readout_widgets)

        # ── Source ────────────────────────────────────────────────────────
        src_box  = QGroupBox("Source")
        src_form = QFormLayout(src_box)

        self.energy_kev = self._dspin(0.1, 100.0, 0.1, 2, 5.9)

        self.dist_random = QCheckBox("Random (uniform over gap)")
        self.dist_random.setChecked(False)
        self.distances  = QLineEdit("0.2,0.5,0.9,1.2")
        self.distances.setToolTip("Comma-separated source y-distances from wire plane [mm]")
        self.dist_random.toggled.connect(lambda on: self.distances.setEnabled(not on))

        self.x_random = QCheckBox("Random (uniform over wire span)")
        self.x_random.setChecked(True)
        self.x_positions = QLineEdit("0.0")
        self.x_positions.setEnabled(False)
        self.x_positions.setToolTip(
            "Comma-separated fixed x-positions [cm] (e.g. 0.0, 0.09, 0.18)")
        self.x_random.toggled.connect(lambda on: self.x_positions.setEnabled(not on))

        src_form.addRow("Energy [keV]",      self.energy_kev)
        src_form.addRow("Distance",          self.dist_random)
        src_form.addRow("  fixed dist [mm]", self.distances)
        src_form.addRow("X position",        self.x_random)
        src_form.addRow("  fixed x [cm]",    self.x_positions)
        root_layout.addWidget(src_box)

        # ── Gas ───────────────────────────────────────────────────────────
        gas_box  = QGroupBox("Gas")
        gas_form = QFormLayout(gas_box)

        # — Composition rows —
        _GAS_LIST = ["ar", "co2", "cf4", "ch4", "c2h6", "n2", "he", "ne"]
        _ION_LIST  = ["ar", "co2", "cf4", "he", "ne"]  # species with IonMobility files

        self.gas1_combo = QComboBox()
        self.gas1_combo.addItems(_GAS_LIST)
        self.gas1_combo.setEditable(True)
        self.gas1_combo.setCurrentText("ar")
        self.gas1_combo.setToolTip("First gas component (Magboltz species name, lowercase)")

        self.frac1_spin = QDoubleSpinBox()
        self.frac1_spin.setRange(1.0, 99.0)
        self.frac1_spin.setSingleStep(1.0)
        self.frac1_spin.setDecimals(1)
        self.frac1_spin.setValue(70.0)
        self.frac1_spin.setSuffix(" %")

        gas1_row = QWidget()
        gas1_h   = QHBoxLayout(gas1_row)
        gas1_h.setContentsMargins(0, 0, 0, 0)
        gas1_h.addWidget(self.gas1_combo)
        gas1_h.addWidget(self.frac1_spin)
        gas_form.addRow("Gas 1 [%]", gas1_row)

        self.gas2_combo = QComboBox()
        self.gas2_combo.addItems(_GAS_LIST)
        self.gas2_combo.setEditable(True)
        self.gas2_combo.setCurrentText("co2")
        self.gas2_combo.setToolTip("Second gas component (Magboltz species name, lowercase)")

        self.gas2_frac_lbl = QLabel("30.0 %")
        self.gas2_frac_lbl.setToolTip("Fraction of gas 2 = 100% − gas 1 fraction (auto-computed)")

        gas2_row = QWidget()
        gas2_h   = QHBoxLayout(gas2_row)
        gas2_h.setContentsMargins(0, 0, 0, 0)
        gas2_h.addWidget(self.gas2_combo)
        gas2_h.addWidget(self.gas2_frac_lbl)
        gas_form.addRow("Gas 2 [%]", gas2_row)

        self.ion_combo = QComboBox()
        self.ion_combo.addItems(_ION_LIST)
        self.ion_combo.setCurrentText("co2")
        self.ion_combo.setToolTip(
            "Ion species for the mobility table (IonMobility_X+_X.txt).\n"
            "Available: ar, co2, cf4, he, ne.\n"
            "Should match the dominant drifting ion in the mixture."
        )
        gas_form.addRow("Ion species", self.ion_combo)

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
        self.e_field_min = self._dspin(10.0, 100_000.0, 100.0, 0, 100.0)
        self.e_field_min.setToolTip(
            "Minimum E-field in the Magboltz table [V/cm].\n"
            "100 V/cm is suitable for most TGC operating conditions."
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
        gas_form.addRow("E-field min [V/cm]", self.e_field_min)
        ef_row = QWidget()
        ef_h   = QHBoxLayout(ef_row)
        ef_h.setContentsMargins(0, 0, 0, 0)
        ef_h.addWidget(self.e_field_max)
        btn_ef_auto = QPushButton("Auto")
        btn_ef_auto.setFixedWidth(44)
        btn_ef_auto.setToolTip(
            "Set to 2× estimated peak near-wire field (Sauli 1977 formula).\n"
            "E_peak = V / (r × (ln(pitch/(2π r)) + π gap/pitch))"
        )
        btn_ef_auto.clicked.connect(self._on_ef_auto)
        ef_h.addWidget(btn_ef_auto)
        gas_form.addRow("E-field max [V/cm]",  ef_row)
        gas_form.addRow("Gas file (auto)",     self.gas_file_label)
        root_layout.addWidget(gas_box)

        # Update gas file label whenever a gas or geometry parameter changes
        self.gas1_combo.currentTextChanged.connect(self._update_gas2_frac_label)
        self.frac1_spin.valueChanged.connect(self._update_gas2_frac_label)
        self.gas1_combo.currentTextChanged.connect(self._update_gas_file_label)
        self.frac1_spin.valueChanged.connect(self._update_gas_file_label)
        self.gas2_combo.currentTextChanged.connect(self._update_gas_file_label)
        self.ion_combo.currentTextChanged.connect(self._update_gas_file_label)
        self.temperature.valueChanged.connect(self._update_gas_file_label)
        self.pressure.valueChanged.connect(self._update_gas_file_label)
        self.penning.toggled.connect(self._update_gas_file_label)
        self.ncoll.valueChanged.connect(self._update_gas_file_label)
        self.max_electron_energy.valueChanged.connect(self._update_gas_file_label)
        self.n_field_pts.valueChanged.connect(self._update_gas_file_label)
        self.e_field_max.valueChanged.connect(self._update_gas_file_label)
        # Geometry changes affect E_peak shown in the label
        self.wire_pitch.valueChanged.connect(self._update_gas_file_label)
        self.wire_diam.valueChanged.connect(self._update_gas_file_label)
        self.gap_cm.valueChanged.connect(self._update_gas_file_label)
        self.wire_volts.valueChanged.connect(self._update_gas_file_label)

        # ── Simulation ────────────────────────────────────────────────────
        sim_box  = QGroupBox("Simulation")
        sim_form = QFormLayout(sim_box)

        self.n_events    = self._spin(1, 100000, 1000)
        self.max_aval    = self._spin(1000, 10000000, 500000)
        self.time_window = self._dspin(10.0, 100000.0, 10.0, 1, 300.0)
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

        self.run_name = QLineEdit()
        self.run_name.setPlaceholderText("auto  (date + voltage + events)")
        self.run_name.setToolTip(
            "Optional label for the run subfolder.\n"
            "Leave blank: yymmdd_hh-mm__VφV__nη  (auto)\n"
            "Filled:      yymmdd_hh-mm__<your label>")
        out_form.addRow("Run name", self.run_name)

        root_layout.addWidget(out_box)

        root_layout.addStretch()
        self.setWidget(container)

        self._update_gas2_frac_label()
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

    def _update_gas2_frac_label(self):
        """Keep the gas-2 fraction label in sync with gas-1 fraction spinner."""
        self.gas2_frac_lbl.setText(f"{100.0 - self.frac1_spin.value():.1f} %")

    def _update_gas_file_label(self):
        gas = {
            "gas1":                   self.gas1_combo.currentText().strip().lower(),
            "gas1_fraction_pct":      self.frac1_spin.value(),
            "gas2":                   self.gas2_combo.currentText().strip().lower(),
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

        e_peak_kvcm = self._compute_peak_field_kvcm()
        e_max_kvcm  = self.e_field_max.value() / 1000.0
        margin      = e_max_kvcm / e_peak_kvcm if e_peak_kvcm > 0 else float("inf")
        if e_max_kvcm < e_peak_kvcm:
            margin_str = " ⚠ below E_peak!"
        elif margin < 1.5:
            margin_str = f" (margin {margin:.1f}×)"
        else:
            margin_str = ""
        self.gas_file_label.setText(
            f"{name}\n[{status}]  E_peak ≈ {e_peak_kvcm:.0f} kV/cm{margin_str}"
        )

    def _compute_peak_field_kvcm(self) -> float:
        """Sauli 1977 peak near-wire field estimate [kV/cm]."""
        r_cm  = self.wire_diam.value() * 0.5e-4          # radius in cm
        pitch = self.wire_pitch.value()
        gap   = self.gap_cm.value()
        volt  = self.wire_volts.value()
        cap   = math.log(pitch / (2 * math.pi * r_cm)) + math.pi * gap / pitch
        return volt / (r_cm * cap) / 1000.0               # V/cm → kV/cm

    def _on_ef_auto(self):
        """Fill e_field_max with 2× E_peak, rounded up to the nearest 50 kV/cm."""
        e_peak_vcm = self._compute_peak_field_kvcm() * 1000.0
        auto_val   = math.ceil(2.0 * e_peak_vcm / 50_000.0) * 50_000.0
        self.e_field_max.setValue(auto_val)

    def _update_readout_widgets(self):
        resistive = self.readout_type.currentText() == "Resistive"
        self.insulator_material.setEnabled(resistive)
        self.insulator_thickness.setEnabled(resistive)
        self.surface_resistivity.setEnabled(resistive)
        self.delayed_signal_cb.setEnabled(resistive)

    # ── file dialogs ─────────────────────────────────────────────────────

    def _browse_out_dir(self):
        path = QFileDialog.getExistingDirectory(self, "Select output directory")
        if path:
            self.out_dir.setText(path)

    # ── serialisation ─────────────────────────────────────────────────────

    def to_config_dict(self) -> dict:
        """Assemble widget values into a config dict suitable for JSON dump."""
        if self.dist_random.isChecked():
            dists = None  # → JSON null → C++ random per-event
        else:
            raw = self.distances.text().strip()
            try:
                dists = [float(x.strip()) for x in raw.split(",") if x.strip()]
            except ValueError:
                dists = [0.7]

        if self.x_random.isChecked():
            x_positions = None
        else:
            raw_x = self.x_positions.text().strip()
            try:
                x_positions = [float(v.strip()) for v in raw_x.split(",") if v.strip()]
            except ValueError:
                x_positions = [0.0]

        if self.all_sense_wires.isChecked():
            sense_wires = None
        else:
            nmax = self.n_wires.value()
            try:
                sense_wires = sorted({
                    int(v.strip())
                    for v in self.sense_wires.text().split(",")
                    if v.strip() and 0 <= int(v.strip()) < nmax
                })
            except ValueError:
                sense_wires = None
            if not sense_wires:
                sense_wires = None

        ro_type = self.readout_type.currentText().lower()
        ins_mat = self.insulator_material.currentText().lower()

        return {
            "geometry": {
                "wire_pitch_cm":    self.wire_pitch.value(),
                "wire_diameter_um": self.wire_diam.value(),
                "gap_cm":           self.gap_cm.value(),
                "n_wires":          self.n_wires.value(),
                "wire_voltage_V":   self.wire_volts.value(),
                "sense_wires":      sense_wires,
            },
            "readout": {
                "type":                       ro_type,
                "insulator_material":         ins_mat,
                "insulator_thickness_um":     self.insulator_thickness.value(),
                "surface_resistivity_ohm_sq": self.surface_resistivity.value() * 1000.0,
                "enable_delayed_signal":      self.delayed_signal_cb.isChecked(),
            },
            "source": {
                "energy_keV":          self.energy_kev.value(),
                "source_distances_mm": dists,
                "x_positions_cm":      x_positions,
            },
            "gas": {
                "gas1":                   self.gas1_combo.currentText().strip().lower(),
                "gas1_fraction_pct":      self.frac1_spin.value(),
                "gas2":                   self.gas2_combo.currentText().strip().lower(),
                "ion_species":            self.ion_combo.currentText().strip().lower(),
                "temperature_K":          self.temperature.value(),
                "pressure_Torr":          self.pressure.value(),
                "enable_penning":         self.penning.isChecked(),
                "n_magboltz_collisions":  self.ncoll.value(),
                "w_value_eV":             self.w_value.value(),
                "max_electron_energy_eV": self.max_electron_energy.value(),
                "n_field_points":         self.n_field_pts.value(),
                "e_field_min_vcm":        self.e_field_min.value(),
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
        sw = g.get("sense_wires", None)
        if sw:
            self.all_sense_wires.setChecked(False)
            self.sense_wires.setText(",".join(str(int(i)) for i in sw))
        else:
            self.all_sense_wires.setChecked(True)
            self.sense_wires.clear()

        ro = d.get("readout", {})
        ro_type = ro.get("type", "conductive")
        self.readout_type.setCurrentIndex(0 if ro_type == "conductive" else 1)
        ins_mat = ro.get("insulator_material", "kapton")
        self.insulator_material.setCurrentIndex(0 if ins_mat == "kapton" else 1)
        self.insulator_thickness.setValue(ro.get("insulator_thickness_um", 100.0))
        self.surface_resistivity.setValue(ro.get("surface_resistivity_ohm_sq", 500000.0) / 1000.0)
        self.delayed_signal_cb.setChecked(ro.get("enable_delayed_signal", True))

        s = d.get("source", {})
        self.energy_kev.setValue(s.get("energy_keV", 5.9))
        dists = s.get("source_distances_mm", [0.2, 0.5, 0.9, 1.2])
        if dists is None:
            self.dist_random.setChecked(True)
        else:
            self.dist_random.setChecked(False)
            self.distances.setText(",".join(str(v) for v in dists))
        x_positions = s.get("x_positions_cm", None)
        if x_positions is None:
            # backward compat: old scalar key
            scalar = s.get("x_position_cm", None)
            x_positions = [float(scalar)] if scalar is not None else None
        if x_positions is None:
            self.x_random.setChecked(True)
        else:
            self.x_random.setChecked(False)
            self.x_positions.setText(",".join(str(v) for v in x_positions))

        gas = d.get("gas", {})
        self.gas1_combo.setCurrentText(gas.get("gas1", "ar"))
        self.frac1_spin.setValue(       gas.get("gas1_fraction_pct", 70.0))
        self.gas2_combo.setCurrentText(gas.get("gas2", "co2"))
        self.ion_combo.setCurrentText( gas.get("ion_species", "co2"))
        self._update_gas2_frac_label()
        self.temperature.setValue(gas.get("temperature_K", 293.15))
        self.pressure.setValue(   gas.get("pressure_Torr", 760.0))
        self.penning.setChecked(  gas.get("enable_penning", True))
        self.ncoll.setValue(      gas.get("n_magboltz_collisions", 10))
        self.w_value.setValue(    gas.get("w_value_eV", 26.0))
        self.max_electron_energy.setValue(gas.get("max_electron_energy_eV", 2000.0))
        self.n_field_pts.setValue(        gas.get("n_field_points", 20))
        self.e_field_min.setValue(        gas.get("e_field_min_vcm", 100.0))
        self.e_field_max.setValue(        gas.get("e_field_max_vcm", 300_000.0))

        sim = d.get("simulation", {})
        self.n_events.setValue(        sim.get("n_events", 1000))
        self.max_aval.setValue(        sim.get("max_avalanche_size", 500000))
        self.time_window.setValue(     sim.get("time_window_ns", 300.0))
        self.time_step.setValue(       sim.get("time_step_ns", 0.5))
        self.enable_ion_drift.setChecked(sim.get("enable_ion_drift", True))
        self.store_drift_lines.setChecked(sim.get("store_drift_lines", True))


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
        _sum_widget = QWidget()
        _sum_vbox   = QVBoxLayout(_sum_widget)
        _sum_vbox.setContentsMargins(4, 4, 4, 4)

        _sum_btn_row = QWidget()
        _sum_btn_h   = QHBoxLayout(_sum_btn_row)
        _sum_btn_h.setContentsMargins(0, 0, 0, 0)
        self.summary_export_btn = QPushButton("Export CSV …")
        self.summary_export_btn.setEnabled(False)
        _sum_btn_h.addWidget(self.summary_export_btn)
        _sum_btn_h.addStretch()
        _sum_vbox.addWidget(_sum_btn_row)

        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        _sum_vbox.addWidget(self.table)

        self.addTab(_sum_widget, "Summary")
        self.summary_export_btn.clicked.connect(self._on_summary_export_csv)

        # ── Plots tab: 2×3 matplotlib figure ─────────────────────────────
        self.plots_canvas = MplCanvas(nrows=2, ncols=3, figsize=(11, 5))
        self.addTab(self.plots_canvas, "Plots")

        # ── Waveforms tab: ROOT TCanvas browser ──────────────────────────
        self._waveform_data: dict = {}
        self._root_canvas  = None    # ROOT TCanvas (kept alive between events)
        self._root_objects: list = []  # TGraph/TLegend objects (prevent Python GC)
        self._charge_canvas  = None   # ROOT TCanvas for charge integrals
        self._charge_objects: list = []
        self._track_data:   dict = {}  # {dist_label: {xpos_label: data_dict}}
        self._track_geom:   dict | None = None
        self._tracks_canvas  = None   # ROOT TCanvas for 3D tracks
        self._tracks_objects: list = []
        self._trk_legend_objects: list = []   # TLegend + proxy objects (prevent Python GC)
        self._trk_zoom_scale: float = 1.0   # <1 zoomed in, >1 zoomed out
        self._trk_view_phi:   float = 32.0  # TPad azimuthal angle (32° avoids Y-label inside box)
        self._trk_view_theta: float = 30.0  # TPad elevation angle (ROOT default)
        self._trk_pan_x:      float = 0.0   # cm offset of X visible centre
        self._trk_pan_y:      float = 0.0   # cm offset of Y visible centre
        self._trk_pan_z:      float = 0.0   # cm offset of Z visible centre
        self._efield_cache: dict | None = None   # {x, y, Ex, Ey} computed arrays
        self._efield_root_canvas  = None   # ROOT TCanvas for E-field maps
        self._efield_objects: list = []

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
        sel_h.addSpacing(8)
        sel_h.addWidget(QLabel("X pos:"))
        self.wave_xpos_combo = QComboBox()
        sel_h.addWidget(self.wave_xpos_combo)
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
        self.wave_xpos_combo.currentIndexChanged.connect(self._on_wave_xpos_changed)
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
        ch_sel_h.addSpacing(8)
        ch_sel_h.addWidget(QLabel("X pos:"))
        self.charge_xpos_combo = QComboBox()
        ch_sel_h.addWidget(self.charge_xpos_combo)
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
        self.charge_xpos_combo.currentIndexChanged.connect(self._on_charge_xpos_changed)
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
        trk_sel_h.addSpacing(12)
        trk_sel_h.addWidget(QLabel("X pos:"))
        self.trk_xpos_combo = QComboBox()
        trk_sel_h.addWidget(self.trk_xpos_combo)
        trk_sel_h.addSpacing(12)
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

        # — view controls (preset orientations + zoom) —
        trk_ctrl_row = QWidget()
        trk_ctrl_h   = QHBoxLayout(trk_ctrl_row)
        trk_ctrl_h.setContentsMargins(0, 0, 0, 0)

        for _label, _phi, _theta, _rz in [
            ("Gap XY",   0,  90, False),  # theta=90 → camera along Z → sees X-Y
            ("Top XZ",   0,   0, False),  # phi=0   → top down (along Y) → sees X-Z
            ("Side YZ",  90,  0, False),  # phi=90  → camera along X → sees Y-Z
            ("3D",      32,  30, True),   # default perspective + reset zoom
        ]:
            _btn = QPushButton(_label)
            _btn.setMaximumWidth(72)
            _btn.clicked.connect(
                (lambda p, t, rz: lambda: self._trk_preset_view(p, t, rz))
                (_phi, _theta, _rz))
            trk_ctrl_h.addWidget(_btn)

        trk_ctrl_h.addSpacing(16)
        trk_zoom_in_btn  = QPushButton("Zoom +")
        trk_zoom_out_btn = QPushButton("Zoom −")
        trk_zoom_in_btn.setMaximumWidth(65)
        trk_zoom_out_btn.setMaximumWidth(65)
        trk_zoom_in_btn.clicked.connect(lambda: self._trk_adjust_zoom(0.7))
        trk_zoom_out_btn.clicked.connect(lambda: self._trk_adjust_zoom(1.0 / 0.7))
        trk_ctrl_h.addWidget(trk_zoom_in_btn)
        trk_ctrl_h.addWidget(trk_zoom_out_btn)
        trk_ctrl_h.addStretch()
        tracks_layout.addWidget(trk_ctrl_row)

        # — pan controls —
        trk_pan_row = QWidget()
        trk_pan_h   = QHBoxLayout(trk_pan_row)
        trk_pan_h.setContentsMargins(0, 0, 0, 0)
        trk_pan_h.addWidget(QLabel("Pan:"))
        for _ax, _d, _lbl in [
            ("x", -1, "X-"), ("x", +1, "X+"),
            ("y", -1, "Y-"), ("y", +1, "Y+"),
            ("z", -1, "Z-"), ("z", +1, "Z+"),
        ]:
            _btn = QPushButton(_lbl)
            _btn.setMaximumWidth(42)
            _btn.clicked.connect(
                (lambda a, d: lambda: self._trk_pan(a, d))(_ax, _d))
            trk_pan_h.addWidget(_btn)
        trk_pan_h.addStretch()
        tracks_layout.addWidget(trk_pan_row)

        trk_hint = QLabel(
            "ROOT canvas opens automatically when results are loaded.\n"
            "Left-click drag: rotate  ·  Zoom / Pan / View buttons above  ·  Right-click: save."
        )
        trk_hint.setWordWrap(True)
        trk_hint.setStyleSheet("color: grey; font-size: 11px;")
        tracks_layout.addWidget(trk_hint)
        tracks_layout.addStretch()

        self.trk_dist_combo.currentIndexChanged.connect(self._on_trk_dist_changed)
        self.trk_xpos_combo.currentIndexChanged.connect(self._on_trk_xpos_changed)
        self.trk_event_slider.sliderReleased.connect(self._update_track_plot)

        self.addTab(tracks_widget, "3D Tracks")

        # ── E-Field tab ───────────────────────────────────────────────────────
        efield_widget = QWidget()
        efield_layout = QVBoxLayout(efield_widget)
        efield_layout.setContentsMargins(8, 6, 8, 6)
        efield_layout.setSpacing(6)

        # — row 1: geometry inputs —
        geom_row = QWidget()
        geom_h   = QHBoxLayout(geom_row)
        geom_h.setContentsMargins(0, 0, 0, 0)

        def _add_geom_spin(label, widget):
            geom_h.addWidget(QLabel(label))
            geom_h.addWidget(widget)
            geom_h.addSpacing(8)

        self.ef_gap       = QDoubleSpinBox()
        self.ef_gap.setRange(0.01, 5.0);  self.ef_gap.setSingleStep(0.01)
        self.ef_gap.setDecimals(3);       self.ef_gap.setValue(0.14)
        self.ef_pitch     = QDoubleSpinBox()
        self.ef_pitch.setRange(0.01, 5.0); self.ef_pitch.setSingleStep(0.01)
        self.ef_pitch.setDecimals(3);      self.ef_pitch.setValue(0.18)
        self.ef_n_wires   = QSpinBox()
        self.ef_n_wires.setRange(2, 200);  self.ef_n_wires.setValue(10)
        self.ef_wire_diam = QDoubleSpinBox()
        self.ef_wire_diam.setRange(1.0, 1000.0); self.ef_wire_diam.setSingleStep(1.0)
        self.ef_wire_diam.setDecimals(1);         self.ef_wire_diam.setValue(50.0)
        self.ef_wire_volt = QDoubleSpinBox()
        self.ef_wire_volt.setRange(100.0, 10000.0); self.ef_wire_volt.setSingleStep(100.0)
        self.ef_wire_volt.setDecimals(0);           self.ef_wire_volt.setValue(1900.0)

        _add_geom_spin("Gap [cm]:",         self.ef_gap)
        _add_geom_spin("Pitch [cm]:",       self.ef_pitch)
        _add_geom_spin("N wires:",          self.ef_n_wires)
        _add_geom_spin("Wire diam [µm]:",   self.ef_wire_diam)
        _add_geom_spin("Wire voltage [V]:", self.ef_wire_volt)
        geom_h.addStretch()
        efield_layout.addWidget(geom_row)

        # — row 2: slice depths + component + colormap + compute button —
        ctrl_row = QWidget()
        ctrl_h   = QHBoxLayout(ctrl_row)
        ctrl_h.setContentsMargins(0, 0, 0, 0)

        ctrl_h.addWidget(QLabel("ZX plane at y [cm]:"))
        self.ef_y_depth = QDoubleSpinBox()
        self.ef_y_depth.setRange(-0.14, 0.14); self.ef_y_depth.setSingleStep(0.01)
        self.ef_y_depth.setDecimals(3);        self.ef_y_depth.setValue(0.0)
        ctrl_h.addWidget(self.ef_y_depth)
        ctrl_h.addSpacing(12)

        ctrl_h.addWidget(QLabel("ZY plane at x [cm]:"))
        self.ef_x_depth = QDoubleSpinBox()
        self.ef_x_depth.setRange(-1.0, 1.0);  self.ef_x_depth.setSingleStep(0.01)
        self.ef_x_depth.setDecimals(3);        self.ef_x_depth.setValue(0.0)
        ctrl_h.addWidget(self.ef_x_depth)
        ctrl_h.addSpacing(12)

        ctrl_h.addWidget(QLabel("Component:"))
        self.ef_component = QComboBox()
        self.ef_component.addItems(["|E|", "Ex", "Ey"])
        ctrl_h.addWidget(self.ef_component)
        ctrl_h.addSpacing(8)

        ctrl_h.addWidget(QLabel("Colormap:"))
        self.ef_cmap = QComboBox()
        self.ef_cmap.addItems(["viridis", "plasma", "hot_r", "RdBu_r"])
        ctrl_h.addWidget(self.ef_cmap)
        ctrl_h.addSpacing(12)

        ctrl_h.addWidget(QLabel("Grid nx:"))
        self.ef_nx = QSpinBox()
        self.ef_nx.setRange(50, 10000)
        self.ef_nx.setSingleStep(100)
        self.ef_nx.setValue(250)
        ctrl_h.addWidget(self.ef_nx)
        ctrl_h.addSpacing(6)
        ctrl_h.addWidget(QLabel("ny:"))
        self.ef_ny = QSpinBox()
        self.ef_ny.setRange(50, 10000)
        self.ef_ny.setSingleStep(100)
        self.ef_ny.setValue(150)
        ctrl_h.addWidget(self.ef_ny)
        ctrl_h.addSpacing(12)

        ef_compute_btn = QPushButton("Compute")
        ef_compute_btn.clicked.connect(lambda: self._update_efield_plots(recompute=True))
        ctrl_h.addWidget(ef_compute_btn)
        ctrl_h.addStretch()
        efield_layout.addWidget(ctrl_row)

        # — hint —
        ef_hint = QLabel(
            "ROOT canvas opens automatically when Compute is clicked.\n"
            "Right-click inside the ROOT window to zoom, change axes, or save."
        )
        ef_hint.setWordWrap(True)
        ef_hint.setStyleSheet("color: grey; font-size: 11px;")
        efield_layout.addWidget(ef_hint)
        efield_layout.addStretch()

        # Depth spinboxes re-slice without recomputing the field
        self.ef_y_depth.valueChanged.connect(
            lambda: self._update_efield_plots(recompute=False))
        self.ef_x_depth.valueChanged.connect(
            lambda: self._update_efield_plots(recompute=False))
        self.ef_component.currentIndexChanged.connect(
            lambda: self._update_efield_plots(recompute=False))
        self.ef_cmap.currentIndexChanged.connect(
            lambda: self._update_efield_plots(recompute=False))

        self.addTab(efield_widget, "E-Field")

        # ── Magboltz tab ──────────────────────────────────────────────────
        self._gas_canvas = None
        self._gas_objects: list = []
        self._gas_props_csv: str | None = None

        gas_widget = QWidget()
        gas_layout = QVBoxLayout(gas_widget)
        gas_layout.setContentsMargins(8, 6, 8, 6)
        gas_layout.setSpacing(6)

        exp_row = QWidget()
        exp_h = QHBoxLayout(exp_row)
        exp_h.setContentsMargins(0, 0, 0, 0)
        self.gas_export_root_btn = QPushButton("Export as ROOT …")
        self.gas_export_csv_btn  = QPushButton("Export as CSV …")
        self.gas_export_root_btn.setEnabled(False)
        self.gas_export_csv_btn.setEnabled(False)
        exp_h.addWidget(self.gas_export_root_btn)
        exp_h.addWidget(self.gas_export_csv_btn)
        exp_h.addStretch()
        gas_layout.addWidget(exp_row)

        gas_hint = QLabel(
            "ROOT canvas opens automatically when gas properties are available.\n"
            "Right-click inside the ROOT window to zoom, change axes, or save."
        )
        gas_hint.setWordWrap(True)
        gas_hint.setStyleSheet("color: grey; font-size: 11px;")
        gas_layout.addWidget(gas_hint)
        gas_layout.addStretch()

        self.gas_export_root_btn.clicked.connect(self._on_gas_export_root)
        self.gas_export_csv_btn.clicked.connect(self._on_gas_export_csv)

        self.addTab(gas_widget, "Magboltz")

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
        self._summary_csv_path = csv_path
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
                if isinstance(val, float) and not pd.isna(val):
                    text = f"{val:.4g}"
                elif isinstance(val, float):  # NaN = random x-position
                    text = "—"
                else:
                    text = str(val)
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(r_idx, c_idx, item)

        self.table.resizeColumnsToContents()
        self.summary_export_btn.setEnabled(True)

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

    def draw_gas_props(self, csv_path: str):
        try:
            df = pd.read_csv(csv_path, comment='#')
        except Exception:  # noqa: BLE001
            return
        df.columns = [c.strip() for c in df.columns]
        self._gas_props_csv = csv_path

        try:
            import ROOT  # noqa: PLC0415
            ROOT.gROOT.SetBatch(False)

            if self._root_canvas_alive(self._gas_canvas, "tgc_magboltz"):
                self._gas_canvas.Close()
            self._gas_canvas = ROOT.TCanvas(
                "tgc_magboltz", "Magboltz Gas Properties", 1200, 900)
            self._gas_canvas.Divide(3, 3)
            self._gas_objects.clear()
            ROOT.TGaxis.SetMaxDigits(0)   # force scientific notation for all values < 1

            e = (df["e_field_Vcm"] / 1000.0).to_numpy().astype("f8")   # kV/cm
            n = len(e)

            # Derive ion species and carrier gas from the "# ion_mobility: <file>"
            # comment written by ExportGasProps at the top of the props CSV.
            # File naming convention: IonMobility_<ION>_<GAS>.txt
            _ion_label, _gas_suffix = "ion", ""
            try:
                with open(csv_path) as _f:
                    _first = _f.readline().strip()
                    if _first.startswith("# ion_mobility:"):
                        _stem = Path(_first.split(":", 1)[1].strip()).stem  # "IonMobility_CO2+_CO2"
                        _parts = _stem[len("IonMobility_"):].split("_", 1)  # ["CO2+", "CO2"]
                        if len(_parts) == 2:
                            _ion_label  = _parts[0]            # e.g. "CO2+"
                            _gas_suffix = " in " + _parts[1]   # e.g. " in CO2"
                    else:
                        # Old CSV (no comment) — binary always uses IonMobility_CO2+_CO2.txt
                        _default = (GARFIELD_INSTALL / "share" / "Garfield"
                                    / "Data" / "IonMobility_CO2+_CO2.txt")
                        if _default.exists():
                            _stem = _default.stem  # "IonMobility_CO2+_CO2"
                            _parts = _stem[len("IonMobility_"):].split("_", 1)
                            if len(_parts) == 2:
                                _ion_label  = _parts[0]
                                _gas_suffix = " in " + _parts[1]
            except Exception:  # noqa: BLE001
                pass
            _ion_title_v  = f"{_ion_label} drift velocity{_gas_suffix}"
            _ion_title_mu = f"{_ion_label} mobility{_gas_suffix}"

            panels = [
                # (pad, col,               take_abs, logy, title, xlabel, ylabel)
                (1, "vd_cm_per_us",     True,  False,
                 "Electron drift velocity",
                 "E [kV/cm]", "|v_{d}| [cm/#mus]"),
                (2, "alpha_per_cm",     False, True,
                 "Townsend #alpha",
                 "E [kV/cm]", "#alpha [cm^{-1}]"),
                (3, "eta_per_cm",       False, True,
                 "Attachment #eta",
                 "E [kV/cm]", "#eta [cm^{-1}]"),
                (4, "dl_sqrtcm",        False, False,
                 "Long. diffusion D_{L}",
                 "E [kV/cm]", "D_{L} [cm^{0.5}]"),
                (5, "dt_sqrtcm",        False, False,
                 "Trans. diffusion D_{T}",
                 "E [kV/cm]", "D_{T} [cm^{0.5}]"),
                (6, None,               False, True,
                 "Effective gain (#alpha-#eta)",
                 "E [kV/cm]", "(#alpha-#eta) [cm^{-1}]"),
                (7, "v_ion_cm_per_us",  True,  False,
                 _ion_title_v,
                 "E [kV/cm]", "|v_{ion}| [cm/#mus]"),
                (8, "mu_ion_cm2_per_Vus", False, False,
                 _ion_title_mu,
                 "E [kV/cm]", "#mu [cm^{2}/(V#upoint#mus)]"),
            ]
            for pad_num, col, take_abs, logy, title, xlabel, ylabel in panels:
                self._gas_canvas.cd(pad_num)
                ROOT.gPad.SetLeftMargin(0.18)
                ROOT.gPad.SetBottomMargin(0.16)
                ROOT.gPad.SetRightMargin(0.04)
                ROOT.gPad.SetTopMargin(0.08)
                ROOT.gPad.SetLogx()
                ROOT.gPad.SetGrid()
                if logy:
                    ROOT.gPad.SetLogy()

                if col is None:
                    # effective gain = alpha - eta
                    if "alpha_per_cm" not in df.columns or "eta_per_cm" not in df.columns:
                        continue
                    yraw = (df["alpha_per_cm"] - df["eta_per_cm"]).to_numpy()
                else:
                    if col not in df.columns:
                        continue
                    yraw = df[col].to_numpy()
                    if take_abs:
                        yraw = np.abs(yraw)

                mask = (yraw > 0) if logy else np.ones(n, dtype=bool)
                xe = e[mask].astype("f8")
                ye = yraw[mask].astype("f8")
                if len(xe) == 0:
                    continue

                g = ROOT.TGraph(len(xe), xe, ye)
                g.SetTitle(f"{title};{xlabel};{ylabel}")
                g.SetLineColor(ROOT.kBlue + 1)
                g.SetMarkerColor(ROOT.kBlue + 1)
                g.SetMarkerStyle(20)
                g.SetMarkerSize(0.5)
                g.SetLineWidth(2)
                g.Draw("ALP")
                for _axis in (g.GetXaxis(), g.GetYaxis()):
                    _axis.SetTitleSize(_axis.GetTitleSize() + 0.03)
                    _axis.SetLabelSize(_axis.GetLabelSize() + 0.03)
                self._gas_objects.append(g)

            self._gas_canvas.cd(9)   # leave pad 9 empty

            self._gas_canvas.Update()
            ROOT.TGaxis.SetMaxDigits(5)   # restore ROOT default so other tabs are unaffected
            self._root_timer.start()

            self.gas_export_root_btn.setEnabled(True)
            self.gas_export_csv_btn.setEnabled(True)

        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Magboltz ROOT canvas error: {exc}")

    def _on_gas_export_root(self):
        if not self._gas_props_csv:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Magboltz as ROOT", str(TGC_DIR), "ROOT files (*.root)")
        if not path:
            return
        try:
            import ROOT  # noqa: PLC0415
            df = pd.read_csv(self._gas_props_csv, comment='#')
            df.columns = [c.strip() for c in df.columns]
            e = (df["e_field_Vcm"] / 1000.0).to_numpy().astype("f8")
            f = ROOT.TFile(path, "RECREATE")
            spec = [
                ("vd",    "vd_cm_per_us",       True,
                 "Electron drift velocity;E [kV/cm];|v_{d}| [cm/#mus]"),
                ("alpha", "alpha_per_cm",        False,
                 "Townsend #alpha;E [kV/cm];#alpha [cm^{-1}]"),
                ("eta",   "eta_per_cm",          False,
                 "Attachment #eta;E [kV/cm];#eta [cm^{-1}]"),
                ("dl",    "dl_sqrtcm",           False,
                 "Long. diffusion;E [kV/cm];D_{L} [cm^{0.5}]"),
                ("dt",    "dt_sqrtcm",           False,
                 "Trans. diffusion;E [kV/cm];D_{T} [cm^{0.5}]"),
                ("v_ion", "v_ion_cm_per_us",     True,
                 "Ion drift velocity;E [kV/cm];|v_{ion}| [cm/#mus]"),
                ("mu",    "mu_ion_cm2_per_Vus",  False,
                 "Ion mobility;E [kV/cm];#mu [cm^{2}/(V#upoint#mus)]"),
            ]
            for gname, col, take_abs, title in spec:
                if col not in df.columns:
                    continue
                y = df[col].to_numpy().astype("f8")
                if take_abs:
                    y = np.abs(y)
                g = ROOT.TGraph(len(e), e, y)
                g.SetTitle(title)
                g.SetName(gname)
                g.Write()
            if "alpha_per_cm" in df.columns and "eta_per_cm" in df.columns:
                eff = (df["alpha_per_cm"] - df["eta_per_cm"]).to_numpy().astype("f8")
                g = ROOT.TGraph(len(e), e, eff)
                g.SetTitle("Effective gain;E [kV/cm];(#alpha-#eta) [cm^{-1}]")
                g.SetName("eff_gain")
                g.Write()
            f.Close()
            self.append_log(f"[GUI] Magboltz data exported to {path}")
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] ROOT export failed: {exc}")
            QMessageBox.warning(self.parent(), "Export failed", str(exc))

    def _on_gas_export_csv(self):
        if not self._gas_props_csv or not Path(self._gas_props_csv).exists():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Magboltz as CSV", str(TGC_DIR), "CSV files (*.csv)")
        if not path:
            return
        import shutil
        shutil.copy2(self._gas_props_csv, path)
        self.append_log(f"[GUI] Magboltz CSV exported to {path}")

    def _on_summary_export_csv(self) -> None:
        csv_path = getattr(self, "_summary_csv_path", None)
        if not csv_path or not Path(csv_path).exists():
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export summary as CSV", "summary.csv", "CSV files (*.csv)")
        if not path:
            return
        import shutil
        shutil.copy2(csv_path, path)
        self.append_log(f"[GUI] Summary CSV exported to {path}")

    def _save_plots_root(self, run_dir) -> None:
        """Write all currently-rendered ROOT canvases to tgc_plots.root in run_dir."""
        try:
            import ROOT  # noqa: PLC0415
            out = Path(run_dir) / "tgc_plots.root"
            tf = ROOT.TFile.Open(str(out), "RECREATE")
            if not tf or tf.IsZombie():
                self.append_log("[GUI] Warning: could not create tgc_plots.root")
                return
            saved = []
            for canvas, key in [
                (self._root_canvas,        "waveforms"),
                (self._charge_canvas,      "charge"),
                (self._tracks_canvas,      "tracks_3d"),
                (self._gas_canvas,         "magboltz"),
                (self._efield_root_canvas, "efield"),
            ]:
                if canvas is None:
                    continue
                if not self._root_canvas_alive(canvas, canvas.GetName()):
                    continue
                tf.cd()
                canvas.Write(key)
                saved.append(key)
            tf.Close()
            if saved:
                self.append_log(
                    f"[GUI] Plots saved → {out.name}  ({', '.join(saved)})")
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Could not save plots ROOT file: {exc}")

    # ── Waveforms ─────────────────────────────────────────────────────────

    @staticmethod
    def _root_canvas_alive(canvas, name: str) -> bool:
        """Return True if the ROOT TCanvas still exists (not closed by the user)."""
        try:
            import ROOT  # noqa: PLC0415
            return (canvas is not None and
                    ROOT.gROOT.GetListOfCanvases().FindObject(name) is not None)
        except Exception:  # noqa: BLE001
            return False

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

        # Force canvas recreation for the new run
        self._root_canvas   = None
        self._charge_canvas = None

        self._waveform_data.clear()
        self.wave_dist_combo.blockSignals(True)
        self.wave_dist_combo.clear()
        self.wave_xpos_combo.blockSignals(True)
        self.wave_xpos_combo.clear()
        self.charge_dist_combo.blockSignals(True)
        self.charge_dist_combo.clear()
        self.charge_xpos_combo.blockSignals(True)
        self.charge_xpos_combo.clear()

        try:
            with uproot.open(root_path) as f:
                dist_keys = sorted({
                    k.split("/")[0] for k in f.keys(cycle=False)
                    if k.split("/")[0].startswith("dist_")
                })
                for key in dist_keys:
                    rest = key.removeprefix("dist_")
                    dist_raw, sep, xpos_raw = rest.partition("_x")
                    dist_label = "—" if dist_raw == "rnd" else dist_raw.replace("p", ".").replace("mm", " mm")
                    xpos_label = xpos_raw.replace("p", ".").replace("mm", " mm") if sep else "—"
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

                        self._waveform_data.setdefault(dist_label, {})[xpos_label] = {
                            "times":   times,
                            "anode":   anode,
                            "cathode": cathode,
                            "mean_a":  mean_a,
                            "mean_c":  mean_c,
                        }
                    except Exception as exc:  # noqa: BLE001
                        self.append_log(f"[GUI] Waveforms: could not read {key}: {exc}")
                # Populate dist combos in sorted order (random "—" first)
                for dl in self._sorted_dists(self._waveform_data.keys()):
                    self.wave_dist_combo.addItem(dl)
                    self.charge_dist_combo.addItem(dl)
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Could not open ROOT file: {exc}")

        self.wave_dist_combo.blockSignals(False)
        self.wave_xpos_combo.blockSignals(False)
        self.charge_dist_combo.blockSignals(False)
        self.charge_xpos_combo.blockSignals(False)
        # Triggering wave syncs charge combos + slider; then explicitly draw charge canvas
        if self.wave_dist_combo.count():
            self._on_wave_dist_changed(0)   # syncs charge combos + creates wave canvas
            self._update_charge_plot()      # create charge canvas on initial load
            self._root_timer.start()

    # ── Waveform/Charge sync helpers ──────────────────────────────────────────

    @staticmethod
    def _sorted_xpos(xpos_dict: dict) -> list:
        """Return x-position labels sorted numerically; '—' (random) sorts first."""
        def _key(s):
            try:
                return float(s.replace(" mm", ""))
            except ValueError:
                return -1e9
        return sorted(xpos_dict.keys(), key=_key)

    @staticmethod
    def _sorted_dists(dists_iterable) -> list:
        """Return distance labels sorted numerically; '—' (random) sorts first."""
        def _key(s):
            try:
                return float(s.replace(" mm", ""))
            except ValueError:
                return -1e9
        return sorted(dists_iterable, key=_key)

    def _rebuild_charge_xpos(self):
        """Repopulate charge_xpos_combo to match the current charge_dist selection."""
        dist_label = self.charge_dist_combo.currentText()
        xpos_dict  = self._waveform_data.get(dist_label, {})
        self.charge_xpos_combo.blockSignals(True)
        self.charge_xpos_combo.clear()
        for xp in self._sorted_xpos(xpos_dict):
            self.charge_xpos_combo.addItem(xp)
        self.charge_xpos_combo.blockSignals(False)

    def _rebuild_wave_xpos(self):
        """Repopulate wave_xpos_combo to match the current wave_dist selection."""
        dist_label = self.wave_dist_combo.currentText()
        xpos_dict  = self._waveform_data.get(dist_label, {})
        self.wave_xpos_combo.blockSignals(True)
        self.wave_xpos_combo.clear()
        for xp in self._sorted_xpos(xpos_dict):
            self.wave_xpos_combo.addItem(xp)
        self.wave_xpos_combo.blockSignals(False)

    def _sync_charge_slider(self):
        """Update charge_event_slider range for the current (charge_dist, charge_xpos)."""
        dist_label = self.charge_dist_combo.currentText()
        xpos_label = self.charge_xpos_combo.currentText()
        data = self._waveform_data.get(dist_label, {}).get(xpos_label)
        if data is None:
            return
        n = len(data["anode"])
        self.charge_event_slider.blockSignals(True)
        self.charge_event_slider.setMaximum(max(0, n - 1))
        self.charge_event_slider.setValue(0)
        self.charge_event_slider.blockSignals(False)
        self.charge_event_label.setText(f"1 / {n}")

    def _sync_wave_slider(self):
        """Update wave_event_slider range for the current (wave_dist, wave_xpos)."""
        dist_label = self.wave_dist_combo.currentText()
        xpos_label = self.wave_xpos_combo.currentText()
        data = self._waveform_data.get(dist_label, {}).get(xpos_label)
        if data is None:
            return
        n = len(data["anode"])
        self.wave_event_slider.blockSignals(True)
        self.wave_event_slider.setMaximum(max(0, n - 1))
        self.wave_event_slider.setValue(0)
        self.wave_event_slider.blockSignals(False)
        self.wave_event_label.setText(f"1 / {n}")

    def _on_wave_dist_changed(self, index: int):
        dist_label = self.wave_dist_combo.currentText()
        xpos_dict  = self._waveform_data.get(dist_label, {})

        # Rebuild xpos combo for this distance
        self.wave_xpos_combo.blockSignals(True)
        self.wave_xpos_combo.clear()
        for xp in self._sorted_xpos(xpos_dict):
            self.wave_xpos_combo.addItem(xp)
        self.wave_xpos_combo.blockSignals(False)

        # Sync charge tab (blocked)
        self.charge_dist_combo.blockSignals(True)
        self.charge_dist_combo.setCurrentIndex(index)
        self.charge_dist_combo.blockSignals(False)
        self._rebuild_charge_xpos()

        self._on_wave_xpos_changed(0)

    def _on_wave_xpos_changed(self, index: int):
        dist_label = self.wave_dist_combo.currentText()
        xpos_label = self.wave_xpos_combo.currentText()
        data = self._waveform_data.get(dist_label, {}).get(xpos_label)
        if data is None:
            return
        n = len(data["anode"])
        self.wave_event_slider.blockSignals(True)
        self.wave_event_slider.setMaximum(max(0, n - 1))
        self.wave_event_slider.setValue(0)
        self.wave_event_slider.blockSignals(False)
        self.wave_event_label.setText(f"1 / {n}")
        # Sync charge xpos (blocked) and update its slider
        self.charge_xpos_combo.blockSignals(True)
        self.charge_xpos_combo.setCurrentIndex(index)
        self.charge_xpos_combo.blockSignals(False)
        self._sync_charge_slider()
        self._update_waveform_plot()

    def _update_waveform_plot(self):
        """Draw the selected event in a ROOT TCanvas (anode top, cathode bottom)."""
        dist_label = self.wave_dist_combo.currentText()
        xpos_label = self.wave_xpos_combo.currentText()
        data  = self._waveform_data.get(dist_label, {}).get(xpos_label)
        label = f"{dist_label}  x={xpos_label}" if xpos_label != "—" else dist_label
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

            if not self._root_canvas_alive(self._root_canvas, "tgc_waveforms"):
                self._root_canvas = None
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
        dist_label = self.charge_dist_combo.currentText()
        xpos_dict  = self._waveform_data.get(dist_label, {})

        # Rebuild xpos combo for this distance
        self.charge_xpos_combo.blockSignals(True)
        self.charge_xpos_combo.clear()
        for xp in self._sorted_xpos(xpos_dict):
            self.charge_xpos_combo.addItem(xp)
        self.charge_xpos_combo.blockSignals(False)

        # Sync wave tab (blocked)
        self.wave_dist_combo.blockSignals(True)
        self.wave_dist_combo.setCurrentIndex(index)
        self.wave_dist_combo.blockSignals(False)
        self._rebuild_wave_xpos()

        self._on_charge_xpos_changed(0)

    def _on_charge_xpos_changed(self, index: int):
        dist_label = self.charge_dist_combo.currentText()
        xpos_label = self.charge_xpos_combo.currentText()
        data = self._waveform_data.get(dist_label, {}).get(xpos_label)
        if data is None:
            return
        n = len(data["anode"])
        self.charge_event_slider.blockSignals(True)
        self.charge_event_slider.setMaximum(max(0, n - 1))
        self.charge_event_slider.setValue(0)
        self.charge_event_slider.blockSignals(False)
        self.charge_event_label.setText(f"1 / {n}")
        # Sync wave xpos (blocked) and update its slider
        self.wave_xpos_combo.blockSignals(True)
        self.wave_xpos_combo.setCurrentIndex(index)
        self.wave_xpos_combo.blockSignals(False)
        self._sync_wave_slider()
        self._update_charge_plot()

    def _update_charge_plot(self):
        """Draw cumulative charge integrals Q(t) in a ROOT TCanvas (anode top, cathode bottom)."""
        dist_label = self.charge_dist_combo.currentText()
        xpos_label = self.charge_xpos_combo.currentText()
        data  = self._waveform_data.get(dist_label, {}).get(xpos_label)
        label = f"{dist_label}  x={xpos_label}" if xpos_label != "—" else dist_label
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

            if not self._root_canvas_alive(self._charge_canvas, "tgc_charges"):
                self._charge_canvas = None
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
        self._tracks_canvas  = None   # force recreation for the new run
        self._trk_zoom_scale = 1.0
        self._trk_view_phi   = 32.0
        self._trk_view_theta = 30.0
        self._trk_pan_x = 0.0
        self._trk_pan_y = 0.0
        self._trk_pan_z = 0.0
        self.trk_dist_combo.blockSignals(True)
        self.trk_dist_combo.clear()
        self.trk_xpos_combo.blockSignals(True)
        self.trk_xpos_combo.clear()

        # Geometry comes from run_config.json written by the binary
        cfg_path = Path(run_dir) / "run_config.json"
        if cfg_path.exists():
            try:
                with open(cfg_path) as f:
                    rc = json.load(f)
                g = rc.get("geometry", {})
                self._track_geom = {
                    "wire_pitch_cm": g.get("wire_pitch_cm",   0.18),
                    "n_wires":       int(g.get("n_wires",     10)),
                    "gap_cm":        g.get("gap_cm",          0.14),
                    "wire_diam_um":  g.get("wire_diameter_um", 50.0),
                    "wire_volt_V":   g.get("wire_voltage_V",  1900.0),
                }
            except Exception as exc:  # noqa: BLE001
                self.append_log(f"[GUI] run_config.json read error: {exc}")

        # Pre-populate E-Field tab geometry spinboxes
        if self._track_geom:
            tg = self._track_geom
            self.ef_gap.setValue(tg["gap_cm"])
            self.ef_pitch.setValue(tg["wire_pitch_cm"])
            self.ef_n_wires.setValue(tg["n_wires"])
            self.ef_wire_diam.setValue(tg["wire_diam_um"])
            self.ef_wire_volt.setValue(tg["wire_volt_V"])
            x_half = (tg["n_wires"] - 1) / 2 * tg["wire_pitch_cm"] + tg["wire_pitch_cm"]
            self.ef_y_depth.setRange(-tg["gap_cm"], tg["gap_cm"])
            self.ef_x_depth.setRange(-x_half, x_half)

        if not Path(root_path).exists():
            self.trk_dist_combo.blockSignals(False)
            return

        try:
            with uproot.open(root_path) as f:
                # Folder names: dist_0p1mm_x0p18mm  (or dist_0p1mm for old files)
                all_keys = sorted({
                    k.split("/")[0] for k in f.keys(cycle=False)
                    if k.split("/")[0].startswith("dist_")
                })
                for key in all_keys:
                    rest = key.removeprefix("dist_")
                    dist_raw, sep, xpos_raw = rest.partition("_x")
                    dist_label = "—" if dist_raw == "rnd" else dist_raw.replace("p", ".").replace("mm", " mm")
                    xpos_label = (xpos_raw.replace("p", ".").replace("mm", " mm")
                                  if sep else "—")
                    try:
                        tree = f[f"{key}/t_signals"]
                        if "primary_x" not in tree.keys():
                            continue   # pre-feature ROOT file — skip silently
                        data_dict = {
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
                        self._track_data.setdefault(dist_label, {})[xpos_label] = data_dict
                    except Exception as exc:  # noqa: BLE001
                        self.append_log(f"[GUI] Track load error for {key}: {exc}")
                # Populate dist combo in sorted order (random "—" first)
                for dl in self._sorted_dists(self._track_data.keys()):
                    self.trk_dist_combo.addItem(dl)
        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Could not open ROOT file for tracks: {exc}")

        self.trk_dist_combo.blockSignals(False)
        self.trk_xpos_combo.blockSignals(False)
        if self.trk_dist_combo.count():
            self._on_trk_dist_changed(0)
            self._root_timer.start()   # keep tracks window responsive

    def _on_trk_dist_changed(self, index: int):
        dist_label = self.trk_dist_combo.currentText()
        xpos_dict  = self._track_data.get(dist_label, {})

        self.trk_xpos_combo.blockSignals(True)
        self.trk_xpos_combo.clear()
        for xpos_label in sorted(
                xpos_dict.keys(),
                key=lambda s: float(s.replace(" mm", "")) if s != "—" else 0.0):
            self.trk_xpos_combo.addItem(xpos_label)
        self.trk_xpos_combo.blockSignals(False)

        self._on_trk_xpos_changed(0)

    def _on_trk_xpos_changed(self, index: int):
        dist_label = self.trk_dist_combo.currentText()
        xpos_label = self.trk_xpos_combo.currentText()
        data = self._track_data.get(dist_label, {}).get(xpos_label)
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
        """Render detector geometry and per-event tracks in a ROOT TCanvas."""
        dist_label = self.trk_dist_combo.currentText()
        xpos_label = self.trk_xpos_combo.currentText()
        data = self._track_data.get(dist_label, {}).get(xpos_label)
        if data is None:
            return
        label = f"{dist_label}  x={xpos_label}" if xpos_label != "—" else dist_label

        ev = self.trk_event_slider.value()
        n  = len(data["primary_x"])
        self.trk_event_label.setText(f"{ev + 1} / {n}")

        geom    = self._track_geom or {}
        pitch   = geom.get("wire_pitch_cm", 0.18)
        n_wires = int(geom.get("n_wires", 10))
        gap     = geom.get("gap_cm", 0.14)
        wire_diam_um = geom.get("wire_diam_um", 50.0)
        r_cm    = wire_diam_um / 2.0 * 1e-4   # µm → cm (actual physical radius)
        z_half  = 0.5
        x_half  = (n_wires - 1) / 2.0 * pitch + pitch
        x_wires = [(i - (n_wires - 1) / 2.0) * pitch for i in range(n_wires)]
        # largest physical half-range across all axes → equal TH3F ranges for
        # correct aspect ratios in ROOT's classic 3D renderer
        max_half = max(x_half, gap * 1.2, z_half)

        try:
            import ROOT  # noqa: PLC0415
            ROOT.gROOT.SetBatch(False)

            if not self._root_canvas_alive(self._tracks_canvas, "tgc_tracks"):
                self._tracks_canvas = None
            if self._tracks_canvas is None:
                self._tracks_canvas = ROOT.TCanvas(
                    "tgc_tracks", "TGC 3D Tracks", 900, 700)

            self._tracks_canvas.cd()
            self._tracks_canvas.Clear()
            self._tracks_objects.clear()
            self._tracks_canvas.SetPhi(self._trk_view_phi)
            self._tracks_canvas.SetTheta(self._trk_view_theta)

            # TH3F frame — defines axes and 3D coordinate range
            s  = self._trk_zoom_scale
            px = self._trk_pan_x
            py = self._trk_pan_y
            pz = self._trk_pan_z
            frame = ROOT.TH3F(
                "trk_frame",
                f"TGC 3D Tracks - {label}, event {ev + 1};"
                "x [cm];y [cm];z [cm]",
                1, px - max_half * s, px + max_half * s,
                1, py - max_half * s, py + max_half * s,
                1, pz - max_half * s, pz + max_half * s,
            )
            frame.SetStats(0)
            # Dynamic margins: make the inner pad area square so ROOT maps the
            # equal-range 3D cube without horizontal stretch, regardless of the
            # actual Qt-determined canvas pixel size.
            _cw = self._tracks_canvas.GetWw() or 900
            _ch = self._tracks_canvas.GetWh() or 700
            _top, _bottom = 0.05, 0.10
            _ih  = _ch * (1.0 - _top - _bottom)         # inner height in px
            _lr  = max(0.18, 1.0 - _ih / _cw)           # left+right fraction needed
            self._tracks_canvas.SetTopMargin(_top)
            self._tracks_canvas.SetBottomMargin(_bottom)
            self._tracks_canvas.SetLeftMargin(_lr * 0.55)   # more left for y-label
            self._tracks_canvas.SetRightMargin(_lr * 0.45)
            # Push all axis titles away from their axis lines.
            frame.GetXaxis().SetTitleOffset(1.6)
            frame.GetYaxis().SetTitleOffset(2.5)
            frame.GetZaxis().SetTitleOffset(1.6)
            frame.Draw()
            self._tracks_objects.append(frame)

            def _clip(xs, ys, zs):
                """Return contiguous sub-segments whose points lie within the
                visible cube [centre ± max_half*s].  Point-mask only — no exact
                boundary intersection, but avoids drawing far outside the frame."""
                hr = max_half * s
                cx_ = self._trk_pan_x
                cy_ = self._trk_pan_y
                cz_ = self._trk_pan_z
                mask = ((xs >= cx_ - hr) & (xs <= cx_ + hr) &
                        (ys >= cy_ - hr) & (ys <= cy_ + hr) &
                        (zs >= cz_ - hr) & (zs <= cz_ + hr))
                segs, i, n = [], 0, len(xs)
                while i < n:
                    if mask[i]:
                        j = i + 1
                        while j < n and mask[j]:
                            j += 1
                        if j - i >= 2:
                            segs.append((xs[i:j], ys[i:j], zs[i:j]))
                        i = j
                    else:
                        i += 1
                return segs

            def _pl3(xs, ys, zs, color, width=1, alpha=1.0):
                ln = ROOT.TPolyLine3D(
                    len(xs), xs.astype("f4"), ys.astype("f4"), zs.astype("f4"))
                if alpha < 1.0:
                    ln.SetLineColorAlpha(color, alpha)
                else:
                    ln.SetLineColor(color)
                ln.SetLineWidth(width)
                ln.Draw("SAME")
                self._tracks_objects.append(ln)

            # ── Cathode planes (rectangular outlines, clipped to visible cube) ──
            _hr = max_half * s
            for y_cath, color in [(-gap, ROOT.kCyan - 7), (gap, ROOT.kYellow - 7)]:
                if not (self._trk_pan_y - _hr <= y_cath <= self._trk_pan_y + _hr):
                    continue
                cx_min = max(-x_half, self._trk_pan_x - _hr)
                cx_max = min(+x_half, self._trk_pan_x + _hr)
                cz_min = max(-z_half, self._trk_pan_z - _hr)
                cz_max = min(+z_half, self._trk_pan_z + _hr)
                if cx_max <= cx_min or cz_max <= cz_min:
                    continue  # cathode not in view
                cx = np.array([cx_min, cx_max, cx_max, cx_min, cx_min], "f4")
                cy = np.full(5, y_cath, "f4")
                cz = np.array([cz_min, cz_min, cz_max, cz_max, cz_min], "f4")
                pl = ROOT.TPolyLine3D(5, cx, cy, cz)
                pl.SetLineColor(color)
                pl.SetLineWidth(2)
                pl.Draw("SAME")
                self._tracks_objects.append(pl)

            # ── Anode wires (cylinder wireframe, actual diameter, semi-transparent) ──
            _wire_alpha = 0.45
            _n_sides = 12                                     # 12-sided polygon approximation
            _angles = np.linspace(0.0, 2.0 * np.pi, _n_sides + 1)  # closed polygon
            # z extent of wires clipped to the visible cube
            _z_lo = max(-z_half, self._trk_pan_z - _hr)
            _z_hi = min(+z_half, self._trk_pan_z + _hr)
            for x_w in x_wires:
                if not (self._trk_pan_x - _hr <= x_w <= self._trk_pan_x + _hr):
                    continue
                xs_c = (x_w + r_cm * np.cos(_angles)).astype("f4")
                ys_c = (r_cm * np.sin(_angles)).astype("f4")
                if _z_hi > _z_lo:
                    # end rings at the visible z boundaries — always within the
                    # frame box, so the cross-section polygon is visible at any zoom
                    for z_ring in (_z_lo, _z_hi):
                        zs_c = np.full(_n_sides + 1, z_ring, "f4")
                        ring = ROOT.TPolyLine3D(_n_sides + 1, xs_c, ys_c, zs_c)
                        ring.SetLineColorAlpha(ROOT.kYellow + 1, _wire_alpha)
                        ring.SetLineWidth(1)
                        ring.Draw("SAME")
                        self._tracks_objects.append(ring)
                # longitudinal edges — z extent clipped to visible range
                if _z_hi > _z_lo:
                    for _a in _angles[:-1]:
                        _xe = float(x_w + r_cm * np.cos(_a))
                        _ye = float(r_cm * np.sin(_a))
                        ln = ROOT.TPolyLine3D(
                            2,
                            np.array([_xe, _xe], "f4"),
                            np.array([_ye, _ye], "f4"),
                            np.array([_z_lo, _z_hi], "f4"),
                        )
                        ln.SetLineColorAlpha(ROOT.kYellow + 1, _wire_alpha)
                        ln.SetLineWidth(1)
                        ln.Draw("SAME")
                        self._tracks_objects.append(ln)

            # ── Primary electron drift (blue) ─────────────────────────────────
            px = np.asarray(data["primary_x"][ev])
            py = np.asarray(data["primary_y"][ev])
            pz = np.asarray(data["primary_z"][ev])
            # No Python-level clipping — always draw the full track so it
            # remains visible at any zoom level. ROOT's 3D→2D projector handles
            # clipping at the pad boundary automatically.
            if len(px) >= 2:
                _pl3(px, py, pz, ROOT.kBlue + 1, 2, alpha=0.65)

            # ── Avalanche cloud (orange markers) ──────────────────────────────
            cx_ = np.asarray(data["cloud_x"][ev])
            cy_ = np.asarray(data["cloud_y"][ev])
            cz_ = np.asarray(data["cloud_z"][ev])
            if len(cx_) > 0:
                p = np.empty(3 * len(cx_), "f4")
                p[0::3] = cx_; p[1::3] = cy_; p[2::3] = cz_
                mrk = ROOT.TPolyMarker3D(len(cx_), p, 7)
                mrk.SetMarkerColorAlpha(ROOT.kOrange + 1, 0.50)
                mrk.SetMarkerSize(0.4)
                mrk.Draw("SAME")
                self._tracks_objects.append(mrk)

            # ── Ion drift paths (colour-coded by destination) ─────────────────
            ion_x = np.asarray(data["ion_x"][ev])
            ion_y = np.asarray(data["ion_y"][ev])
            ion_z = np.asarray(data["ion_z"][ev])
            inpts = np.asarray(data["ion_npts"][ev])
            tol   = 0.05 * gap
            off   = 0
            for n_seg in inpts:
                xs = ion_x[off:off + n_seg]
                ys = ion_y[off:off + n_seg]
                zs = ion_z[off:off + n_seg]
                off += n_seg
                if len(ys) < 2:
                    # single start-point: ion immediately absorbed — draw as marker
                    if len(ys) == 1:
                        pt = np.array([float(xs[0]), float(ys[0]), float(zs[0])], "f4")
                        dot = ROOT.TPolyMarker3D(1, pt, 7)
                        dot.SetMarkerColorAlpha(ROOT.kGray + 2, 0.5)
                        dot.SetMarkerSize(0.3)
                        dot.Draw("SAME")
                        self._tracks_objects.append(dot)
                    continue
                y_end = float(ys[-1])
                if abs(y_end - (-gap)) <= tol:
                    col = ROOT.kGreen + 2    # → readout cathode
                elif abs(y_end - gap) <= tol:
                    col = ROOT.kMagenta      # → non-readout cathode
                else:
                    col = ROOT.kGray + 1     # absorbed / out of window
                for _seg in _clip(xs, ys, zs):
                    _pl3(*_seg, col, 1, alpha=0.55)

            # ── Legend ────────────────────────────────────────────────────────
            self._trk_legend_objects.clear()

            def _mk_line(color: int, width: int = 1) -> object:
                ln = ROOT.TLine()
                ln.SetLineColor(color)
                ln.SetLineWidth(width)
                return ln

            def _mk_marker(color: int) -> object:
                mk = ROOT.TMarker()
                mk.SetMarkerColor(color)
                mk.SetMarkerStyle(7)
                mk.SetMarkerSize(1.2)
                return mk

            _leg_proxies = [
                (_mk_line(ROOT.kBlue + 1, 2),  "Primary e^{-}",          "L"),
                (_mk_marker(ROOT.kOrange + 1), "Avalanche",               "P"),
                (_mk_line(ROOT.kGreen + 2),    "Ion #rightarrow readout", "L"),
                (_mk_line(ROOT.kMagenta),      "Ion #rightarrow other",   "L"),
                (_mk_line(ROOT.kGray + 1),     "Ion (absorbed)",          "L"),
            ]
            leg = ROOT.TLegend(0.70, 0.76, 0.99, 0.97)
            leg.SetBorderSize(0)
            leg.SetFillColorAlpha(ROOT.kWhite, 0.75)
            leg.SetTextSize(0.028)
            for _proxy, _label, _opt in _leg_proxies:
                leg.AddEntry(_proxy, _label, _opt)
            leg.Draw()
            self._trk_legend_objects = [leg] + [p for p, _, _ in _leg_proxies]
            # ──────────────────────────────────────────────────────────────────

            self._tracks_canvas.Update()

        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] ROOT 3D tracks error: {exc}")

    def _trk_adjust_zoom(self, factor: float) -> None:
        """Scale the visible axis range and redraw (factor < 1 = zoom in)."""
        self._trk_zoom_scale = max(0.005, min(20.0, self._trk_zoom_scale * factor))
        self._update_track_plot()

    def _trk_pan(self, axis: str, direction: int) -> None:
        """Shift the visible centre by 30 % of the current visible half-range."""
        geom    = self._track_geom or {}
        pitch   = geom.get("wire_pitch_cm", 0.18)
        n_wires = int(geom.get("n_wires", 10))
        gap     = geom.get("gap_cm", 0.14)
        x_half  = (n_wires - 1) / 2.0 * pitch + pitch
        max_half = max(x_half, gap * 1.2, 0.5)   # matches equal-range TH3F
        step = max_half * self._trk_zoom_scale * 0.3 * direction
        if   axis == "x": self._trk_pan_x += step
        elif axis == "y": self._trk_pan_y += step
        else:             self._trk_pan_z += step
        self._update_track_plot()

    def _trk_preset_view(self, phi: float, theta: float,
                         reset_zoom: bool = False) -> None:
        """Set pad view angles and redraw."""
        self._trk_view_phi   = phi
        self._trk_view_theta = theta
        if reset_zoom:
            self._trk_zoom_scale = 1.0
            self._trk_pan_x = self._trk_pan_y = self._trk_pan_z = 0.0
        self._update_track_plot()

    # ── E-Field ────────────────────────────────────────────────────────────────

    @staticmethod
    def _compute_efield(gap, pitch, n_wires, wire_diam_um, wire_volt,
                        nx=250, ny=150, n_images=12):
        """
        Compute TGC electric field on a 2D grid (XY plane) using the image-charge
        method with Sauli's capacitance approximation.

        Returns (x, y, Ex, Ey) — all in cm / (V/cm).
        """
        a       = wire_diam_um * 1e-4 / 2           # wire radius [cm]
        x_wires = [(i - (n_wires - 1) / 2) * pitch for i in range(n_wires)]
        x_half  = (n_wires - 1) / 2 * pitch + pitch
        eps0    = 8.854e-14                          # F/cm

        # Sauli 1977 capacitance coefficient (eq. 2.3)
        C_inv  = np.log(pitch / (2 * np.pi * a)) + np.pi * gap / pitch
        lam_w  = 2 * np.pi * eps0 * wire_volt / C_inv  # C/cm

        x = np.linspace(-x_half, x_half, nx)
        y = np.linspace(-gap * 0.999, gap * 0.999, ny)
        X, Y = np.meshgrid(x, y)

        Ex = np.zeros((ny, nx))
        Ey = np.zeros((ny, nx))
        n_arr = np.arange(-n_images, n_images + 1)

        for x_w in x_wires:
            y_imgs = (2 * n_arr * gap)[:, None, None]          # (N, 1, 1)
            signs  = ((-1) ** np.abs(n_arr))[:, None, None]    # (N, 1, 1)
            dx     = X[None] - x_w                              # (1, ny, nx)
            dy     = Y[None] - y_imgs                           # (N, ny, nx)
            r2     = np.maximum(dx ** 2 + dy ** 2, a ** 2)
            fac    = signs * lam_w / (2 * np.pi * eps0 * r2)
            Ex    += np.sum(fac * dx, axis=0)
            Ey    += np.sum(fac * dy, axis=0)

        return x, y, Ex, Ey

    def _update_efield_plots(self, recompute: bool = True):
        """Draw electric-field maps for XY, ZX-slice, and ZY-slice planes."""
        if recompute:
            gap       = self.ef_gap.value()
            pitch     = self.ef_pitch.value()
            n_wires   = self.ef_n_wires.value()
            diam_um   = self.ef_wire_diam.value()
            wire_volt = self.ef_wire_volt.value()
            try:
                x, y, Ex, Ey = self._compute_efield(
                    gap, pitch, n_wires, diam_um, wire_volt,
                    nx=self.ef_nx.value(), ny=self.ef_ny.value())
            except Exception as exc:  # noqa: BLE001
                self.append_log(f"[GUI] E-field computation error: {exc}")
                return
            self._efield_cache = {"x": x, "y": y, "Ex": Ex, "Ey": Ey,
                                  "gap": gap, "pitch": pitch,
                                  "n_wires": n_wires}
        elif self._efield_cache is None:
            return  # nothing cached yet

        c       = self._efield_cache
        x, y    = c["x"], c["y"]
        Ex, Ey  = c["Ex"], c["Ey"]
        gap     = c["gap"]
        pitch   = c["pitch"]
        n_wires = c["n_wires"]
        x_wires = [(i - (n_wires - 1) / 2) * pitch for i in range(n_wires)]
        x_half  = (n_wires - 1) / 2 * pitch + pitch

        comp = self.ef_component.currentText()
        if   comp == "|E|": field = np.sqrt(Ex ** 2 + Ey ** 2)
        elif comp == "Ex":  field = Ex
        else:               field = Ey

        y0 = self.ef_y_depth.value()
        x0 = self.ef_x_depth.value()
        iy = int(np.argmin(np.abs(y - y0)))
        ix = int(np.argmin(np.abs(x - x0)))
        X, Y = np.meshgrid(x, y)

        try:
            import ROOT  # noqa: PLC0415
            ROOT.gROOT.SetBatch(False)

            if not self._root_canvas_alive(self._efield_root_canvas, "tgc_efield"):
                self._efield_root_canvas = None
            if self._efield_root_canvas is None:
                self._efield_root_canvas = ROOT.TCanvas(
                    "tgc_efield", "TGC E-Field", 1300, 500)

            self._efield_root_canvas.Clear()
            self._efield_root_canvas.Divide(3, 1)
            self._efield_objects.clear()

            # ── Pad 1: XY plane (TH2F COLZ) ──────────────────────────────────
            self._efield_root_canvas.cd(1)
            ROOT.gPad.SetLeftMargin(0.12)
            ROOT.gPad.SetRightMargin(0.17)   # extra room for COLZ colour bar
            ROOT.gPad.SetTopMargin(0.12)
            ROOT.gPad.SetBottomMargin(0.14)
            nx_v, ny_v = len(x), len(y)
            dx = (x[-1] - x[0]) / max(nx_v - 1, 1)
            dy = (y[-1] - y[0]) / max(ny_v - 1, 1)
            h2 = ROOT.TH2F(
                "h_efield",
                f"XY plane - {comp} [V/cm];x [cm];y [cm]",
                nx_v, x[0] - dx / 2, x[-1] + dx / 2,
                ny_v, y[0] - dy / 2, y[-1] + dy / 2,
            )
            h2.SetStats(0)
            h2.FillN(
                nx_v * ny_v,
                X.flatten().astype("f8"),
                Y.flatten().astype("f8"),
                field.flatten().astype("f8"),
            )
            h2.Draw("COLZ")
            for x_w in x_wires:
                m = ROOT.TMarker(x_w, 0.0, 5)
                m.SetMarkerColor(ROOT.kWhite)
                m.SetMarkerSize(1.5)
                m.Draw()
                self._efield_objects.append(m)
            for y_c in [-gap, gap]:
                ln = ROOT.TLine(-x_half, y_c, x_half, y_c)
                ln.SetLineStyle(2)
                ln.SetLineColor(ROOT.kGray + 1)
                ln.Draw()
                self._efield_objects.append(ln)
            self._efield_objects.append(h2)

            # ── Pad 2: E(x) profile at y = y0 ────────────────────────────────
            self._efield_root_canvas.cd(2)
            ROOT.gPad.SetLeftMargin(0.16)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetTopMargin(0.12)
            ROOT.gPad.SetBottomMargin(0.14)
            g_xz = ROOT.TGraph(len(x), x.astype("f8"), field[iy, :].astype("f8"))
            g_xz.SetTitle(
                f"{comp} at y = {y[iy]:.3f} cm;x [cm];{comp} [V/cm]")
            g_xz.SetLineColor(ROOT.kBlue + 1)
            g_xz.SetLineWidth(2)
            ROOT.gPad.SetGrid()
            g_xz.Draw("AL")
            self._efield_objects.append(g_xz)

            # ── Pad 3: E(y) profile at x = x0 ────────────────────────────────
            self._efield_root_canvas.cd(3)
            ROOT.gPad.SetLeftMargin(0.16)
            ROOT.gPad.SetRightMargin(0.05)
            ROOT.gPad.SetTopMargin(0.12)
            ROOT.gPad.SetBottomMargin(0.14)
            g_zy = ROOT.TGraph(len(y), y.astype("f8"), field[:, ix].astype("f8"))
            g_zy.SetTitle(
                f"{comp} at x = {x[ix]:.3f} cm;y [cm];{comp} [V/cm]")
            g_zy.SetLineColor(ROOT.kRed + 1)
            g_zy.SetLineWidth(2)
            ROOT.gPad.SetGrid()
            g_zy.Draw("AL")
            self._efield_objects.append(g_zy)

            self._efield_root_canvas.Update()
            self._root_timer.start()   # ensure timer is running

        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] ROOT E-field error: {exc}")


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

        # Build the run subfolder name: yymmdd_hh-mm__ + user tag or auto params
        date_pfx = datetime.now().strftime("%y%m%d_%H-%M")
        tag      = self.config_panel.run_name.text().strip()
        if tag:
            subdir = f"{date_pfx}__{tag}"
        else:
            v = int(cfg["geometry"]["wire_voltage_V"])
            n = cfg["simulation"]["n_events"]
            subdir = f"{date_pfx}__V{v}V__n{n}"

        self.results_panel.clear_log()
        self.results_panel.setCurrentIndex(0)   # show Log tab while running
        self.statusBar().showMessage("Running…")

        self._runner = SimRunner(cfg, str(out_path), run_name=subdir)
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
        self._try_load_gas_props()
        self.results_panel.setCurrentIndex(1)   # switch to Summary tab
        self.results_panel._save_plots_root(run_dir)

    def _try_load_gas_props(self):
        """Load Magboltz properties CSV if it exists for the current gas config."""
        gas_cfg = self.config_panel.to_config_dict().get("gas", {})
        props_path = TGC_DIR / derive_gas_props_filename(gas_cfg)
        if props_path.exists():
            self.results_panel.draw_gas_props(str(props_path))

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
            self._try_load_gas_props()
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

        # Stop ROOT timer and close all ROOT TCanvas windows before Qt tears
        # down its macOS Cocoa layer — prevents "drawable not found" crash.
        rp = self.results_panel
        rp._root_timer.stop()
        try:
            import ROOT  # noqa: PLC0415
            for _canvas in [rp._root_canvas, rp._charge_canvas,
                            rp._tracks_canvas, rp._efield_root_canvas,
                            rp._gas_canvas]:
                try:
                    if (_canvas is not None and
                            ROOT.gROOT.GetListOfCanvases()
                                      .FindObject(_canvas.GetName())):
                        _canvas.Close()
                except Exception:  # noqa: BLE001
                    pass
        except Exception:  # noqa: BLE001
            pass

        super().closeEvent(event)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("TGC Simulation")
    win = MainWindow()
    win.show()
    code = app.exec_()
    # os._exit bypasses Python/ROOT atexit destructors that crash on macOS
    # when ROOT's Cocoa layer outlives Qt's autorelease pool.
    os._exit(code)


if __name__ == "__main__":
    main()
