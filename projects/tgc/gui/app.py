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
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QFontDatabase
from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
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

SCRIPT_DIR = Path(__file__).parent.resolve()   # …/projects/tgc/gui/
TGC_DIR    = (SCRIPT_DIR / "..").resolve()     # …/projects/tgc/
BINARY     = TGC_DIR / "build" / "tgc_sim"


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

        try:
            self._proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                cwd=str(TGC_DIR),   # relative gas-file paths resolve from here
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

        gas_row  = QWidget()
        gas_h    = QHBoxLayout(gas_row)
        gas_h.setContentsMargins(0, 0, 0, 0)
        self.gas_file = QLineEdit("ar_70_co2_30.gas")
        self.gas_file.setToolTip("Path to the Magboltz gas table file")
        btn_gas = QPushButton("…")
        btn_gas.setFixedWidth(28)
        btn_gas.clicked.connect(self._browse_gas_file)
        gas_h.addWidget(self.gas_file)
        gas_h.addWidget(btn_gas)

        self.penning = QCheckBox()
        self.penning.setChecked(True)
        self.ncoll = self._spin(1, 100, 10)
        self.ncoll.setToolTip("Magboltz collision cycles per field point (higher = more accurate)")

        gas_form.addRow("Temperature [K]",     self.temperature)
        gas_form.addRow("Pressure [Torr]",     self.pressure)
        gas_form.addRow("Gas file",            gas_row)
        gas_form.addRow("Penning transfer",    self.penning)
        gas_form.addRow("Magboltz ncoll",      self.ncoll)
        root_layout.addWidget(gas_box)

        # ── Simulation ────────────────────────────────────────────────────
        sim_box  = QGroupBox("Simulation")
        sim_form = QFormLayout(sim_box)

        self.n_events    = self._spin(1, 100000, 1000)
        self.max_aval    = self._spin(1000, 10000000, 500000)
        self.time_window = self._dspin(10.0, 2000.0, 10.0, 1, 300.0)
        self.time_step   = self._dspin(0.1, 10.0, 0.1, 2, 0.5)

        sim_form.addRow("Events",             self.n_events)
        sim_form.addRow("Max avalanche size", self.max_aval)
        sim_form.addRow("Time window [ns]",   self.time_window)
        sim_form.addRow("Time step [ns]",     self.time_step)
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

    # ── file dialogs ─────────────────────────────────────────────────────

    def _browse_gas_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select gas file", "",
            "Gas files (*.gas);;All files (*)"
        )
        if path:
            self.gas_file.setText(path)

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

        return {
            "geometry": {
                "wire_pitch_cm":    self.wire_pitch.value(),
                "wire_diameter_um": self.wire_diam.value(),
                "gap_cm":           self.gap_cm.value(),
                "n_wires":          self.n_wires.value(),
                "wire_voltage_V":   self.wire_volts.value(),
            },
            "source": {
                "energy_keV":          self.energy_kev.value(),
                "source_distances_mm": dists,
                "x_position_cm":       x_pos,
            },
            "gas": {
                "temperature_K":         self.temperature.value(),
                "pressure_Torr":         self.pressure.value(),
                "gas_file":              self.gas_file.text(),
                "enable_penning":        self.penning.isChecked(),
                "n_magboltz_collisions": self.ncoll.value(),
            },
            "simulation": {
                "n_events":           self.n_events.value(),
                "max_avalanche_size": self.max_aval.value(),
                "time_window_ns":     self.time_window.value(),
                "time_step_ns":       self.time_step.value(),
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
        self.gas_file.setText(    gas.get("gas_file", "ar_70_co2_30.gas"))
        self.penning.setChecked(  gas.get("enable_penning", True))
        self.ncoll.setValue(      gas.get("n_magboltz_collisions", 10))

        sim = d.get("simulation", {})
        self.n_events.setValue(   sim.get("n_events", 1000))
        self.max_aval.setValue(   sim.get("max_avalanche_size", 500000))
        self.time_window.setValue(sim.get("time_window_ns", 300.0))
        self.time_step.setValue(  sim.get("time_step_ns", 0.5))


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
        self.table = QTableWidget()
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.addTab(self.table, "Summary")

        # ── Plots tab: 2×2 matplotlib figure ─────────────────────────────
        self.plots_canvas = MplCanvas(nrows=2, ncols=2, figsize=(7, 5))
        self.addTab(self.plots_canvas, "Plots")

        # ── Waveforms tab: 1×2 matplotlib figure ─────────────────────────
        self.wave_canvas = MplCanvas(nrows=1, ncols=2, figsize=(7, 3))
        self.addTab(self.wave_canvas, "Waveforms")

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

        _errplot(axes[0], "mean_anode_charge_fC",   "sem_anode_charge_fC",
                 "⟨Q_anode⟩ [fC]",        "Anode charge vs distance")
        _errplot(axes[1], "mean_cathode_charge_fC", "sem_cathode_charge_fC",
                 "⟨Q_cathode⟩ [fC]",      "Cathode charge vs distance")
        _errplot(axes[2], "mean_charge_ratio",      "sem_charge_ratio",
                 "Q_cathode / Q_anode",   "Charge ratio vs distance")
        _errplot(axes[3], "mean_avalanche_size",    None,
                 "⟨Avalanche size⟩",      "Avalanche size vs distance")

        self.plots_canvas.figure.tight_layout()
        self.plots_canvas.draw()

    # ── Waveforms ─────────────────────────────────────────────────────────

    def draw_waveforms(self, root_path: str):
        try:
            import uproot  # noqa: PLC0415 — optional import
        except ImportError:
            self.append_log("[GUI] uproot not available — waveform tab will be empty")
            return

        if not Path(root_path).exists():
            self.append_log(f"[GUI] ROOT file not found: {root_path}")
            return

        ax_a, ax_c = self.wave_canvas.axes
        ax_a.cla()
        ax_c.cla()

        try:
            with uproot.open(root_path) as f:
                # Collect unique top-level dist_* directory names
                dist_keys = []
                seen: set[str] = set()
                for raw_key in f.keys(cycle=False):
                    top = raw_key.split("/")[0]
                    if top.startswith("dist_") and top not in seen:
                        dist_keys.append(top)
                        seen.add(top)

                for key in sorted(dist_keys):
                    # Build a human-readable label from "dist_0p7mm" → "0.7 mm"
                    label = key.removeprefix("dist_").replace("p", ".").replace("mm", " mm")

                    try:
                        pa = f[f"{key}/p_anode_signal"]
                        pc = f[f"{key}/p_cathode_signal"]
                        times   = pa.axis(0).centers()
                        anode   = pa.values()
                        cathode = pc.values()
                    except (KeyError, AttributeError, ValueError) as exc:
                        self.append_log(f"[GUI] Could not read waveforms for {key}: {exc}")
                        continue

                    ax_a.plot(times, anode,   lw=1.5, label=label)
                    ax_c.plot(times, cathode, lw=1.5, label=label)

        except Exception as exc:  # noqa: BLE001
            self.append_log(f"[GUI] Could not read ROOT file: {exc}")
            return

        for ax, title in [(ax_a, "Mean anode signal"), (ax_c, "Mean cathode signal")]:
            ax.set_xlabel("Time [ns]")
            ax.set_ylabel("Current [fC/ns]")
            ax.set_title(title)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        self.wave_canvas.figure.tight_layout()
        self.wave_canvas.draw()


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._runner: SimRunner | None = None

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
        self.results_panel.draw_waveforms(root_path)
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
