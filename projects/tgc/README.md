# TGC Detector Simulation

A Garfield++ simulation of a Thin Gap Chamber (TGC) multi-wire proportional detector.
The simulation models primary ionisation from 5.9 keV Fe-55 X-rays, electron drift and
avalanche multiplication, and the induced charge on both the wire plane (anode) and one
cathode plane.  Results can be explored via a PyQt5 desktop GUI or the command-line
binary directly.

---

## Project structure

```
projects/tgc/
├── src/
│   └── tgc_sim.cc              ← simulation binary (C++20, ~500 lines)
├── config/
│   ├── default_tgc.json        ← production config (4 distances, 1000 events each)
│   └── smoke_tgc.json          ← fast smoke-test config (1 distance, 10 events)
├── gui/
│   └── app.py                  ← PyQt5 desktop GUI
├── third_party/
│   └── nlohmann/json.hpp       ← vendored single-header JSON library
├── CMakeLists.txt
├── build/                      ← cmake output (gitignored)
└── results/                    ← simulation output (gitignored)
```

---

## Detector geometry

A TGC consists of an array of thin anode wires stretched between two grounded cathode
planes.  The small cathode-anode gap (here 1.4 mm) produces a steep electric-field
gradient that enables high gas gain at modest applied voltages.

```
  y = +1.4 mm  ───────────────────────────────  cathode_top (0 V, ground)
               ||||||||    gas gap    |||||||||
  y =  0.0 mm  ─ ○ ─ ○ ─ ○ ─ ○ ─ ○ ─ ○ ─ ○ ─  anode wires (+1900 V)
               ||||||||    gas gap    |||||||||
  y = -1.4 mm  ───────────────────────────────  cathode (0 V, readout)

               ← 1.8 mm →
               wire pitch
```

| Parameter         | Value        | Notes                               |
|-------------------|--------------|-------------------------------------|
| Wire count        | 10           |                                     |
| Wire diameter     | 50 μm        | radius = 25 μm                      |
| Wire pitch        | 1.8 mm       | centre-to-centre spacing            |
| Cathode-anode gap | 1.4 mm       | distance from wire plane to cathode |
| Wire voltage      | +1900 V      | configurable via `wire_voltage_V`   |
| Cathode voltage   | 0 V          | both planes grounded                |
| Gas               | Ar:CO2 70:30 | 1 atm, 20 °C                        |

---

## Physics

### 1. Gas transport coefficients — Magboltz

Electron drift velocity, diffusion, attachment, and Townsend coefficients are computed
by the Magboltz Monte Carlo code (via `MediumMagboltz`) over a logarithmic electric-field
grid from 100 V/cm to 300 kV/cm.  Results are cached to a `.gas` file and reloaded on
subsequent runs.

**Penning transfer** is enabled by default.  In Ar:CO2 the lowest Ar metastable levels
(Ar\*(³P₀) and Ar\*(³P₂) at 11.55–11.72 eV) lie *below* the CO2 ionisation potential
(13.78 eV) and cannot directly ionise CO2 molecules.  The effective Penning enhancement
observed in practice arises from higher Ar excited states (3p⁵4p and above) whose
energies reach and exceed 13.78 eV, allowing them to ionise CO2 via

> Ar\* + CO2 → Ar + CO2⁺ + e⁻

Garfield++ models this with an effective transfer fraction *r* (tabulated from measured
absolute gas-gain data for Ar:CO2 70:30) applied inside
`MediumMagboltz::EnablePenningTransfer()`.

**Ion mobility** — after avalanche multiplication, positive ions must drift back to the
cathode to complete the Ramo-theorem induced-charge calculation.  In Ar:CO2 the dominant
drifting species is CO2⁺ rather than Ar⁺: the lower ionisation potential of CO2
(13.78 eV vs. Ar 15.76 eV) means Ar⁺ rapidly charge-transfers to CO2⁺ on nanosecond
timescales.  The simulation loads the `IonMobility_CO2+_CO2.txt` table from the
Garfield++ data directory.  Garfield++ stores one positive-ion mobility table per gas
object; CO2⁺ is the best single-species approximation for this mixture.

### 2. Primary ionisation — TrackHeed

The 5.9 keV photon is transported by `TrackHeed::TransportPhoton`.  Heed uses detailed
cross-section tables to simulate:

* **Photoelectric absorption** in the gas (dominant at 5.9 keV in Ar:CO2).  The mean
  free path at 1 atm is several centimetres, so only a small fraction of photons interact
  in the 1.4 mm gap; the simulation skips non-interacting events and reports the
  interaction fraction in the summary.
* **Delta-electron cascade** — the photoelectron slows down and produces secondary
  ionisation along its track.  With W ≈ 26 eV/pair for Ar:CO2, a 5.9 keV photon yields
  roughly 220 primary electron-ion pairs.  All conduction electrons are returned with
  position and kinetic energy, ready for `AvalancheMicroscopic`.
* **Fluorescence photons** (Ar K-alpha, ~2.96 keV) may appear alongside the primary
  cluster.  They are currently not re-transported; adding a second `TransportPhoton` call
  for each fluorescence photon is a straightforward extension.

### 3. Electron avalanche — AvalancheMicroscopic

Each primary electron is transported individually by `AvalancheMicroscopic`, which steps
electrons through the gas using the Runge-Kutta-Fehlberg algorithm and samples elastic,
inelastic, ionising, and attachment collisions from the Magboltz cross-section tables.
All avalanche electrons are tracked until they are collected by a wire or absorbed.

A `max_avalanche_size` cap prevents runaway events from consuming excessive CPU.

### 4. Signal induction — Shockley-Ramo theorem

`Sensor` computes the induced current on each electrode at every step using the weighting
(Ramo) field:

> i(t) = q · v(t) · **E**_w(**x**(t))

where **E**_w is the weighting field of the electrode (the field that would exist if that
electrode were at 1 V and all others at 0 V, with all space charges removed).

Two readout channels are defined:

* **`anode`** — all 10 wires share this label; their weighting fields are summed
  automatically by `Sensor`.
* **`cathode`** — the bottom cathode plane at y = −1.4 mm.

---

## Gas file

The first run generates `ar_70_co2_30.gas` in the working directory.  Magboltz runs
`n_magboltz_collisions` collision cycles per field point (default: 10; smoke test: 5).
This takes roughly **5–15 minutes** for the default settings.

To skip regeneration on subsequent runs, place the `.gas` file at the path given by
`gas.gas_file` in the config (default: `ar_70_co2_30.gas` in the current working
directory).  You can move or rename it and update the config key accordingly.

To regenerate with higher accuracy, increase `n_magboltz_collisions` (20–50 is
recommended for publication-quality results) and delete the existing `.gas` file.

---

## Build

**Requirements:** Garfield++ installed at `../../local/garfield`, ROOT 6 (via conda),
CMake ≥ 3.20, a C++20-capable compiler.

From the `projects/tgc/` directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DVDT_INCLUDE_DIR=<conda-root>/include \
  -DVDT_LIBRARY=<conda-root>/lib/libvdt.dylib \
  -DCMAKE_PREFIX_PATH="$(pwd)/../../local/garfield;<conda-root>"
cmake --build build -j4
```

Replace `<conda-root>` with your conda environment root
(e.g. `$(conda info --base)` or the output of `conda info --json | python3 -m json.tool | grep active_prefix`).

**Why the explicit flags are needed:**  `source ../../local/garfield/share/Garfield/setupGarfield.sh`
sets `GARFIELD_INSTALL` and `HEED_DATABASE` (needed at *runtime*) but CMake
find-package resolution requires an explicit `CMAKE_PREFIX_PATH`.  Additionally,
ROOT's `FindVdt.cmake` does not search conda paths automatically, so `VDT_INCLUDE_DIR`
and `VDT_LIBRARY` must be supplied explicitly.

**Concrete invocation for this machine** (miniforge3 at `/Users/luca/miniforge3`):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DVDT_INCLUDE_DIR=/Users/luca/miniforge3/include \
  -DVDT_LIBRARY=/Users/luca/miniforge3/lib/libvdt.dylib \
  -DCMAKE_PREFIX_PATH="$(pwd)/../../local/garfield;/Users/luca/miniforge3"
cmake --build build -j4
```

---

## Usage

### Command-line binary

```bash
# Run with default config (4 source distances, 1000 events each):
./build/tgc_sim

# Custom config:
./build/tgc_sim --config config/default_tgc.json --out results/

# Single source distance (quick check):
./build/tgc_sim --distance 0.7

# All options:
./build/tgc_sim --help
```

The smoke CTest validates the build without generating a gas file from scratch
(assumes `ar_70_co2_30.gas` already exists):

```bash
ctest --test-dir build --verbose
```

---

## GUI

`gui/app.py` is a PyQt5 desktop application that wraps `tgc_sim` with a
point-and-click interface: edit parameters, launch, and inspect results — all
without touching JSON or the terminal.

### Prerequisites

```bash
conda install pyqt matplotlib pandas    # if not already present
pip install uproot                       # ROOT-file reading (pure Python, no ROOT install needed)
```

### Launch

```bash
# From the repo root or any directory:
python3 projects/tgc/gui/app.py
```

The binary is resolved automatically as `projects/tgc/build/tgc_sim`.  If it has not
been built yet, the window opens with a warning in the title bar.

### Window layout

```
┌──────────────────────────────────────────────────────────────────────────┐
│  TGC Simulation               [▶ Run]  [■ Stop]  [Load Config]  [Save Config] │
├────────────────────────────────┬─────────────────────────────────────────┤
│  ▼ Geometry                    │  [ Log | Summary | Plots | Waveforms ]  │
│    Wire pitch [cm]   0.18      │                                         │
│    Wire diameter [μm] 50       │  ← live output / results shown here     │
│    Gap [cm]          0.14      │                                         │
│    N wires           10        │                                         │
│    Wire voltage [V]  1900      │                                         │
│  ▼ Source                      │                                         │
│    Energy [keV]      5.9       │                                         │
│    Distances [mm]  0.2,0.5,…  │                                         │
│    X position  [✓ Random]      │                                         │
│  ▼ Gas                         │                                         │
│    Temperature [K]   293.15    │                                         │
│    Pressure [Torr]   760       │                                         │
│    Gas file  [ar_70…gas] […]  │                                         │
│    Penning  [✓]  ncoll  10     │                                         │
│  ▼ Simulation                  │                                         │
│    Events         1000         │                                         │
│    Max aval. size 500000       │                                         │
│    Time window [ns] 300        │                                         │
│    Time step [ns]   0.5        │                                         │
│  ▼ Output                      │                                         │
│    Directory  [results/] […]   │                                         │
└────────────────────────────────┴─────────────────────────────────────────┘
```

| Tab | Contents |
|---|---|
| **Log** | Live stdout stream from `tgc_sim`, auto-scrolling |
| **Summary** | Table from `summary.csv` — one row per source distance |
| **Plots** | 2 × 2 matplotlib figure: ⟨Q_anode⟩, ⟨Q_cathode⟩, charge ratio, and avalanche size vs source distance (with SEM error bars) |
| **Waveforms** | Mean anode and cathode current waveforms overlaid per distance, read directly from the ROOT file via uproot |

**▶ Run** starts the simulation in a background thread (the window stays fully
responsive).  **■ Stop** sends SIGTERM.  **Load Config** / **Save Config** read and
write `.json` files that are fully compatible with the CLI `--config` flag.

---

## Configuration reference

All parameters live in a JSON file (default: `config/default_tgc.json`).

### `geometry`

| Key                | Type  | Unit | Default | Description                         |
|--------------------|-------|------|---------|-------------------------------------|
| `wire_pitch_cm`    | float | cm   | 0.18    | Centre-to-centre wire spacing       |
| `wire_diameter_um` | float | μm   | 50.0    | Wire outer diameter                 |
| `gap_cm`           | float | cm   | 0.14    | Distance from wire plane to cathode |
| `n_wires`          | int   | —    | 10      | Number of anode wires               |
| `wire_voltage_V`   | float | V    | 1900.0  | Voltage applied to all wires        |

### `source`

| Key                   | Type          | Unit | Default              | Description                                               |
|-----------------------|---------------|------|----------------------|-----------------------------------------------------------|
| `energy_keV`          | float         | keV  | 5.9                  | Photon energy (Fe-55 K-alpha line)                        |
| `source_distances_mm` | float array   | mm   | [0.2, 0.5, 0.9, 1.2] | Source y-distances from wire plane to scan                |
| `x_position_cm`       | float or null | cm   | null                 | Fixed photon x-position; `null` = uniform random over wires |

`source_distances_mm` is the scan axis.  Each value is the y-coordinate (in mm) at which
the photon is placed, measured from the wire plane (y = 0).  Values must lie in the range
(0, `gap_cm` × 10).  The photon always travels in the −y direction.

### `gas`

| Key                     | Type   | Unit | Default            | Description                                      |
|-------------------------|--------|------|--------------------|--------------------------------------------------|
| `temperature_K`         | float  | K    | 293.15             | Gas temperature                                  |
| `pressure_Torr`         | float  | Torr | 760.0              | Gas pressure                                     |
| `gas_file`              | string | —    | "ar_70_co2_30.gas" | Path to gas table (generated if absent)          |
| `enable_penning`        | bool   | —    | true               | Enable Penning transfer (recommended for Ar:CO2) |
| `n_magboltz_collisions` | int    | —    | 10                 | Magboltz collision cycles per field point        |

### `simulation`

| Key                  | Type  | Unit | Default  | Description                                     |
|----------------------|-------|------|----------|-------------------------------------------------|
| `n_events`           | int   | —    | 1000     | Number of photon events per source distance     |
| `max_avalanche_size` | int   | —    | 500000   | Electron count cap per `AvalancheElectron` call |
| `time_window_ns`     | float | ns   | 300.0    | Signal collection window                        |
| `time_step_ns`       | float | ns   | 0.5      | Signal time bin width                           |

---

## Output

All output is written to `<out_dir>/V<voltage>V__n<n_events>/`.

### ROOT file (`tgc_sim.root`)

The ROOT file contains one subdirectory per source distance (e.g. `dist_0p7mm/`) and a
`summary/` directory.

**Per-distance histograms:**

| Object                  | Type     | Description                                     |
|-------------------------|----------|-------------------------------------------------|
| `h_anode_charge`        | TH1D     | Total induced charge on all wires [fC]          |
| `h_cathode_charge`      | TH1D     | Total induced charge on bottom cathode [fC]     |
| `h_ratio_charge`        | TH1D     | Q_cathode / Q_anode per event                   |
| `h_n_clusters`          | TH1D     | Photoabsorption clusters per event (typically 1)|
| `h_n_primary_electrons` | TH1D     | Primary electrons from TrackHeed per event      |
| `h_avalanche_size`      | TH1D     | Total avalanche electrons per event             |
| `p_anode_signal`        | TProfile | Mean induced current on anode vs time [fC/ns]   |
| `p_cathode_signal`      | TProfile | Mean induced current on cathode vs time [fC/ns] |

**Summary graphs (in `summary/`):**

| Object           | Type         | Description                                   |
|------------------|--------------|-----------------------------------------------|
| `g_anode_charge` | TGraphErrors | ⟨Q_anode⟩ ± SEM vs source distance [fC]      |
| `g_cathode_charge` | TGraphErrors | ⟨Q_cathode⟩ ± SEM vs source distance [fC]  |
| `g_charge_ratio` | TGraphErrors | ⟨Q_cathode/Q_anode⟩ ± SEM vs source distance |

### CSV file (`summary.csv`)

One row per source distance with columns:

```
source_distance_mm, n_events, n_interacted, interaction_fraction,
mean_anode_charge_fC, rms_anode_charge_fC, sem_anode_charge_fC,
mean_cathode_charge_fC, rms_cathode_charge_fC, sem_cathode_charge_fC,
mean_charge_ratio, rms_charge_ratio, sem_charge_ratio,
mean_primary_electrons, mean_avalanche_size
```

### Config echo (`run_config.json`)

The resolved configuration used for the run, serialised to JSON for reproducibility.

### Summary PNG (`summary/tgc_summary.png`)

Three-panel figure: ⟨Q_anode⟩, ⟨Q_cathode⟩, and charge ratio vs source distance.

---

## Interpreting results

### Anode signal shape (`p_anode_signal`)

The average anode waveform has two components:

* **Fast electron component** (0–20 ns): electrons drift from their production point to
  the nearest wire.  The anode weighting potential peaks sharply near the wire surface,
  so the fast component carries the bulk of the induced charge.
* **Slow ion component** (up to ~300 ns): positive ions (CO2⁺) drift away from the
  wire toward both cathodes.  Their contribution to the anode signal has the same
  (positive) polarity as the electron component, producing a long positive tail.

The unprocessed simulation output (`p_anode_signal`) therefore shows a positive pulse that
rises sharply in the first ~20 ns and decays slowly over several hundred nanoseconds — it
is not bipolar.  A bipolar shape would appear after applying a differentiating RC filter
(as real front-end amplifiers do), but no such shaping is modelled here.

### Cathode signal shape (`p_cathode_signal`)

The cathode signal is dominated by the slow ion component: as CO2⁺ ions drift from the
wire plane toward the readout cathode, the cathode weighting potential rises monotonically
and the induced current accumulates over ~100–300 ns.  The resulting signal is a slow,
monotonically rising positive pulse.

### Charge ratio vs source distance

By the Ramo-theorem identity (the weighting potentials of all electrodes sum to unity),
the total induced charge across all electrodes is zero for any fully collected charge:

```
Q_anode + Q_cathode_bottom + Q_cathode_top = 0
```

In magnitude this means the two cathodes together collect as much charge as the wires.
Because the TGC has equal cathode-anode gaps on both sides of the wire plane (both
1.4 mm), and because the avalanche always occurs close to the wire, the resulting ion
cloud splits approximately equally between the top and bottom cathodes.  With only the
bottom cathode instrumented, the expected charge ratio is:

```
Q_cathode / Q_anode ≈ 0.5
```

Deviations from 0.5 can arise from slight asymmetry in the ion path (the source
illuminates from below, so the primary electrons arrive at the wire from the −y half of
the gap) and from event-to-event fluctuations in the avalanche.  The primary-electron
drift contribution to the ratio is negligible because the avalanche amplifies the charge
by a factor ~10⁴ before collection.

### Interaction fraction

Because the 5.9 keV photon mean free path (~3–5 cm in Ar:CO2 at 1 atm) is much longer
than the 1.4 mm gap, only ~3 % of events produce an interaction.  The remaining ~97 %
are skipped and do not contribute to histograms.  To increase statistics, increase
`n_events`; the histograms already contain only interacting events.

---

## Performance notes

| Config                   | Gas generation     | Events/distance | Typical wall time    |
|--------------------------|--------------------|-----------------|----------------------|
| Smoke (`smoke_tgc.json`) | ~2 min (if needed) | 10              | ~5 min               |
| Default (1000 events)    | ~10 min (once)     | 1000            | 1–4 h per distance   |
| Fast check               | cached             | 50              | ~15 min per distance |

The dominant cost is `AvalancheMicroscopic`: each of the ~220 primary electrons runs a
full microscopic avalanche.  To speed up:

* Reduce `n_events` (e.g. 100 for exploratory runs).
* Reduce `max_avalanche_size` (e.g. 10 000); this caps gain fluctuations but gives
  faster mean estimates.
* Use `DriftLineRKF` instead of `AvalancheMicroscopic` for drift-only simulations
  (no multiplication, much faster).
* Run multiple `--distance` jobs in parallel on separate cores.

Gas generation (Magboltz) is a one-time cost; the `.gas` file is reused on every
subsequent run.
