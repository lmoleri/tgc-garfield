# TGC Detector Simulation

A Garfield++ simulation of a Thin Gap Chamber (TGC) multi-wire proportional detector.
The simulation models primary ionisation from 5.9 keV Fe-55 X-rays, electron drift and
avalanche multiplication, and the induced charge on both the wire plane (anode) and one
cathode plane.

---

## Detector geometry

A TGC consists of an array of thin anode wires stretched between two grounded cathode
planes. The small cathode-anode gap (here 1.4 mm) produces a steep electric field
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

| Parameter             | Value          | Notes                              |
|-----------------------|----------------|------------------------------------|
| Wire count            | 10             |                                    |
| Wire diameter         | 50 μm          | radius = 25 μm                     |
| Wire pitch            | 1.8 mm         | centre-to-centre spacing           |
| Cathode-anode gap     | 1.4 mm         | distance from wire plane to cathode|
| Wire voltage          | +1900 V        | configurable via `wire_voltage_V`  |
| Cathode voltage       | 0 V            | both planes grounded               |
| Gas                   | Ar:CO2 70:30   | 1 atm, 20 °C                       |

---

## Physics

### 1. Gas transport coefficients — Magboltz

Electron drift velocity, diffusion, attachment, and Townsend coefficients are computed
by the Magboltz Monte Carlo code (via `MediumMagboltz`) over a logarithmic electric-field
grid from 100 V/cm to 300 kV/cm.  Results are cached to a `.gas` file and loaded on
subsequent runs.

**Penning transfer** is enabled by default.  In Ar:CO2 mixtures the metastable argon
levels (11.5–11.7 eV) can ionise CO2 (ionisation potential 13.78 eV — wait, that is
above Ar metastables, so the primary Penning channel is Ar* → Ar+ + CO2). The exact
Penning rate is taken from Garfield++'s built-in tabulation.

**Ion mobility** — after avalanche multiplication, positive ions must drift back to the
cathode to complete the Ramo-theorem induced-charge calculation.  In Ar:CO2 mixtures the
dominant drifting ion species is CO2+ rather than Ar+: the lower ionisation potential of
CO2 (13.78 eV vs. Ar 15.76 eV) means that Ar+ rapidly charge-transfers to CO2+ on
nanosecond timescales.  The simulation loads the `IonMobility_CO2+_CO2.txt` table from
the Garfield++ data directory.  Garfield++ accepts a single positive-ion mobility table
per gas object; CO2+ is the best single-species approximation for this mixture.

### 2. Primary ionisation — TrackHeed

The 5.9 keV photon is transported by `TrackHeed::TransportPhoton`.  Heed uses detailed
cross-section tables to simulate:

* **Photoelectric absorption** in the gas (dominant at 5.9 keV in Ar:CO2).  The mean
  free path at 1 atm is several centimetres, so only a small fraction of photons interact
  in the 1.4 mm gap; the simulation skips events where no interaction occurs and reports
  the interaction fraction in the summary.
* **Delta-electron cascade** — the photoelectron slows down and produces secondary
  ionisation along its track.  With W ≈ 26 eV/pair for Ar:CO2, a 5.9 keV photon yields
  roughly 220 primary electron-ion pairs.  All conduction electrons are returned with
  position and kinetic energy ready for `AvalancheMicroscopic`.
* **Fluorescence photons** (Ar K-alpha, ~2.96 keV) may appear in `cluster.photons`.
  They are currently not re-transported; adding a second `TransportPhoton` call for each
  fluorescence photon is a straightforward extension.

### 3. Electron avalanche — AvalancheMicroscopic

Each primary electron is transported individually by `AvalancheMicroscopic`, which steps
electrons through the gas using the Runge-Kutta-Fehlberg algorithm and samples elastic,
inelastic, ionising, and attachment collisions from the Magboltz cross-section tables.
All avalanche electrons are tracked until they are collected by a wire or absorbed.

An `avalanche_size_limit` cap prevents runaway events from consuming excessive CPU.

### 4. Signal induction — Shockley-Ramo theorem

`Sensor` computes induced current on each electrode at every step using the weighting
(Ramo) field:

> i(t) = q · v(t) · E_w(x(t))

where `E_w` is the weighting field of the electrode (the field that would exist if that
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

To skip regeneration on subsequent runs, the `.gas` file must be on the path given by
`gas.gas_file` in the config (default: `ar_70_co2_30.gas` in the current working
directory).  You can move or rename it and update the config key accordingly.

To regenerate with higher accuracy, increase `n_magboltz_collisions` (20–50 is
recommended for publication-quality results) and delete the existing `.gas` file.

---

## Build

From the `projects/tgc/` directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DVDT_INCLUDE_DIR=/Users/luca/miniforge3/include \
  -DVDT_LIBRARY=/Users/luca/miniforge3/lib/libvdt.dylib \
  -DCMAKE_PREFIX_PATH="/path/to/garfield++;/path/to/miniforge3"
cmake --build build -j4
```

Substitute the actual paths for your installation.  The concrete paths in this workspace are:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DVDT_INCLUDE_DIR=/Users/luca/miniforge3/include \
  -DVDT_LIBRARY=/Users/luca/miniforge3/lib/libvdt.dylib \
  -DCMAKE_PREFIX_PATH="$(pwd)/../../local/garfield;/Users/luca/miniforge3"
cmake --build build -j4
```

**Note:** `source ../../local/garfield/share/Garfield/setupGarfield.sh` sets the
`GARFIELD_INSTALL` and `HEED_DATABASE` environment variables needed at *runtime*, but
CMake find-package resolution requires the explicit `CMAKE_PREFIX_PATH`.  Also, the ROOT
conda package depends on `Vdt` which must be pointed to explicitly via `VDT_INCLUDE_DIR`
and `VDT_LIBRARY`.

Requirements: Garfield++ installed at `../../local/garfield`, ROOT 6 (via conda),
CMake ≥ 3.20, a C++20-capable compiler.

---

## Usage

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

The smoke CTest (50× faster) validates the build without generating a gas file from scratch
(assumes `ar_70_co2_30.gas` already exists):

```bash
ctest --test-dir build --verbose
```

---

## Configuration reference

All parameters live in a JSON file (default: `config/default_tgc.json`).

### `geometry`

| Key                | Type   | Unit | Default | Description                          |
|--------------------|--------|------|---------|--------------------------------------|
| `wire_pitch_cm`    | float  | cm   | 0.18    | Centre-to-centre wire spacing        |
| `wire_diameter_um` | float  | μm   | 50.0    | Wire outer diameter                  |
| `gap_cm`           | float  | cm   | 0.14    | Distance from wire plane to cathode  |
| `n_wires`          | int    | —    | 10      | Number of anode wires                |
| `wire_voltage_V`   | float  | V    | 1900.0  | Voltage applied to all wires         |

### `source`

| Key                    | Type         | Unit | Default              | Description                                              |
|------------------------|--------------|------|----------------------|----------------------------------------------------------|
| `energy_keV`           | float        | keV  | 5.9                  | Photon energy (Fe-55 K-alpha line)                       |
| `source_distances_mm`  | float array  | mm   | [0.2, 0.5, 0.9, 1.2] | Source y-distances from wire plane to scan               |
| `x_position_cm`        | float or null| cm   | null                 | Fixed photon x-position; null = uniform random over wires|

`source_distances_mm` is the scan axis.  Each value is the y-coordinate (in mm) at which
the photon is placed, measured from the wire plane (y = 0).  Values must be in the range
(0, `gap_cm` × 10).  The photon always travels in the −y direction.

### `gas`

| Key                      | Type   | Unit | Default             | Description                                       |
|--------------------------|--------|------|---------------------|---------------------------------------------------|
| `temperature_K`          | float  | K    | 293.15              | Gas temperature                                   |
| `pressure_Torr`          | float  | Torr | 760.0               | Gas pressure                                      |
| `gas_file`               | string | —    | "ar_70_co2_30.gas"  | Path to gas table (generated if absent)           |
| `enable_penning`         | bool   | —    | true                | Enable Penning transfer (recommended for Ar:CO2)  |
| `n_magboltz_collisions`  | int    | —    | 10                  | Magboltz collision cycles per field point         |

### `simulation`

| Key                  | Type   | Unit | Default  | Description                                              |
|----------------------|--------|------|----------|----------------------------------------------------------|
| `n_events`           | int    | —    | 1000     | Number of photon events per source distance              |
| `max_avalanche_size` | int    | —    | 500000   | Electron count cap per `AvalancheElectron` call          |
| `time_window_ns`     | float  | ns   | 300.0    | Signal collection window                                 |
| `time_step_ns`       | float  | ns   | 0.5      | Signal time bin width                                    |

---

## Output

All output is written to `<out_dir>/V<voltage>V__n<n_events>/`.

### ROOT file (`tgc_sim.root`)

The ROOT file contains one subdirectory per source distance (e.g. `dist_0p7mm/`) and a
`summary/` directory.

**Per-distance histograms:**

| Object               | Type      | Description                                          |
|----------------------|-----------|------------------------------------------------------|
| `h_anode_charge`     | TH1D      | Total induced charge on all wires [fC]               |
| `h_cathode_charge`   | TH1D      | Total induced charge on bottom cathode [fC]          |
| `h_ratio_charge`     | TH1D      | Q_cathode / Q_anode per event                        |
| `h_n_clusters`       | TH1D      | Number of photoabsorption clusters (always 1)        |
| `h_n_primary_electrons` | TH1D   | Primary electrons from TrackHeed per event           |
| `h_avalanche_size`   | TH1D      | Total avalanche electrons per event                  |
| `p_anode_signal`     | TProfile  | Mean induced current on anode vs time [fC/ns]        |
| `p_cathode_signal`   | TProfile  | Mean induced current on cathode vs time [fC/ns]      |

**Summary graphs (in `summary/`):**

| Object                 | Type         | Description                                     |
|------------------------|--------------|-------------------------------------------------|
| `g_anode_charge`       | TGraphErrors | ⟨Q_anode⟩ ± SEM vs source distance [fC]        |
| `g_cathode_charge`     | TGraphErrors | ⟨Q_cathode⟩ ± SEM vs source distance [fC]      |
| `g_charge_ratio`       | TGraphErrors | ⟨Q_cathode/Q_anode⟩ ± SEM vs source distance   |

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
  so the fast component carries most of the induced charge.
* **Slow ion component** (up to ~300 ns): positive ions (CO2+) drift toward the cathodes.
  Their contribution to the anode signal is opposite in sign to the electron component
  (electrons moving toward the anode and ions moving away both induce the same polarity
  on the anode), producing the characteristic bipolar shape after differentiation.

### Cathode signal shape (`p_cathode_signal`)

The cathode signal is dominated by the slow ion component (ions drifting toward the
cathode where the weighting potential rises monotonically).  It rises over ~100–300 ns
and is of the same polarity as the anode electron signal.

### Charge ratio vs source distance

A photon absorbed close to the wire plane (small `source_distance_mm`) produces primary
electrons that are immediately accelerated into the avalanche region.  The avalanche
occurs very close to the wire, so the electron induction is fast but the ion path to the
cathode is long: **Q_cathode / Q_anode is smaller** for small distances.

A photon absorbed close to the cathode (large `source_distance_mm`) produces primary
electrons that drift the full gap before avalanching near the wire.  The avalanche still
occurs close to the wire, but the ion path is essentially the same.  The primary-electron
drift path changes, which affects the timing of the signal but not the integrated induced
charge significantly.  The dominant effect on the charge ratio comes from the spread of
the avalanche position distribution relative to the cathode weighting field.

### Interaction fraction

Because the 5.9 keV photon mean free path (~3–5 cm in Ar:CO2 at 1 atm) is much longer
than the 1.4 mm gap, only ~3 % of events produce an interaction.  The remaining ~97 %
are skipped and do not contribute to histograms.  To increase statistics, increase
`n_events`; the histograms already contain only interacting events.

---

## Performance notes

| Config                  | Gas generation | Events/distance | Typical wall time |
|-------------------------|----------------|-----------------|-------------------|
| Smoke (`smoke_tgc.json`)| ~2 min (if needed) | 10           | ~5 min            |
| Default (1000 events)   | ~10 min (once) | 1000            | 1–4 h per point   |
| Fast check              | cached         | 50              | ~15 min per point |

The dominant cost is `AvalancheMicroscopic`: each of the ~220 primary electrons runs a
full microscopic avalanche.  To speed up:

* Reduce `n_events` (e.g. 100 for exploratory runs).
* Reduce `max_avalanche_size` (e.g. 10000); this caps gain fluctuations but gives
  faster mean estimates.
* Use `DriftLineRKF` instead of `AvalancheMicroscopic` for drift-only simulations
  (no multiplication, much faster).
* Run multiple `--distance` jobs in parallel on separate cores.

Gas generation (Magboltz) is a one-time cost; the `.gas` file is reused on every
subsequent run.
