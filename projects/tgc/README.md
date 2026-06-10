# TGC Detector Simulation

A Garfield++ simulation of a Thin Gap Chamber (TGC) multi-wire proportional detector.
The simulation deposits primary ionisation electrons at a configurable depth in the gas
gap, transports them through avalanche multiplication, and records the induced charge on
both the wire plane (anode) and one cathode plane.  Results can be explored via a PyQt5
desktop GUI or the command-line binary directly.

---

## Project structure

```
projects/tgc/
├── src/
│   └── tgc_sim.cc              ← simulation binary (C++20, ~1350 lines)
├── config/
│   ├── default_tgc.json        ← production config
│   ├── smoke_tgc_2.json        ← fast smoke test (ncoll=2, 10 events)
│   └── smoke_tgc_5.json        ← medium smoke test (ncoll=5)
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
  y = +1.4 mm  ───────────────────────────────  cathode_top (0 V, ground / non-readout)
               ||||||||    gas gap    |||||||||
  y =  0.0 mm  ─ ○ ─ ○ ─ ○ ─ ○ ─ ○ ─ ○ ─ ○ ─  anode wires (+1900 V)
               ||||||||    gas gap    |||||||||
  y = -1.4 mm  ───────────────────────────────  cathode (0 V, readout pad)

               ← 1.8 mm →
               wire pitch

  Source distance sign convention
  ─────────────────────────────────────────────
   0 mm → wire plane centre (y = 0)
  +d mm → readout pad side  (y = −d/10 cm)
  −d mm → cathode_top side  (y = +d/10 cm)
```

### Resistive readout option

When `readout.type = "resistive"`, the bottom cathode plane is replaced by a
layered structure.  The gas boundary condition is unchanged (the resistive layer
acts as a grounded conductor for DC fields), but the Ramo weighting potential
and the cathode signal shape are modified.

```
  y = −gap     ─── resistive layer (infinitely thin, ρ_s [Ω/sq]) ───
               ███████  insulator (Kapton/FR4, thickness d)  ███████
  y = −gap−d   ────────────── conductive readout pads ──────────────
```

The resistive layer is grounded at its four edges.  Deposited charge remains at
its landing point but the local surface potential decays with time constant
τ = ε₀ ε_r ρ_s L²/(π² d), where L = nWires × wirePitch / 2.

| Parameter         | Value        | Notes                               |
|-------------------|--------------|-------------------------------------|
| Wire count        | 10           |                                     |
| Wire diameter     | 50 μm        | radius = 25 μm                      |
| Wire pitch        | 1.8 mm       | centre-to-centre spacing            |
| Cathode-anode gap | 1.4 mm       | distance from wire plane to cathode |
| Wire voltage      | +1900 V      | configurable via `wire_voltage_V`   |
| Cathode voltage   | 0 V          | both planes grounded                |
| Gas               | Ar:CO2 70:30 | 750 Torr, 20 °C                     |

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

### 2. Primary ionisation — W-value model

For each simulated event, **N** primary electrons are deposited at the source position:

```
N = round( E_photon [eV] / W [eV/pair] )
```

For 5.9 keV Fe-55 and W = 26 eV/pair (Ar:CO2 70:30) this gives N ≈ 227.

Rather than transporting all N electrons individually (which would multiply the CPU cost
by ~227), the simulation runs **one representative avalanche** and scales the resulting
induced charge and avalanche size by N.  This is exact for the mean Q_cathode/Q_anode
ratio: all electrons start at the same position, so the Shockley-Ramo weighting is
identical for each, and scaling is equivalent to superposing N independent avalanches.
The event-to-event fluctuation in the charge ratio is dominated by the single-avalanche
(Polya) variance rather than by the Poisson variance in N.

The W-value is configurable via `gas.w_value_eV` (default 26 eV).

### 3. Electron avalanche — AvalancheMicroscopic

One representative electron is transported by `AvalancheMicroscopic`, which steps
electrons through the gas using the Runge-Kutta-Fehlberg algorithm and samples elastic,
inelastic, ionising, and attachment collisions from the Magboltz cross-section tables.
All avalanche electrons are tracked until they are collected by a wire or absorbed.

A `max_avalanche_size` cap prevents runaway events from consuming excessive CPU.

### 3b. Ion drift — DriftLineRKF

`AvalancheMicroscopic` counts the ions produced (`ni`) but does not transport them.
After each avalanche, the simulation iterates over every electron track endpoint
returned by `GetElectronEndpoints()`: the **start position** of track 0 is the primary
photoionisation ion; the start positions of tracks 1…n are the positions of avalanche
ions created in ionising collisions near the wire.  Each ion is transported by
`DriftLineRKF::DriftIon()` using the CO2⁺ mobility loaded from the Garfield++ data
directory, and its Ramo-theorem induced current is added to the sensor.

This is the dominant CPU cost for events with large avalanches (see Performance notes).

### 4. Signal induction — Shockley-Ramo theorem

`Sensor` computes the induced current on each electrode at every step using the weighting
(Ramo) field:

> i(t) = q · v(t) · **E**_w(**x**(t))

where **E**_w is the weighting field of the electrode (the field that would exist if that
electrode were at 1 V and all others at 0 V, with all space charges removed).

Two readout channels are defined:

* **`anode`** — the sense wires (all wires by default; configurable via
  `geometry.sense_wires`) share this label; their weighting fields are summed
  automatically by `Sensor`.  Non-sense wires are labelled `"field"` — they
  shape the electrostatic field but are not added as readout electrodes.
* **`cathode`** — the bottom cathode plane at y = −1.4 mm.

#### Conductive mode (default)

Standard Ramo induction: the cathode weighting potential is computed by
`ComponentAnalyticField` (1 V on the cathode plane, 0 V on wires and top).

#### Resistive mode

Two corrections apply when `readout.type = "resistive"`:

1. **Dielectric attenuation** — the conductive pad sits behind an insulating
   substrate of permittivity ε_r and thickness d.  The 1-D Poisson solution
   gives a reduced weighting potential in the gas:

   > W(y) = α (y + gap) / gap,   α = ε_r · gap / (d + ε_r · gap)

   For Kapton (ε_r = 3.5) with d = 100 μm and gap = 1.4 mm, α ≈ 0.98.

2. **Delayed signal** — the deposited surface charge remains at its landing
   point but the grounded edges pull the local resistive-layer potential toward
   0 V with time constant τ = ε₀ ε_r ρ_s L²/(π² d).  This causes the
   weighting potential to decay as W(y,t) = W(y) · exp(−t/τ), which contributes
   a time-distributed signal on two timescales:
   - *During drift in the gas*: the time-varying weighting potential modifies
     how much each drifting charge induces on the pad at each moment.
   - *After collection*: the fixed surface charge couples to the pad through a
     decaying potential (exponential tail with characteristic time τ).

   Both contributions are computed automatically by Garfield++'s
   `ComponentUser::SetDelayedWeightingPotential` framework and are included in
   `sensor.GetSignal("cathode", k)` once `sensor.EnableDelayedSignal()` is active.
   The τ printed to stdout at startup can be used to set `time_window_ns` long
   enough to capture the desired fraction of the delayed charge.

---

## Gas file

The gas file name is derived automatically from the gas configuration parameters —
there is no `gas_file` key in the config.  The naming scheme is:

```
ar70_co2_30_T{T}_P{P}_Ee{Ee}_Ef{Ef}k_n{n}_c{c}_{pen|nopen}.gas
```

For the default config this produces:
`ar70_co2_30_T293_P750_Ee2000_Ef400k_n50_c10_pen.gas`

The file is written to (and looked up from) the working directory of the binary,
which is `projects/tgc/` when run via the GUI or with the standard CMake invocation.

**On first run**: if the file does not exist, Magboltz generates it.  This takes
roughly **5–15 minutes** for the default settings (`n_magboltz_collisions = 10`).
The GUI shows `[will be generated]` next to the derived name before you click Run.

**On subsequent runs**: the file is loaded instantly.

**To regenerate with higher accuracy**: increase `n_magboltz_collisions` (20–50 is
recommended for publication-quality results).  Because the name encodes `c{n}`, a
new filename is derived automatically and the old cached file is left untouched.

**To regenerate from scratch**: delete the `.gas` file from `projects/tgc/`.

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
# Run with default config (2 source distances, 1 event each):
./build/tgc_sim

# Custom config:
./build/tgc_sim --config config/default_tgc.json --out results/

# Single source distance (quick check):
./build/tgc_sim --distance 0.7

# All options:
./build/tgc_sim --help
```

The smoke CTest validates the build without generating a gas file from scratch.
Two smoke configs are provided: `smoke_tgc_2.json` (ncoll=2, fastest) and
`smoke_tgc_5.json` (ncoll=5, slightly higher accuracy).

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
┌─────────────────────────────────────────────────────────────────────────────────────────┐
│  TGC Simulation          [▶ Run]  [■ Stop]  [Load Config]  [Save Config]                │
├─────────────────────────────────┬───────────────────────────────────────────────────────┤
│  ▼ Geometry                     │  [ Log | Summary | Plots | Waveforms |                │
│    Wire pitch [cm]   0.18       │    Charge | E-Field | 3D Tracks | Magboltz ]           │
│    Wire diameter [μm] 50        │                                                       │
│    Gap [cm]          0.14       │  ← live output / results shown here                   │
│    N wires           10         │                                                       │
│    Wire voltage [V]  1900       │                                                       │
│    Sense wires  [✓ All wires]   │                                                       │
│               [e.g. 4,5]        │                                                       │
│  ▼ Readout                      │                                                       │
│    Type        [Conductive ▼]   │                                                       │
│    (Insulator  [Kapton ▼])      │                                                       │
│    (Thickness  100 μm)          │                                                       │
│    (Resistivity 500 kΩ/sq)      │                                                       │
│    (Delayed signal [✓])         │                                                       │
│  ▼ Source                       │                                                       │
│    Energy [keV]      5.9        │                                                       │
│    Distances [mm]  -0.7,0.7     │                                                       │
│    X position  [ Random]        │                                                       │
│      fixed x [cm]  [0.0,0.9]   │  ← comma-separated; one run per value                 │
│  ▼ Gas                          │                                                       │
│    Gas 1  [ar  ▼]  70.0 %      │                                                       │
│    Gas 2  [co2 ▼]  30.0 %      │                                                       │
│    Ion species  [co2 ▼]         │                                                       │
│    Temperature [K]   293.15     │                                                       │
│    Pressure [Torr]   750        │                                                       │
│    Gas file  [ar_70…gas] […]   │                                                       │
│    Penning  [✓]  ncoll  10      │                                                       │
│    W-value [eV]  26.0           │                                                       │
│  ▼ Simulation                   │                                                       │
│    Events         1            │                                                       │
│    Max aval. size 500000        │                                                       │
│    Time window [ns] 40000       │                                                       │
│    Time step [ns]   0.5         │                                                       │
│    Ion transport  [✓]           │                                                       │
│    Store drift lines [✓]        │                                                       │
│  ▼ Output                       │                                                       │
│    Directory  [results/] […]    │                                                       │
│    Run name   [auto (date+V+n)] │                                                       │
└─────────────────────────────────┴───────────────────────────────────────────────────────┘
```

| Tab | Contents |
|---|---|
| **Log** | Live stdout stream from `tgc_sim`, auto-scrolling |
| **Summary** | Table from `summary.csv` — one row per source distance |
| **Plots** | 2 × 3 matplotlib figure: ⟨Q_anode⟩, ⟨Q_cathode⟩, ⟨Q_cathode_top⟩, charge ratio, and avalanche size vs source distance (with SEM error bars). Sixth cell empty |
| **Waveforms** | Mean anode and cathode current waveforms overlaid per (distance, x-position) combination. A distance selector and (when fixed x-positions were simulated) an x-position dropdown choose the folder to display. Read directly from the ROOT file via uproot |
| **Charge** | Cumulative charge integrals Q(t) — running integral of each waveform — for anode and cathode, per (distance, x-position) pair. An event slider selects individual events. ROOT TCanvas opens separately (PyROOT required) |
| **E-Field** | Interactive 2D electric field map in any of the XY, XZ, or YZ planes at a configurable depth; binning configurable from 50 to 10 000 bins per axis (PyROOT required) |
| **3D Tracks** | Per-event 3D detector view in a ROOT TCanvas showing detector geometry and drift lines with correct aspect ratios. Controls: preset view buttons (Gap XY / Top XZ / Side YZ / 3D reset), zoom ± (down to 0.5 % of full range), pan X/Y/Z. Distance and x-position selectors mirror the simulated folder structure. Wires rendered as semi-transparent 12-sided tube wireframes at actual diameter (clipped to the visible frame); cathode planes clipped to visible cube. Primary electron and ion drift lines colour-coded (blue / green / magenta / grey) and semi-transparent (PyROOT required) |
| **Magboltz** | Gas transport-property viewer: reads the `_props.csv` sidecar file exported automatically alongside the `.gas` file. Displays electron drift velocity, Townsend α, attachment η, longitudinal and transverse diffusion coefficients, effective gain (α − η), ion drift velocity, and ion mobility as a function of E-field — all in one ROOT TCanvas (8 panels). Export buttons save the plots to a ROOT file or copy the CSV. Requires PyROOT |

**▶ Run** starts the simulation in a background thread (the window stays fully
responsive).  **■ Stop** sends SIGTERM.  **Load Config** / **Save Config** read and
write `.json` files that are fully compatible with the CLI `--config` flag.

> **Note:** The E-Field, 3D Tracks, and Magboltz tabs open ROOT TCanvas windows and require
> PyROOT (ROOT importable from Python — available automatically when using the
> conda ROOT installation described in the Build section above).

---

## Configuration reference

All parameters live in a JSON file (default: `config/default_tgc.json`).

### `geometry`

| Key                | Type  | Unit | Default | Description                         |
|--------------------|-------|------|---------|-------------------------------------|
| `wire_pitch_cm`    | float | cm   | 0.18    | Centre-to-centre wire spacing. Sets the periodicity of the analytic field solution; changing this alters the electric-field map and weighting fields computed by `ComponentAnalyticField` |
| `wire_diameter_um` | float | μm   | 50.0    | Wire outer diameter (radius = value/2 μm). Sets the inner boundary of the avalanche region; thinner wires produce a higher peak field and larger gas gain for the same voltage |
| `gap_cm`           | float | cm   | 0.14    | Distance from the wire plane to **each** cathode (the geometry is symmetric: both gaps equal this value). Increasing the gap reduces the average drift field, lowering gain |
| `n_wires`          | int   | —    | 10      | Number of anode wires in the simulation cell. More wires increase the sensitive area but do not change single-wire physics |
| `wire_voltage_V`   | float | V    | 1900.0  | High voltage applied to all anode wires (cathodes grounded). Primary handle for tuning gas gain; a ~100 V change shifts gain by roughly one order of magnitude |
| `sense_wires`      | int array or null | — | `null` | 0-based indices (0 = leftmost wire) of the wires summed into the anode readout channel. `null` or absent → all wires read out (backward-compatible default). Non-listed wires remain at full HV and contribute to the electric field but are excluded from the Ramo weighting calculation. Validated: non-empty, every index in `[0, n_wires)` |

### `readout`

| Key                            | Type   | Unit  | Default         | Description |
|--------------------------------|--------|-------|-----------------|-------------|
| `type`                         | string | —     | `"conductive"`  | Cathode readout model. `"conductive"`: standard grounded plane (Ramo theorem only). `"resistive"`: adds insulating substrate and resistive layer — see Physics section |
| `insulator_material`           | string | —     | `"kapton"`      | Insulating substrate material. Sets the relative permittivity used for the dielectric correction and τ calculation. `"kapton"` → ε_r = 3.5; `"fr4"` → ε_r = 4.6. Ignored when `type = "conductive"` |
| `insulator_thickness_um`       | float  | μm    | 100.0           | Thickness of the insulating substrate between the resistive layer and the conductive pads. Affects both the dielectric correction factor α and the time constant τ. Ignored when `type = "conductive"` |
| `surface_resistivity_ohm_sq`   | float  | Ω/sq  | 500000.0        | Sheet resistance of the resistive layer. Enters only the time constant τ (does not affect the static field or α). Ignored when `type = "conductive"` |
| `enable_delayed_signal`        | bool   | —     | `true`          | When `false`, skips `SetDelayedWeightingPotential` — the static α-corrected weighting potential is still applied, but `sensor.EnableDelayedSignal()` is not called. Removes the ~200× per-step overhead for exploratory runs. Ignored when `type = "conductive"` |

### `source`

| Key                   | Type          | Unit | Default              | Description                                               |
|-----------------------|---------------|------|----------------------|-----------------------------------------------------------|
| `energy_keV`          | float         | keV  | 5.9                  | Photon energy of the simulated X-ray source (Fe-55 K-alpha line). Determines the number of primary electrons via `N = round(E_photon[eV] / w_value_eV)` (≈227 for Fe-55 at 5.9 keV and W = 26 eV) |
| `source_distances_mm` | float array   | mm   | [-0.7, 0.7] | List of signed y-positions (mm) at which primary electrons are placed, measured from the wire plane. Positive → readout cathode side (y < 0); negative → cathode_top side (y > 0). Each distance is a separate simulation run. Values are clamped to (−`gap_cm`×10, +`gap_cm`×10) |
| `x_positions_cm`      | float array or null | cm | [0.0, 0.9] | Fixed lateral (x) positions [cm] for the photon interaction point. `null` → uniform random over the wire array each event; the shipped default specifies two fixed x-positions (wire centre and midpoint). One or more values → one simulation run per distance × x-position pair; ROOT directory named `dist_Nmm_xMmm/`. Backward-compatible with the old scalar `x_position_cm` key (wrapped into a one-element list) |

### `gas`

| Key                     | Type   | Unit  | Default                   | Description                                      |
|-------------------------|--------|-------|---------------------------|--------------------------------------------------|
| `temperature_K`         | float  | K     | 293.15                    | Gas temperature passed to Magboltz. Affects gas number density (n ∝ 1/T at fixed pressure), which shifts drift velocity and Townsend coefficients. Must match physical detector conditions; 293.15 K = 20 °C |
| `pressure_Torr`         | float  | Torr  | 750.0                     | Gas pressure passed to Magboltz. Together with temperature, sets gas density. 750 Torr ≈ 0.987 atm. Reducing pressure increases mean free path and electron energy, raising gain |
| `enable_penning`        | bool   | —     | true                      | Activates Penning transfer via `MediumMagboltz::EnablePenningTransfer()`. In Ar:CO2 70:30 this raises the effective Townsend coefficient by ~20–40 %, bringing simulated gain closer to measured values. Should be left on for this mixture |
| `n_magboltz_collisions` | int    | —     | 10                        | Monte Carlo collision cycles Magboltz runs per E-field grid point. Higher values reduce statistical uncertainty in transport coefficients at the cost of longer gas-file generation. 2–5: smoke test; 10: default; 20–50: publication quality |
| `max_electron_energy_eV`| float  | eV    | 2000.0                    | Upper energy bound for the Magboltz cross-section look-up table. Must exceed the maximum kinetic energy electrons reach near the wire (typically 500–1000 eV here). Too low a ceiling causes Magboltz to extrapolate, producing unphysical transport coefficients |
| `n_field_points`        | int    | —     | 50                        | Number of logarithmically spaced E-field values at which Magboltz computes transport coefficients (from ~100 V/cm to `e_field_max_vcm`). More points give smoother interpolation; fewer points speed up gas-file generation |
| `e_field_max_vcm`       | float  | V/cm  | 400000.0                  | Upper E-field limit for the Magboltz transport table. Must exceed the peak near-wire field E_peak = V_wire / (r × (ln(pitch/(2π r)) + π gap/pitch)). For the default geometry E_peak ≈ 156 kV/cm (2.6× margin at the default 400 kV/cm). The binary prints E_peak at startup and warns when `e_field_max_vcm < E_peak`; recommends ≥ 1.5× margin. The GUI "Auto" button fills in 2× E_peak rounded to the next 50 kV/cm |
| `w_value_eV`            | float  | eV    | 26.0                      | Mean energy to create one electron–ion pair in the gas mixture (W-value). Determines primary electron count: `N = round(energy_keV × 1000 / w_value_eV)`. The measured value for Ar:CO2 70:30 is ~26 eV |

### `simulation`

| Key                  | Type  | Unit | Default  | Description                                     |
|----------------------|-------|------|----------|-------------------------------------------------|
| `n_events`           | int   | —    | 1        | Number of avalanche simulations per source distance. The shipped default is 1 (single-event inspection); set to 1000 for ~3 % SEM on gain. More events reduce statistical uncertainty on mean charge and charge ratio (SEM ∝ 1/√N) |
| `max_avalanche_size` | int   | —    | 500000   | Maximum number of electrons tracked per `AvalancheMicroscopic` call. Truncates runaway avalanches to prevent excessive CPU use. Reduce to ~10 000 for fast exploratory runs (biases the high-gain tail). Note: each avalanche electron corresponds to one `DriftLineRKF::DriftIon` call, so smaller values also reduce ion-drift CPU cost |
| `time_window_ns`     | float | ns   | 40000.0  | Duration of the induced-current waveform recorded on each electrode. The shipped default of 40 μs captures the full ion drift (~5–8 μs). Use ~300 ns for electron-signal-only studies: 300 ns captures the full electron component (≲20 ns) and the first ~34 % of the slow ion tail |
| `time_step_ns`       | float | ns   | 0.5      | Width of each time bin in the `TProfile` waveforms (`p_anode_signal`, `p_cathode_signal`). Finer bins give better time resolution but larger ROOT histograms. 0.5 ns is sufficient to resolve the fast electron peak (~5–10 ns FWHM) |
| `enable_ion_drift`   | bool  | —    | true     | Drift positive ions after each electron avalanche using `DriftLineRKF`. When enabled, every ion created during the avalanche is transported to a cathode and its Ramo-theorem induced current is added to the waveform. Disabling skips ion signal computation entirely, greatly reducing CPU time for large avalanches at the cost of losing the cathode signal and ion tail |
| `store_drift_lines`  | bool  | —    | true     | When `true`, `AvalancheMicroscopic` records every intermediate collision step in the primary electron drift line (not just start and end), producing denser 3D path data for the GUI 3D Tracks viewer at the cost of larger ROOT files |

---

## Output

All output is written to a date-stamped subdirectory inside `<out_dir>`:

| Mode | Folder name |
|------|-------------|
| Auto (no run name) | `yymmdd_hh-MM__V<voltage>V__n<n_events>/` |
| Custom run name | `yymmdd_hh-MM__<label>/` |

The date prefix (`yymmdd_hh-MM`) is generated at launch time.  From the GUI, set
"Run name" in the Output group to use a custom label; leave it blank for the auto
format.  When invoking the binary directly, pass `--run-name <label>` to override;
omitting it generates the auto format via `BuildRunFolderName` in `tgc_sim.cc`.

### Gas properties sidecar (`<gasfile>_props.csv`)

When the gas is set up, `tgc_sim` automatically writes a `<gasfile>_props.csv` sidecar
file next to the `.gas` file.  It contains 8 Magboltz transport coefficients (electron
drift velocity, Townsend α, attachment η, longitudinal and transverse diffusion, effective
gain α − η, ion drift velocity and ion mobility) sampled at each E-field grid point.
This file is read by the **Magboltz** tab in the GUI.

### ROOT file (`tgc_sim.root`)

The ROOT file contains one subdirectory per (distance, x-position) combination and a
`summary/` directory.  Directory naming:

| Condition | Example name |
|-----------|-------------|
| `x_positions_cm: null` (random x) | `dist_0p7mm/` |
| `x_positions_cm: [0.0]` | `dist_0p7mm_x0mm/` |
| `x_positions_cm: [0.0, 0.18]` | `dist_0p7mm_x0mm/`, `dist_0p7mm_x1p8mm/` |

The x-position suffix uses millimetres (× 10 relative to the cm config value) with
decimal points replaced by `p` (e.g. 0.18 cm → 1.8 mm → `x1p8mm`).

**Per-distance histograms:**

| Object                  | Type     | Description                                     |
|-------------------------|----------|-------------------------------------------------|
| `h_anode_charge`        | TH1D     | Induced charge on the sense wires integrated over the time window [fC] |
| `h_cathode_charge`      | TH1D     | Induced charge on bottom (readout) cathode integrated over the time window [fC] |
| `h_cathode_top_charge`  | TH1D     | Induced charge on top (non-readout) cathode integrated over the time window [fC] |
| `h_ratio_charge`        | TH1D     | Q_cathode / Q_anode per event                   |
| `h_n_primary_electrons` | TH1D     | Primary electrons per event (= round(E/W), constant) |
| `h_avalanche_size`      | TH1D     | Total avalanche electrons per event             |
| `p_anode_signal`        | TProfile | Mean induced current on anode vs time [fC/ns]   |
| `p_cathode_signal`      | TProfile | Mean induced current on cathode vs time [fC/ns] |
| `p_cathode_top_signal`  | TProfile | Mean induced current on cathode_top vs time [fC/ns] |

**Summary graphs (in `summary/`):**

| Object           | Type         | Description                                   |
|------------------|--------------|-----------------------------------------------|
| `g_anode_charge`       | TGraphErrors | ⟨Q_anode⟩ ± SEM vs source distance [fC]        |
| `g_cathode_charge`     | TGraphErrors | ⟨Q_cathode⟩ ± SEM vs source distance [fC]      |
| `g_cathode_top_charge` | TGraphErrors | ⟨Q_cathode_top⟩ ± SEM vs source distance [fC]  |
| `g_charge_ratio`       | TGraphErrors | ⟨Q_cathode/Q_anode⟩ ± SEM vs source distance   |

**Per-event tree (`t_signals`, one entry per simulated event):**

| Branch | Type | Description |
|--------|------|-------------|
| `event` | int | Event index (0-based) |
| `anode_charge_fC` | float | Integrated anode charge for this event [fC] |
| `cathode_charge_fC` | float | Integrated cathode charge for this event [fC] |
| `anode` | vector\<float\> | Per-bin anode current waveform [fC/ns], length = `time_window_ns / time_step_ns` |
| `cathode` | vector\<float\> | Per-bin cathode current waveform [fC/ns] |
| `primary_x/y/z` | vector\<float\> | 3D points along the primary electron drift path [cm]. 2 points (start + end) when `store_drift_lines = false`; full collision-step trajectory when `true` |
| `cloud_x/y/z` | vector\<float\> | Start positions of secondary (avalanche) electron tracks [cm], up to 500 points |
| `ion_x/y/z` | vector\<float\> | Flattened ion drift paths [cm], up to 100 ions |
| `ion_npts` | vector\<int\> | Number of drift-line points stored per ion path |

This tree is read by the GUI Waveforms, Charges, and 3D Tracks tabs via `uproot`.

---

### CSV file (`summary.csv`)

One row per (source distance, x-position) combination.  Columns:

| Column | Unit | Description |
|--------|------|-------------|
| `source_distance_mm` | mm | Signed y-distance of the photon interaction from the wire plane. Positive → readout cathode side (y < 0); negative → cathode_top side (y > 0). See sign convention in *Detector geometry* |
| `x_position_cm` | cm | Fixed lateral x-position of the interaction. Empty (blank field) when `x_positions_cm: null` — i.e. the x-position was drawn uniformly at random each event |
| `n_events` | — | Number of avalanche simulations run at this (distance, x-position) combination |
| `n_interacted` | — | Events that produced at least one primary electron. Currently always equals `n_events` — every event interacts by construction |
| `interaction_fraction` | — | `n_interacted / n_events`. Always 1.0 |
| `mean_anode_charge_fC` | fC | Mean induced charge on the sense wires, integrated over the time window. With the default 40 μs window the full ion tail is captured; with a shorter window (e.g. 300 ns) only ~34 % of the ion contribution is included |
| `rms_anode_charge_fC` | fC | RMS (σ) of the per-event anode charge distribution. Reflects Polya/exponential avalanche fluctuations |
| `sem_anode_charge_fC` | fC | Standard error of the mean (σ/√N): statistical uncertainty on `mean_anode_charge_fC` |
| `mean_cathode_charge_fC` | fC | Mean induced charge on the **readout (bottom) cathode** plane at y = −gap. Dominated by the slow CO2⁺ ion drift component; expected to be ~50 % of `mean_anode_charge_fC` by Ramo-theorem charge conservation |
| `rms_cathode_charge_fC` | fC | RMS of the per-event bottom-cathode charge distribution |
| `sem_cathode_charge_fC` | fC | SEM for `mean_cathode_charge_fC` |
| `mean_cathode_top_charge_fC` | fC | Mean induced charge on the **non-readout (top) cathode** plane at y = +gap. Together with the bottom cathode it completes the Ramo identity: Q_anode + Q_cathode + Q_cathode_top ≈ 0 |
| `rms_cathode_top_charge_fC` | fC | RMS of the per-event top-cathode charge distribution |
| `sem_cathode_top_charge_fC` | fC | SEM for `mean_cathode_top_charge_fC` |
| `mean_charge_ratio` | — | Mean of Q_cathode / Q_anode computed per event. Expected ≈ 0.5 for symmetric gaps; see *Charge ratio vs source distance* |
| `rms_charge_ratio` | — | RMS of the per-event charge-ratio distribution |
| `sem_charge_ratio` | — | SEM for `mean_charge_ratio` |
| `mean_primary_electrons` | — | Average primary electron count per event: `round(energy_keV × 1000 / w_value_eV)`. Constant for a fixed source energy and W-value (≈ 227 for Fe-55 at W = 26 eV) |
| `mean_avalanche_size` | — | Mean total electrons produced in one representative single-electron avalanche (≈ gas gain). Multiply by `mean_primary_electrons` to estimate the total collected electrons per Fe-55 event |

### Config echo (`run_config.json`)

The resolved configuration used for the run, serialised to JSON for reproducibility.

### Summary PNG (`summary/tgc_summary.png`)

Three-panel figure: ⟨Q_anode⟩, ⟨Q_cathode⟩, and charge ratio vs source distance.

### GUI plots snapshot (`tgc_plots.root`)

Written by the GUI at the end of each run (requires PyROOT).  Contains the
currently-rendered ROOT TCanvas objects, one per tab, keyed as follows:

| Key | Source tab | Contents |
|-----|------------|----------|
| `waveforms` | Waveforms  | Mean anode/cathode signal waveforms |
| `charge`    | Charges    | Cumulative charge integrals Q(t) |
| `tracks_3d` | 3D Tracks  | 3D drift-line and geometry view |
| `magboltz`  | Magboltz   | Gas transport-coefficient panels |
| `efield`    | E-Field    | 2D electric-field maps |

Only keys for canvases that were open and alive at run completion are written.
The file is omitted (and a warning logged) if PyROOT is unavailable.

---

## Interpreting results

### Anode signal shape (`p_anode_signal`)

The average anode waveform has two components:

* **Fast electron component** (0–20 ns): electrons drift from their production point to
  the nearest wire.  The anode weighting potential peaks sharply near the wire surface,
  so the fast component carries the bulk of the induced charge.
* **Slow ion component** (~5–8 μs physically): positive ions (CO2⁺) drift away from
  the wire toward both cathodes.  Their contribution to the anode signal has the same
  (positive) polarity as the electron component, producing a long positive tail.

The full ion tail extends over ~5–8 μs (estimated from the CO2⁺ reduced mobility
K₀ ≈ 1.7 cm²/(V·s) and the average gap field of ~13 600 V/cm).  With the default
40 μs window the ion tail is fully captured.  For a shorter window (e.g. 300 ns) only
the first ~34 % is collected: in 300 ns the ions travel roughly 75 μm out of the
1.4 mm gap, but the Ramo weighting potential changes most rapidly near the wire, so the
first fraction is disproportionately large.  `p_anode_signal` shows a fast peak from
electrons followed by a slowly decaying positive tail — it is not bipolar.  A bipolar
shape would appear after a differentiating RC filter, which is not modelled here.

### Cathode signal shape (`p_cathode_signal`)

The cathode signal is dominated by the slow ion component: as CO2⁺ ions drift from the
wire plane toward the readout cathode, the cathode weighting potential rises monotonically.
The full induction extends over ~5–8 μs as the ions travel the 1.4 mm gap.  With the
default 40 μs window the cathode signal reaches its plateau.  With a shorter 300 ns
window it is still rising, having reached roughly one-third of its final value.

In **resistive mode**, the cathode waveform includes both the prompt component
(from charge drifting in the gas, attenuated by factor α) and the delayed
component (exponential tail from the decaying surface potential).  The delayed
tail has characteristic time τ printed to stdout at startup; set
`time_window_ns` ≥ 5τ to capture most of the delayed charge.

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

Every event deposits primary electrons at the configured distance, so the interaction
fraction is always 100 %.  The `n_interacted` and `interaction_fraction` fields in
`summary.csv` equal `n_events` and 1.0 respectively.

---

## Performance notes

| Config                        | Gas generation     | Events/distance | Typical wall time   |
|-------------------------------|--------------------|-----------------|---------------------|
| Smoke (`smoke_tgc_2.json`)    | ~2 min (if needed) | 10              | ~1–2 min            |
| Default (`default_tgc.json`)  | ~10 min (once)     | 1 (shipped)     | seconds per distance|
| Production (1000 events)      | reuse cached file  | 1000            | ~20–60 min/distance |

The dominant cost is one `AvalancheMicroscopic` simulation per event (the result is then
scaled by N_primary ≈ 227).  To speed up:

* Reduce `n_events` (e.g. 50 for exploratory runs).
* Reduce `max_avalanche_size` (e.g. 10 000); this caps gain fluctuations but gives
  faster mean estimates.
* Run multiple `--distance` jobs in parallel on separate cores.

Gas generation (Magboltz) is a one-time cost; the `.gas` file is reused on every
subsequent run.
