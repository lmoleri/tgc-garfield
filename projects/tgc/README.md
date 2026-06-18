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
│   └── tgc_sim.cc              ← simulation binary (C++20, ~1650 lines)
├── config/
│   ├── default_tgc.json            ← production config
│   ├── default_tgc_fast-current.json ← fast preset (short window, narrow E-grid)
│   ├── smoke_tgc_2.json            ← fast smoke test (ncoll=2, 10 events)
│   └── smoke_tgc_5.json            ← medium smoke test (ncoll=5)
├── gui/
│   └── app.py                  ← PyQt5 desktop GUI
├── tools/
│   └── plot_ion_tail.py        ← fits i₀/(1+t/t₀) to the anode ion tail
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
  y = −gap     ── resistive layer (20 cm, ρ_s [Ω/sq], ∥-wire edges grounded) ──
               ███████  insulator (Kapton/FR4, thickness d)  ███████
  y = −gap−d   ──────────── conductive readout pad (current size) ───────────
               ▒▒▒▒▒▒▒  insulator (1 mm FR4)  ▒▒▒▒▒▒▒   ← optional, ground_plane_enabled
  y = −gap−d−1mm ──────────────── ground plane ──────────────────────────────
```

The resistive layer is a square sheet (default 20 × 20 cm, `resistive_layer_size_cm`)
grounded along its two edges **parallel to the wires**; the conductive readout pad behind it
keeps its usual (smaller) size.  Because only those two edges are grounded, deposited charge
relaxes *across* the wire array toward them, and the local surface potential decays with time
constant τ = ε₀ ε_r ρ_s L²/(π² d), where L = `resistive_layer_size_cm` / 2 (half the
across-wire span between the grounded edges).

Optionally (`ground_plane_enabled`) a **grounded plane** sits a further 1 mm (FR4) below the
readout pad.  It adds a pad-to-ground capacitance that reduces the pad signal (see Physics
§4 below); it does not change the DC drift field or τ, and has no effect in conductive
readout (a solid grounded cathode fully shields it).

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
grid from `e_field_min_vcm` to `e_field_max_vcm` (default 100 V/cm to 400 kV/cm).
Results are cached to a `.gas` file and reloaded on
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
timescales.  The simulation loads the `IonMobility_{X}+_{X}.txt` table from the
Garfield++ data directory, where `X` is `gas.ion_species` (default `co2`).
Garfield++ stores one positive-ion mobility table per gas object; CO2⁺ is the best
single-species approximation for the default Ar:CO2 mixture.

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
`DriftLineRKF::DriftIon()` using the configured ion mobility (`gas.ion_species`,
default CO2⁺) loaded from the Garfield++ data directory, and its Ramo-theorem induced
current is added to the sensor.

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
   substrate of permittivity ε_r and thickness d.  The spatial shape of the pad
   weighting potential is the **wire-screened analytic cathode weighting
   potential** of `ComponentAnalyticField` (same as conductive mode: 1 V on the
   cathode plane, 0 V on wires and top, so the grounded wires suppress it near
   the wire plane and in the top half-gap), scaled by the dielectric
   transparency factor:

   > W(x,y,z) = α · W_cathode(x,y,z),   α = ε_r · gap / (d + ε_r · gap)

   For Kapton (ε_r = 3.5) with d = 100 μm and gap = 1.4 mm, α ≈ 0.98.
   (An earlier version used a 1-D linear W(y) that ignored wire screening; it
   overstated the prompt electron spike and the top-going-ion contribution on
   the pad.)

   Equivalently α is a capacitive divider, α = C_ins / (C_ins + C_gap), with
   C_ins = ε₀ε_r/d (pad ↔ resistive layer) and C_gap = ε₀/gap (pad ↔ gas
   return).  When `ground_plane_enabled`, a grounded plane d₂ below the pad adds
   a third arm C_gnd = ε₀ε_r2/d₂ (pad ↔ ground) in parallel with the gas return:

   > α = C_ins / (C_ins + C_gap + C_gnd)

   so the backplane reduces the pad signal.  For the defaults (Kapton upper
   insulator, 1 mm FR4 below) α drops 0.98 → ≈ 0.87 (pad signal −11 %); a 1 mm
   *air* gap gives ≈ 0.95.  τ and the DC field are unchanged (the 1 mm backplane
   is too far to alter the sheet's relaxation capacitance, and the grounded
   cathode shields it from the gas), and the option has no effect in conductive
   readout.

   **Accuracy and validity of the α model.**  The α factor is the analytic
   *dielectric-transparency* approximation: it assumes the field is uniform across
   the gas/insulator interface, so the gas weighting potential keeps its exact
   wire-screened shape and is merely rescaled by the scalar α.  Only that scalar is
   approximate — and it is a good approximation here: the wire-induced ripple has
   decayed to ~exp(−2π·gap/pitch) ≈ 0.7 % by the time it reaches the interface
   (a full gap below the wires), so α·W_cathode is accurate to **~1 % for the
   default geometry**.  The error grows with d/gap and with ε_r.  Crucially, the
   weighting field is only ever evaluated at charge positions **inside the gas** —
   the gas medium is bounded by the analytic cathode plane at y = −gap, so every
   electron and ion is absorbed at the resistive layer and none enter the insulator
   (`ComponentAnalyticField` returns no medium behind the plane; the microscopic
   avalanche and the RKF ion drift terminate there).  The insulator's internal
   weighting potential therefore never enters the induced-signal calculation; the
   GUI **Weighting Field** tab fills it in (the α→1 ramp up to the pad) for
   visualisation only.  A fully rigorous layered solve — for thick-insulator or
   high-ε_r regimes where the ~1 % no longer holds — would require a 2-D boundary
   element model of the wires + dielectric + pad (Garfield `ComponentNeBem2d`) or an
   external FEM solution imported via `ComponentElmer` / `ComponentComsol`.

2. **Delayed signal** — the deposited surface charge remains at its landing
   point but the grounded edges pull the local resistive-layer potential toward
   0 V with time constant τ = ε₀ ε_r ρ_s L²/(π² d).  This causes the
   weighting potential to decay as W(x,y,z,t) = W(x,y,z) · exp(−t/τ), which
   contributes a time-distributed signal on two timescales:
   - *During drift in the gas*: the time-varying weighting potential modifies
     how much each drifting charge induces on the pad at each moment.
   - *After collection*: the fixed surface charge couples to the pad through a
     decaying potential (exponential tail with characteristic time τ).

   Because the delayed weighting potential is separable, W(x,y,z,t) =
   W(x,y,z)·exp(−t/τ), both contributions reduce *exactly* to a causal
   exponential filter applied to the binned prompt pad current:

   > i_pad(t) = i_prompt(t) − (1/τ) ∫₀ᵗ i_prompt(t′) e^{−(t−t′)/τ} dt′

   The simulation records the cheap prompt signal and applies this filter once
   per event (`ApplyResistiveRelaxation`).  This is mathematically identical to
   Garfield++'s per-drift-step `SetDelayedWeightingPotential` machinery for this
   model, but has no per-step cost (the old approach evaluated the weighting
   potential at 200 delayed times per drift step and spread each into the time
   bins — resistive runs now cost the same as conductive ones) and is exact at
   the time-bin resolution instead of a 200-point sampling.
   The τ printed to stdout at startup can be used to set `time_window_ns` long
   enough to capture the desired fraction of the delayed charge (≥ 5τ for the
   full relaxation; pair a long window with a coarser `time_step_ns` to keep the
   bin count reasonable).

#### Why the electron spike appears on the simulated cathode

The cathode waveform shows a sharp electron spike in **both** readout modes, and
this is expected.  In conductive mode it is direct Ramo induction: the avalanche
electrons deliver their (wire-screened, hence small) share of pad charge within
~1 ns, so even a few fC arriving that fast produces a tall current spike
(~4 fC/ns peak at default settings) — prominent in *current*, minor in *charge*.

In resistive mode the spike survives because a resistive sheet is **transparent
to fast signals**.  The sheet can only screen a lateral field disturbance of
scale λ on a time scale τ_screen(λ) ≈ ρ_s ε₀ ε_r λ²/d.  For the defaults
(500 kΩ/sq, 100 μm Kapton) and the ~1.4 mm induction footprint of an avalanche,
τ_screen ≈ 100–300 ns — the sheet cannot react within the ~1 ns spike, so the
induction passes through to the pad unattenuated.  (The same physics lets ATLAS
TGC pickup strips behind ~1 MΩ/sq graphite cathodes deliver 25 ns-class trigger
signals.)  Consistently, the exp(−t/τ) post-filter with τ = 157 μs removes only
~10⁻⁵ of a 1 ns spike.  Note the single-τ filter models global drainage to the
grounded edges, not this scale-dependent screening; that approximation is safe
for the spike whenever τ_screen(footprint) ≫ the spike duration, i.e. unless
ρ_s·ε_r/d is ~50× smaller than the defaults.

If a measured pad waveform shows no spike, the cause is usually the measurement
chain rather than the chamber: (a) the pad–sheet coupling is ~31 pF/cm²
(ε₀ε_r/d for 100 μm Kapton), so a large pad into 50 Ω gives a 15–50 ns input RC
that attenuates a 1 ns spike ×15–50; (b) a digitizer sampling at 10 ns averages
it down another ×10; (c) if the real coating's ρ_s is much lower than configured,
τ_screen drops into the few-ns range and the sheet itself starts eating the
spike; (d) the spike carries only a small fraction of the pad charge (wire
screening), so after any of the above it sits below the noise floor.

### 5. Front-end electronics — fast current amplifier

An **opt-in** model of the **CIVIDEC C2-TCT broadband current amplifier** can be
applied to the two physical readout channels — the **anode** (wires) and the
**cathode** (readout pad) — to produce the amplifier output voltage that a scope
would record, for direct comparison with measured waveforms.  It is **off by
default** (`amplifier.enable = false`), so a default run is byte-for-byte
unchanged; the raw induced-current waveforms are always kept, and the amplifier
output is written to **new** branches/profiles in mV.

Datasheet parameters: current amplifier, gain **40 dB** (×100), analog bandwidth
**10 kHz – 2 GHz**, input impedance **50 Ω**, AC-coupled input (**1 nF**),
non-inverting bipolar, ±1 V linear output.  Per the hardware setup the **wire**
input carries an **additional 470 pF capacitor in series**.

The amplifier is linear and time-invariant, so its effect is a filter cascade on
the per-bin induced current (note **i [fC/ns] ≡ i [µA]**, since fC/ns = 10⁻⁶ A):

> V_out(t) [mV] = G · R_in · ( LP_{2 GHz} ∘ HP_{coupling} ∘ HP_{10 kHz} )[ i(t) ] · 10⁻³

- **HP_coupling** — AC coupling from the series input capacitor and the 50 Ω input,
  a one-pole high-pass with τ = R_in · C_series.  This is the element that
  distinguishes the channels:
  - cathode (pad): C = 1 nF → τ ≈ **50 ns**
  - anode (wire): C = 470 pF ⊕ 1 nF = 320 pF → τ ≈ **16 ns** (faster baseline
    restoration / more differentiation than the pad)
- **HP_10 kHz** — the amplifier's intrinsic lower band edge, τ = 1/(2π·10 kHz) ≈
  15.9 µs (a slow droop; negligible in short windows, visible in long ones).
- **LP_2 GHz** — the upper band edge, τ = 1/(2π·2 GHz) ≈ 0.08 ns (smooths sub-ns
  features; only matters at very fine `time_step_ns`).
- **Scaling** — G = 10^(40/20) = 100 and R_in = 50 Ω give **V_out [mV] = 5 · i [fC/ns]**.

A one-pole high-pass is exactly the recursion already used for the resistive
relaxation (`i_out = i − (1/τ)·lowpass_τ(i)`); the low-pass is the standard
`y[k] = b·y[k−1] + (1−b)·x[k]`, `b = e^{−Δt/τ}`.  The chain is applied to the
post-relaxation pad current and the raw wire current (the electronics sit after
the detector).  No noise term is included (the model is deterministic; the ~4 µA
spike dwarfs the 0.4 µA rms input noise anyway).

Note that this **fast** current amplifier **keeps the electron spike**: a high-pass
passes fast signals and a 2 GHz roll-off does not smear a ~1 ns feature, so the
amplifier output still shows the spike on both channels.  The absence of the spike
in measured data therefore comes from elsewhere (input loading, grounding, the
actual ρ_s) — see *Why the electron spike appears on the simulated cathode* above.

---

## Gas file

The gas file name is derived automatically from the gas configuration parameters —
there is no `gas_file` key in the config.  The naming scheme is:

```
{gas1}{f1}_{gas2}_{f2}_T{T}_P{P}_Ee{Ee}_Ef{Efmin}v-{Efmax}k_n{n}_c{c}_{pen|nopen}.gas
```

For the default config this produces:
`ar70_co2_30_T293_P750_Ee2000_Ef100v-400k_n50_c10_pen.gas`

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
│  ▼ Amplifier                    │                                                       │
│    Enable        [ ]            │                                                       │
│    (Gain [dB]    40 …)          │  ← rows enabled when Enable is checked                │
│  ▼ Source                       │                                                       │
│    Energy [keV]      5.9        │                                                       │
│    Distance     [ Random]       │                                                       │
│      fixed dist [mm] -0.7,0.7   │  ← comma-separated; one run per value                 │
│    X position   [ Random]       │                                                       │
│      fixed x [cm]  [0.0,0.09]  │  ← Random unchecks → uniform over wire span           │
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
| **Waveforms** | Mean anode and cathode current waveforms overlaid per (distance, x-position) combination. A distance selector and (when fixed x-positions were simulated) an x-position dropdown choose the folder to display (a random distance/x-position appears as `—`). The **e⁻/ion components** checkbox overlays the separate electron and ion contributions to each induced current (requires a ROOT file with the component-split branches). The **Amplifier output [mV]** checkbox switches the traces to the front-end amplifier output voltage (requires a file produced with `amplifier.enable = true`). Read directly from the ROOT file via uproot |
| **Charge** | Cumulative charge integrals Q(t) — running integral of each waveform — for anode and cathode, per (distance, x-position) pair. An event slider selects individual events. ROOT TCanvas opens separately (PyROOT required) |
| **E-Field** | Interactive 2D electric field map in any of the XY, XZ, or YZ planes at a configurable depth; binning configurable from 50 to 10 000 bins per axis (PyROOT required) |
| **Weighting Field** | Exact Shockley–Ramo weighting field/potential of a selected electrode (`anode`, `cathode`, `cathode_top`), computed with Garfield's `ComponentAnalyticField` — the same geometry the simulation uses, so the wire screening is faithful. Quantity selectable (W potential, \|E_w\|, E_w,x, E_w,y); XY colour map plus X/Y profile slices. Interactive from the geometry spinboxes — no simulation run needed. In resistive mode the cathode map is α-scaled and the insulator region shows the weighting-potential ramp (α→1 up to the readout pad), with the resistive layer and pad marked; the time-domain exp(−t/τ) relaxation is not shown on the static map. Requires PyROOT **and** a loadable Garfield library |
| **3D Tracks** | Per-event 3D detector view in a ROOT TCanvas showing detector geometry and drift lines with correct aspect ratios. Controls: preset view buttons (Gap XY / Top XZ / Side YZ / 3D reset), zoom ± (down to 0.5 % of full range), pan X/Y/Z. Distance and x-position selectors mirror the simulated folder structure. Wires rendered as semi-transparent 12-sided tube wireframes at actual diameter (clipped to the visible frame); cathode planes clipped to visible cube. Primary electron and ion drift lines colour-coded (blue / green / magenta / grey) and semi-transparent (PyROOT required) |
| **Magboltz** | Gas transport-property viewer: reads the `_props.csv` sidecar file exported automatically alongside the `.gas` file. Displays electron drift velocity, Townsend α, attachment η, longitudinal and transverse diffusion coefficients, effective gain (α − η), ion drift velocity, and ion mobility as a function of E-field — all in one ROOT TCanvas (8 panels). Export buttons save the plots to a ROOT file or copy the CSV. Requires PyROOT |

**▶ Run** starts the simulation in a background thread (the window stays fully
responsive).  **■ Stop** sends SIGTERM.  **Load Config** / **Save Config** read and
write `.json` files that are fully compatible with the CLI `--config` flag.

> **Note:** The E-Field, Weighting Field, 3D Tracks, and Magboltz tabs open ROOT TCanvas
> windows and require PyROOT (ROOT importable from Python — available automatically when
> using the conda ROOT installation described in the Build section above). The Weighting
> Field tab additionally loads the Garfield shared library into the process; if it cannot
> be found the tab logs a message and is disabled (the rest of the GUI is unaffected).

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
| `enable_delayed_signal`        | bool   | —     | `true`          | When `true`, the exp(−t/τ) resistive relaxation is applied to the pad waveform as an exact exponential post-filter (negligible cost). When `false`, only the static α-corrected weighting potential is used — no relaxation tail. Ignored when `type = "conductive"` |
| `ground_plane_enabled`         | bool   | —     | `false`         | Add a grounded plane below the readout pad (resistive only). It adds a pad-to-ground capacitance `C_gnd` that lowers the weighting-potential factor α (`α = C_ins/(C_ins+C_gap+C_gnd)`), reducing the pad signal. Does not change the DC field or τ. No effect when `type = "conductive"` (the solid grounded cathode shields it — a note is printed) |
| `ground_plane_insulator_um`    | float  | μm    | 1000.0          | Thickness of the pad ↔ ground-plane insulator (default 1 mm). Used only when `ground_plane_enabled` |
| `ground_plane_insulator_material` | string | —  | `"fr4"`         | Dielectric of the pad ↔ ground-plane gap: `"kapton"` (ε_r=3.5), `"fr4"` (ε_r=4.6), or `"air"` (ε_r=1.0). Sets `C_gnd`, hence the signal reduction. Used only when `ground_plane_enabled` |

### `amplifier`

Opt-in front-end model (CIVIDEC C2-TCT) applied to the anode and cathode waveforms; see
Physics § 5. All keys are ignored when `enable = false`.

| Key                    | Type   | Unit | Default      | Description |
|------------------------|--------|------|--------------|-------------|
| `enable`               | bool   | —    | `false`      | When `true`, produce the amplifier output voltage [mV] for the anode (wire) and cathode (pad) channels as the `anode_amp` / `cathode_amp` branches and `p_anode_amp` / `p_cathode_amp` profiles. When `false`, those outputs are zero and the run is byte-for-byte identical to one built without this feature |
| `gain_db`              | float  | dB   | 40.0         | Voltage gain (40 dB = ×100). Sets the output scale V_out[mV] = 10^(gain_db/20)·R_in·i[µA]·10⁻³ |
| `input_impedance_ohm`  | float  | Ω    | 50.0         | Amplifier input impedance. Sets the AC-coupling high-pass τ = R_in·C and the output scale |
| `bandwidth_high_hz`    | float  | Hz   | 2.0e9        | Upper −3 dB band edge → input low-pass τ = 1/(2π·f). Smooths sub-ns features; only significant at very small `time_step_ns` |
| `bandwidth_low_hz`     | float  | Hz   | 1.0e4        | Lower −3 dB band edge → high-pass τ = 1/(2π·f) ≈ 15.9 µs. Slow baseline droop, visible only in long time windows |
| `coupling_cap_nf`      | float  | nF   | 1.0          | AC-coupling capacitor at the input. With R_in sets the pad high-pass τ (1 nF, 50 Ω → 50 ns) |
| `wire_series_cap_pf`   | float  | pF   | 470.0        | Extra series capacitor on the **anode (wire)** input only; in series with the coupling cap (470 pF ⊕ 1 nF = 320 pF) it gives the wire channel a shorter high-pass τ ≈ 16 ns |

### `source`

| Key                   | Type          | Unit | Default              | Description                                               |
|-----------------------|---------------|------|----------------------|-----------------------------------------------------------|
| `energy_keV`          | float         | keV  | 5.9                  | Photon energy of the simulated X-ray source (Fe-55 K-alpha line). Determines the number of primary electrons via `N = round(E_photon[eV] / w_value_eV)` (≈227 for Fe-55 at 5.9 keV and W = 26 eV) |
| `source_distances_mm` | float array or null | mm | [-0.7, 0.7] | List of signed y-positions (mm) at which primary electrons are placed, measured from the wire plane. Positive → readout cathode side (y < 0); negative → cathode_top side (y > 0). Each distance is a separate simulation run. Values are clamped to (−`gap_cm`×10, +`gap_cm`×10). `null` → a single run with the distance drawn uniformly over the gap each event (ROOT directory `dist_rnd/`) |
| `x_positions_cm`      | float array or null | cm | [0.0, 0.9] | Fixed lateral (x) positions [cm] for the photon interaction point. `null` → uniform random over the wire array each event; the shipped default specifies two fixed x-positions (wire centre and midpoint). One or more values → one simulation run per distance × x-position pair; ROOT directory named `dist_Nmm_xMmm/`. Backward-compatible with the old scalar `x_position_cm` key (wrapped into a one-element list) |

### `gas`

| Key                     | Type   | Unit  | Default                   | Description                                      |
|-------------------------|--------|-------|---------------------------|--------------------------------------------------|
| `gas1`                  | string | —     | `"ar"`                    | First gas species (Magboltz name, lowercase). Passed to the `MediumMagboltz(gas1, frac1, gas2, 100−frac1)` constructor and used as the `.gas` filename prefix |
| `gas1_fraction_pct`     | float  | %     | 70.0                      | Volume fraction of `gas1`; `gas2` gets the remaining `100 − gas1_fraction_pct`. Validated to lie in (0, 100) |
| `gas2`                  | string | —     | `"co2"`                   | Second (complementary) gas species |
| `ion_species`           | string | —     | `"co2"`                   | Which single-component ion-mobility table to load: `IonMobility_{X}+_{X}.txt` from the Garfield++ data directory (X uppercased). Should match the dominant drifting ion of the mixture. Files ship for ar, co2, cf4, he, ne |
| `temperature_K`         | float  | K     | 293.15                    | Gas temperature passed to Magboltz. Affects gas number density (n ∝ 1/T at fixed pressure), which shifts drift velocity and Townsend coefficients. Must match physical detector conditions; 293.15 K = 20 °C |
| `pressure_Torr`         | float  | Torr  | 750.0                     | Gas pressure passed to Magboltz. Together with temperature, sets gas density. 750 Torr ≈ 0.987 atm. Reducing pressure increases mean free path and electron energy, raising gain |
| `enable_penning`        | bool   | —     | true                      | Activates Penning transfer via `MediumMagboltz::EnablePenningTransfer()`. In Ar:CO2 70:30 this raises the effective Townsend coefficient by ~20–40 %, bringing simulated gain closer to measured values. Should be left on for this mixture |
| `n_magboltz_collisions` | int    | —     | 10                        | Monte Carlo collision cycles Magboltz runs per E-field grid point. Higher values reduce statistical uncertainty in transport coefficients at the cost of longer gas-file generation. 2–5: smoke test; 10: default; 20–50: publication quality |
| `max_electron_energy_eV`| float  | eV    | 2000.0                    | Upper energy bound for the Magboltz cross-section look-up table. Must exceed the maximum kinetic energy electrons reach near the wire (typically 500–1000 eV here). Too low a ceiling causes Magboltz to extrapolate, producing unphysical transport coefficients |
| `n_field_points`        | int    | —     | 50                        | Number of logarithmically spaced E-field values at which Magboltz computes transport coefficients (from `e_field_min_vcm` to `e_field_max_vcm`). More points give smoother interpolation; fewer points speed up gas-file generation |
| `e_field_min_vcm`       | float  | V/cm  | 100.0                     | Lower E-field limit for the Magboltz transport table. Sets the gentle-drift end of the logarithmic grid; encoded in the `.gas` filename as `Ef{min}v-…`. Lower it only if the gas-gap field can fall below 100 V/cm |
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
| `ion_max_step_um`    | float | μm   | 5.0      | Cap on the `DriftLineRKF` integration step.  The stepper's steps otherwise grow geometrically (×10 per step) and the induced current is sampled only at drift-line points, so an uncapped surface-born ion's ~10 μm step spans ~5–8 ns right where i(t) varies fastest — producing an artificial flat shelf with a sharp kink ~8 ns after the electron spike.  5 μm resolves the early ion signal to ~1–2 ns at roughly 5× the ion-drift CPU; `0` disables the cap |
| `random_seed`        | int   | —    | 0        | Seed for the random-number generators.  Seeds **both** ROOT's `gRandom` (which drives source-position sampling) and Garfield's own transport engine (the `TRandom3` behind `AvalancheMicroscopic` / `DriftLineRKF`, installed via `Garfield::Random::SetEngine`) — both are required for reproducibility, since Garfield does not draw from `gRandom`.  `0` (default) self-seeds both randomly each run; any positive integer fixes the avalanche, ion-drift, and source-position sequence so runs are bit-for-bit reproducible — useful for A/B comparisons and debugging |

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
| `source_distances_mm: null` (random distance) | `dist_rnd/` (or `dist_rnd_x0mm/` with fixed x) |

The x-position suffix uses millimetres (× 10 relative to the cm config value) with
decimal points replaced by `p` (e.g. 0.18 cm → 1.8 mm → `x1p8mm`).  A negative
distance uses an `m` prefix (e.g. −0.7 mm → `dist_m0p7mm`); a random distance uses
the literal `dist_rnd`.

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
| `p_anode_electron`      | TProfile | Electron-only component of the anode current vs time [fC/ns] (prompt + delayed) |
| `p_anode_ion`           | TProfile | Ion-only component of the anode current vs time [fC/ns] |
| `p_cathode_electron`    | TProfile | Electron-only component of the cathode current vs time [fC/ns] |
| `p_cathode_ion`         | TProfile | Ion-only component of the cathode current vs time [fC/ns] |
| `p_anode_amp`           | TProfile | Mean anode amplifier output vs time [mV] (zero unless `amplifier.enable`) |
| `p_cathode_amp`         | TProfile | Mean cathode amplifier output vs time [mV] (zero unless `amplifier.enable`) |

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
| `anode_e` / `anode_i` | vector\<float\> | Electron / ion component of the anode current [fC/ns]; the two sum to `anode` bin-by-bin |
| `cathode_e` / `cathode_i` | vector\<float\> | Electron / ion component of the cathode current [fC/ns]; the two sum to `cathode` |
| `anode_amp` / `cathode_amp` | vector\<float\> | Amplifier output voltage [mV] for the anode (wire) and cathode (pad) channels; present in every file but zero unless `amplifier.enable` (see Physics § 5) |
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
| `source_distance_mm` | mm | Signed y-distance of the photon interaction from the wire plane. Positive → readout cathode side (y < 0); negative → cathode_top side (y > 0). See sign convention in *Detector geometry*. The literal `random` when `source_distances_mm: null` (distance drawn per event) |
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
| `wfield`    | Weighting Field | Per-electrode weighting field/potential maps |

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

**Electron spike width.**  The simulated electron spike is extremely narrow
(FWHM ≈ 0.2 ns — the intrinsic duration of a single avalanche, which multiplies within
the last few wire radii where electrons move at >100 µm/ns).  Real measured pulses are
~3 ns wide.  The difference is the W-value point-deposit approximation: all N ≈ 227
primaries start at one point and are represented by one scaled avalanche, so they share
a single arrival time.  In reality the primaries are distributed along the
photoelectron + Auger track (~150–300 µm for 5.9 keV in Ar → 1.8–3.6 ns arrival spread
at the ~83 µm/ns gap drift velocity) and each diffuses independently
(σ_L ≈ 40 µm over a 0.7 mm drift → ~0.5 ns), and the readout electronics adds ~1 ns.
The approximation is exact for the charge observables (Q_anode, Q_cathode, ratio) but
discards the per-event time structure of the electron peak.

**Why the ion tail opens with a plateau.**  Although the near-wire field is radial
(E ∝ 1/r), the *induced* current scales as 1/r²: the ion drift velocity (v = μE) and
the wire's Ramo weighting field both fall as 1/r.  Solving the ion motion in the 1/r
field gives r²(t) = r₀²(1 + t/t₀), so

> i(t) = i₀ / (1 + t/t₀),   with   t₀ = r₀² / (2 μ k)

where r₀ is the ion's birth radius (≈ the 25 μm wire radius) and k = E·r is the field
constant.  Because the induced current depends on r², which barely changes while the ion
is still near r₀, i(t) is flat — a plateau — for t < t₀, then decays as 1/t once r²
grows.  For ions born at the wire surface t₀ = r₀²/(2μk) ≈ 4–7 ns; a fit over a wider
window returns a larger effective t₀ (~13 ns) because ions are born over a range of
radii.  This is the differential form of the classic logarithmic ion charge
Q(t) ∝ ln(1 + t/t₀).  `tools/plot_ion_tail.py` fits i₀/(1+t/t₀) to a run and overlays
it in log-log and linear scale.  Beyond the near-wire region the current falls faster
than 1/t as ions cross into the more uniform bulk field.

**Drift-line discretization artifact (and the `ion_max_step_um` cap).**  Without a step
cap, the first ~8 ns after the spike show an *artificially* flat shelf (constant to ~5
digits) ending in a sharp kink: `DriftLineRKF` grows its integration steps geometrically
(×10 per step) and the induced current is sampled only at the drift-line points, so a
surface-born ion's ~10 μm step spans ~5–8 ns — exactly where i(t) should already have
fallen by ~half.  All avalanche ions are born at nearly the same spot, so their step
boundaries kink in lockstep instead of averaging out.  The `simulation.ion_max_step_um`
cap (default 5 μm) shortens these segments so the early rollover is genuinely resolved;
set it to 0 to recover the old (faster, under-sampled) behaviour.

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

**The fast electron spike on the cathode is physical, not a bug.**  It is tempting to
expect no electron component on the cathode because the cathode weighting potential is
small near the wire (the wire screens it — see the **Weighting Field** tab).  But what
sets the induced charge is the *change* in weighting potential along the electron's path,
ΔW_cathode, not its absolute value.  Just outside the wire surface (~30 μm) the weighting
potentials are W_anode ≈ 0.36 and W_cathode ≈ W_cathode_top ≈ 0.32 (they sum to 1).  By the
up–down symmetry of the two equidistant cathode planes, the avalanche electrons collapsing
onto the wire split their induced charge almost equally between the two cathodes: each gets
Wc/(1−Wa) ≈ 0.5 of the anode electron signal (the simulation gives cathode_e/anode_e ≈
0.50).  This is genuine Shockley–Ramo induction — the three electrodes satisfy
Q_anode + Q_cathode + Q_cathode_top ≈ 0 to machine precision.  The spike carries only ~8 %
of the cathode *charge* but dominates the peak *current* because it is sub-nanosecond, so a
real cathode readout usually does not resolve it (limited amplifier bandwidth; and the
point-deposit approximation makes the simulated spike artificially narrow, ~0.2 ns vs ~3 ns
in reality — see *Anode signal shape*).  The per-electrode electron/ion split
(`anode_e`/`anode_i`/`cathode_e`/`cathode_i` branches and the Waveforms-tab **e⁻/ion
components** overlay) lets you inspect this directly.

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
