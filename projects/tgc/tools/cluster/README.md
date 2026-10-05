# Gain-vs-position scan — running it (local cores or an HPC cluster)

Maps the **single-electron avalanche gain** as a function of the primary-electron
position (x across the wire pitch, y depth in the gap), using `tgc_sim`'s
microscopic avalanche. The config [`config/scan_gain.json`](../../config/scan_gain.json)
turns ion drift **off** and sets `energy_keV=0.026` (→ `nPrimary=1`), so
`mean_avalanche_size` in `summary.csv` is the single-electron gain and the
`t_signals` tree's `avalanche_size` branch holds the per-event gain.

> **Why not a fast DriftLineRKF Townsend-integral map?** It was tried and diverges for
> this thin-wire geometry (gain ~1e38 vs the microscopic ~4e4) — the line integral
> doesn't model collection at the wire. The microscopic avalanche is the correct path.

**Cost:** the microscopic avalanche is **~10 s/event** here (gain ~4e4). Budget
accordingly: `n_points × n_events × 10 s`, divided by the cores/tasks you run.
Gain is ion-drift-independent, so **no `GARFIELD_INSTALL` / ion-mobility table is
needed** for the scan — only the Garfield/ROOT runtime libraries.

## Build (once, on the cluster)

```bash
cd projects/tgc
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$GARFIELD_INSTALL;$ROOTSYS"      # + -DVDT_* if ROOT needs Vdt
cmake --build build -j
```
Ensure the gas table is present in `projects/tgc/` (it is committed:
`ar70_co2_30_T293_P760_Ee2000_Ef100v-400k_n10_c2_pen.gas`). If absent, generate it
once (serially) before launching a parallel scan, so jobs don't all regenerate it:
```bash
./build/tgc_sim --config config/scan_gain.json --distance 0.5 --run-name warmup --out results/gain_scan
```

## Option 1 — one multi-core machine (no scheduler)

```bash
python3 tools/gain_scan.py run --config config/scan_gain.json   # one process per depth, across cores
python3 tools/gain_scan.py plot --csv results/gain_scan/gain_map.csv --rms
```
`run` launches one `tgc_sim` per depth (`--distance`) with the shared x-list and
merges the `summary.csv` rows. For finer control, edit the x/y grids in the config.

## Option 2 — SLURM cluster (array job)

```bash
# 1. Emit per-job configs + a manifest (prints the array size N). --xchunks splits
#    the x-list for more, smaller tasks (good when you have many cores):
python3 tools/gain_scan.py emit --config config/scan_gain.json \
        --out $SCRATCH/gain_scan --xchunks 4

# 2. Submit the array (N from step 1). Edit submit_gain_scan.sbatch first:
#    set TGC_REPO / the module/source line / account / partition.
sbatch --array=0-$((N-1)) tools/cluster/submit_gain_scan.sbatch \
       $SCRATCH/gain_scan/manifest.txt $SCRATCH/gain_scan

# 3. After the array finishes, merge + plot (login node):
python3 tools/gain_scan.py merge --out $SCRATCH/gain_scan
python3 tools/gain_scan.py plot  --csv  $SCRATCH/gain_scan/gain_map.csv --rms
```

Each array task `cd`s to the repo (for the gas cache), reads its manifest line, and
runs `tgc_sim --config <job.json> --out $SCRATCH/gain_scan --run-name <tag>`. Jobs
never clobber each other (distinct run-names); results merge by the self-describing
`source_distance_mm` / `x_position_cm` columns.

## Other schedulers / no SLURM

The `emit` manifest is scheduler-agnostic (`<abs config path> <run-name>` per line).
- **PBS/HTCondor:** mirror `submit_gain_scan.sbatch` — index the manifest by the array
  index and run the same `tgc_sim` line.
- **Plain parallel on one node:**
  ```bash
  OUT=$SCRATCH/gain_scan
  awk '{print $1, $2}' $OUT/manifest.txt | \
    xargs -P "$(nproc)" -L1 sh -c 'cd <repo>; ./build/tgc_sim --config "$1" --out '"$OUT"' --run-name "$2"' _
  python3 tools/gain_scan.py merge --out $OUT
  ```

## Outputs
- `summary.csv` per run-dir and merged `gain_map.csv` — `mean_avalanche_size` per (x,y).
- `tgc_sim.root` per run-dir — the `t_signals` tree's `avalanche_size` branch gives the
  full per-event gain distribution (Polya) at each point; `plot --rms` reports its width.
- `gain_map.png` — 2D heatmap + gain-vs-x and gain-vs-depth slices.
