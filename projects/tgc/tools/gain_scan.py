#!/usr/bin/env python3
"""
Drive and plot a 2D avalanche-gain scan over the primary-electron position
(x transverse across the wire pitch, y depth in the gap), using tgc_sim's
MICROSCOPIC avalanche (the physically correct gain for this thin-wire geometry;
the DriftLineRKF Townsend-integral shortcut was found to diverge here).

Facts it relies on (see src/tgc_sim.cc):
  * one run does a full CROSS PRODUCT of source_distances_mm x x_positions_cm
    (one summary.csv row per (x, y) cell);
  * config/scan_gain.json turns ion drift OFF and sets energy_keV=0.026 so
    nPrimary=1 and `mean_avalanche_size` IS the single-electron gain;
  * the binary is single-threaded, so we parallelise by launching one process
    per depth (--distance <y>) and merging the self-describing summary.csv rows.

Cost: the microscopic avalanche is ~10 s/event here (gain ~4e4), so a fine map is
many CPU-hours -> parallelise across cores (`run`) or export to a cluster
(`emit` + tools/cluster/submit_gain_scan.sbatch).

Subcommands
-----------
    run        --config C [--out DIR] [--jobs N]  launch locally in parallel, then pool
    emit       --config C --out DIR [--xchunks K] write per-job configs + manifest (scheduler)
    accumulate --out DIR                          pool per-event gains of all batches -> gain_map.csv
    plot       [--csv F] [--config C] [--rms]     2D gain heatmap + 1D slices from the pooled CSV
    waveform   --out DIR [--tcut 100]             mean peak-normalized anode current shape vs
                                                  position (needs an ion-drift run, e.g.
                                                  config/scan_waveform.json); T_cut selectable here
"""

import argparse
import concurrent.futures as cf
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).parent.resolve()
TGC_DIR = (SCRIPT_DIR / "..").resolve()
BINARY = TGC_DIR / "build" / "tgc_sim"
DEFAULT_OUT = TGC_DIR / "results" / "gain_scan"


def _ytag(y: float) -> str:
    return "y_" + f"{y:g}".replace(".", "p").replace("-", "m")


def _derive_gas_filename(gas: dict) -> str:
    """Mirror DeriveGasFileName in src/tgc_sim.cc (for the cache-hit pre-check)."""
    f1 = round(gas["gas1_fraction_pct"])
    return (f"{gas['gas1']}{f1}_{gas['gas2']}_{100 - f1}"
            f"_T{round(gas['temperature_K'])}_P{round(gas['pressure_Torr'])}"
            f"_Ee{round(gas['max_electron_energy_eV'])}"
            f"_Ef{round(gas['e_field_min_vcm'])}v-{round(gas['e_field_max_vcm'] / 1000)}k"
            f"_n{gas['n_field_points']}_c{gas['n_magboltz_collisions']}"
            f"_{'pen' if gas.get('enable_penning', True) else 'nopen'}.gas")


def _wire_positions_cm(geom: dict) -> np.ndarray:
    n, pitch = geom["n_wires"], geom["wire_pitch_cm"]
    return (np.arange(n) - (n - 1) / 2.0) * pitch


def _warn_ion_drift(cfg: dict) -> None:
    """Ion drift is slow and needs the CO2+ mobility table via GARFIELD_INSTALL."""
    if not cfg["simulation"].get("enable_ion_drift", False):
        return
    print("note: enable_ion_drift is true — the waveform includes the ion tail (needed "
          "for the waveform test) but each event is slower.", file=sys.stderr)
    if not os.environ.get("GARFIELD_INSTALL"):
        print("ERROR-RISK: GARFIELD_INSTALL is not set; tgc_sim will abort (no ion-mobility "
              "table). Export GARFIELD_INSTALL=<...>/local/garfield before running.",
              file=sys.stderr)


def _check_gas(cfg: dict, config_name: str) -> bool:
    gas_file = TGC_DIR / _derive_gas_filename(cfg["gas"])
    if gas_file.exists():
        return True
    y0 = cfg["source"]["source_distances_mm"][0]
    print(f"error: gas table {gas_file.name} not found at the project root.\n"
          f"       Generate it once first:\n"
          f"       cd {TGC_DIR} && ./build/tgc_sim --config {config_name} "
          f"--distance {y0} --run-name warmup --out results/gain_scan",
          file=sys.stderr)
    return False


def _parse_tag(tag: str):
    """Invert FileSafeNumber: 'dist_0p1mm_x0p3mm' -> (dist_mm=0.1, x_cm=0.03).

    'p'->'.', 'm'->'-'; the tag's x is in mm (tgc_sim stores x*10), so /10 for cm.
    Returns (None, None) for un-parseable tags (e.g. random 'dist_rnd').
    """
    m = re.match(r"dist_(.+?)mm(?:_x(.+?)mm)?$", tag)
    if not m:
        return None, None

    def num(s):
        if s is None:
            return None
        try:
            return float(s.replace("p", ".").replace("m", "-"))
        except ValueError:
            return None

    dist = num(m.group(1))
    x_mm = num(m.group(2))
    return dist, (x_mm / 10.0 if x_mm is not None else None)


def _pool_events(out_dir: Path) -> Path:
    """Pool the per-event `avalanche_size` of EVERY tgc_sim.root under out_dir,
    grouped by point tag, into out_dir/gain_map.csv — the incremental, no-waste
    aggregation. Running more batches under out_dir and re-pooling only adds events.
    """
    import uproot
    from collections import defaultdict
    pools = defaultdict(list)
    for root in sorted(Path(out_dir).glob("**/tgc_sim.root")):
        try:
            with uproot.open(root) as f:
                for name in f.keys(cycle=False):
                    if name.rsplit("/", 1)[-1] != "t_signals":
                        continue
                    arr = f[name]["avalanche_size"].array(library="np")
                    if len(arr):
                        pools[name.rsplit("/", 1)[0]].append(arr)
        except Exception as e:  # noqa: BLE001 — skip a corrupt/partial file, keep going
            print(f"warning: skipping {root}: {e}", file=sys.stderr)
    if not pools:
        raise SystemExit(f"no t_signals/avalanche_size found under {out_dir}")

    rows = []
    for tag, arrs in pools.items():
        dist, x = _parse_tag(tag)
        if dist is None or x is None:
            continue
        g = np.concatenate(arrs)
        n = int(len(g))
        rms = float(g.std())
        rows.append({
            "source_distance_mm": dist, "x_position_cm": x, "n_events": n,
            "mean_avalanche_size": float(g.mean()),
            "rms_avalanche_size": rms,
            "sem_avalanche_size": rms / np.sqrt(n) if n else float("nan"),
        })
    df = (pd.DataFrame(rows)
          .sort_values(["source_distance_mm", "x_position_cm"]).reset_index(drop=True))
    path = Path(out_dir) / "gain_map.csv"
    df.to_csv(path, index=False)
    return path


# ─────────────────────────────── run (local) ────────────────────────────────

def run(args) -> int:
    if not BINARY.exists():
        print(f"error: {BINARY} not built.", file=sys.stderr)
        return 1
    cfg = json.loads(Path(TGC_DIR / args.config if not Path(args.config).is_absolute()
                          else args.config).read_text())
    _warn_ion_drift(cfg)
    if not _check_gas(cfg, args.config):
        return 1

    ys = cfg["source"]["source_distances_mm"]
    nx, nev = len(cfg["source"]["x_positions_cm"]), cfg["simulation"]["n_events"]
    out_dir = Path(args.out)
    label = args.label or datetime.now().strftime("batch_%Y%m%d_%H%M%S")
    batch_dir = out_dir / label
    batch_dir.mkdir(parents=True, exist_ok=True)
    jobs = args.jobs or len(ys)
    print(f"Batch '{label}': {len(ys)} depth jobs ({nx} x-points x {nev} events each) "
          f"on up to {jobs} workers -> {batch_dir}")

    def _one(y: float):
        cmd = [str(BINARY), "--config", args.config, "--distance", f"{y:g}",
               "--out", str(batch_dir), "--run-name", _ytag(y)]
        t0 = time.time()
        p = subprocess.run(cmd, cwd=TGC_DIR, capture_output=True, text=True)
        return y, p.returncode, time.time() - t0, p.stderr[-400:]

    t_start = time.time()
    rcs = []
    with cf.ThreadPoolExecutor(max_workers=jobs) as ex:
        for y, rc, dt, err in ex.map(_one, ys):
            print(f"  [{'ok ' if rc == 0 else 'FAIL'}] y={y:>5} mm  ({dt:5.1f} s)"
                  + ("" if rc == 0 else f"  {err}"))
            rcs.append(rc)
    if any(rc != 0 for rc in rcs):
        print("error: one or more jobs failed.", file=sys.stderr)
        return 1
    path = _pool_events(out_dir)
    print(f"\nDone in {time.time() - t_start:.1f} s -> {path}  "
          f"(pooled over every batch under {out_dir})")
    print(f"Plot with:  python3 tools/gain_scan.py plot --csv {path}")
    return 0


# ─────────────────────────── emit (for a scheduler) ─────────────────────────

def emit(args) -> int:
    cfg = json.loads(Path(TGC_DIR / args.config if not Path(args.config).is_absolute()
                          else args.config).read_text())
    _warn_ion_drift(cfg)
    if not _check_gas(cfg, args.config):
        return 1
    out_dir = Path(args.out)
    # Each batch is self-contained under out_dir/<label>/; accumulate pools over out_dir.
    label = args.label or datetime.now().strftime("batch_%Y%m%d_%H%M%S")
    batch_dir = out_dir / label
    jobs_dir = batch_dir / "jobs"
    jobs_dir.mkdir(parents=True, exist_ok=True)

    ys = cfg["source"]["source_distances_mm"]
    xs = cfg["source"]["x_positions_cm"]
    xchunks = max(1, args.xchunks)
    x_splits = [list(c) for c in np.array_split(xs, min(xchunks, len(xs)))]

    manifest = []
    for y in ys:
        for ci, xsub in enumerate(x_splits):
            sub = json.loads(json.dumps(cfg))           # deep copy
            sub["source"]["source_distances_mm"] = [y]
            sub["source"]["x_positions_cm"] = xsub
            tag = _ytag(y) + (f"_xc{ci}" if xchunks > 1 else "")
            cfg_path = jobs_dir / f"{tag}.json"
            cfg_path.write_text(json.dumps(sub, indent=2))
            # manifest line: <absolute config path> <run-name>; the sbatch cd's to
            # the repo root (for the gas cache) and passes these to tgc_sim.
            manifest.append(f"{cfg_path.resolve()} {tag}")

    (batch_dir / "manifest.txt").write_text("\n".join(manifest) + "\n")
    print(f"Emitted {len(manifest)} jobs -> {jobs_dir}\n"
          f"Manifest: {batch_dir / 'manifest.txt'}  (array size = {len(manifest)})\n"
          f"Run tgc_sim with --out {batch_dir} (sbatch/xargs), then\n"
          f"`gain_scan.py accumulate --out {out_dir}` to pool all batches.")
    return 0


def accumulate(args) -> int:
    """Pool every batch under --out into a combined gain_map.csv (incremental, no waste)."""
    path = _pool_events(Path(args.out))
    df = pd.read_csv(path)
    print(f"pooled {len(df)} (x,y) points, {int(df['n_events'].sum())} events total "
          f"({int(df['n_events'].min())}-{int(df['n_events'].max())} per point) -> {path}")
    return 0


# ─────────────────────────────────── plot ───────────────────────────────────

def _point_rms(root_path: Path):
    """Yield (point_name, mean, rms) of per-event gain from each point's t_signals tree.

    Reads the `avalanche_size` branch (not the h_avalanche_size histogram, whose
    ROOT auto-bin buffer is not flushed below ~1000 entries and reads as empty).
    """
    import uproot
    with uproot.open(root_path) as f:
        for name in f.keys(cycle=False):
            if name.rsplit("/", 1)[-1] != "t_signals":
                continue
            arr = f[name]["avalanche_size"].array(library="np")
            if len(arr) == 0:
                continue
            yield name.rsplit("/", 1)[0], float(arr.mean()), float(arr.std())


def plot(args) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    csv = Path(args.csv)
    if not csv.exists():
        print(f"error: {csv} not found. Run the scan first.", file=sys.stderr)
        return 1
    df = pd.read_csv(csv)
    cfg = json.loads((TGC_DIR / args.config).read_text())
    wires = _wire_positions_cm(cfg["geometry"])

    piv = df.pivot_table(index="source_distance_mm", columns="x_position_cm",
                         values="mean_avalanche_size")
    xs, ys, Z = piv.columns.to_numpy(), piv.index.to_numpy(), piv.to_numpy()

    fig, (axm, axx, axy) = plt.subplots(1, 3, figsize=(16, 5))
    im = axm.pcolormesh(xs * 10.0, ys, Z, shading="nearest", cmap="viridis")
    for xw in wires:
        if xs.min() - 1e-9 <= xw <= xs.max() + 1e-9:
            axm.axvline(xw * 10.0, color="white", lw=0.8, ls="--", alpha=0.7)
    axm.set_xlabel("x  [mm]   (dashed = wire)")
    axm.set_ylabel("source depth from wire plane  [mm]")
    axm.set_title("Single-electron gain  vs (x, depth)")
    fig.colorbar(im, ax=axm, label="mean avalanche size  (gain)")

    def dist_from_wire_mm(xc):
        return np.min(np.abs(np.asarray(xc)[:, None] - wires[None, :]), axis=1) * 10.0
    dfw = dist_from_wire_mm(xs)
    order = np.argsort(dfw)
    for yv in ys[:: max(1, len(ys) // 4)]:
        axx.plot(dfw[order], piv.loc[yv].to_numpy()[order], "o-", ms=3, label=f"y={yv:g} mm")
    axx.set_xlabel("distance from nearest wire  [mm]")
    axx.set_ylabel("mean avalanche size  (gain)")
    axx.set_title("Gain vs transverse position")
    axx.legend(fontsize=8)

    for xv in xs[:: max(1, len(xs) // 4)]:
        axy.plot(ys, piv[xv].to_numpy(), "o-", ms=3,
                 label=f"x={xv*10:g} mm (d_w={dist_from_wire_mm([xv])[0]:.2g})")
    axy.set_xlabel("source depth from wire plane  [mm]")
    axy.set_ylabel("mean avalanche size  (gain)")
    axy.set_title("Gain vs depth")
    axy.legend(fontsize=8)

    fig.tight_layout()
    out_png = csv.with_name("gain_map.png")
    fig.savefig(out_png, dpi=130)
    print(f"wrote {out_png}")

    if args.rms:
        cap = cfg["simulation"]["max_avalanche_size"]
        if "rms_avalanche_size" in df.columns:  # pooled CSV — use it directly
            print(f"\nPer-point gain (pooled over all batches); cap={cap}:")
            for _, r in df.iterrows():
                flag = "  <-- near cap!" if r.mean_avalanche_size > 0.8 * cap else ""
                print(f"  dist={r.source_distance_mm:>5g} mm  x={r.x_position_cm*10:>4.1f} mm  "
                      f"N={int(r.n_events):>4}  mean={r.mean_avalanche_size:10.1f}  "
                      f"rms={r.rms_avalanche_size:10.1f}  sem={r.sem_avalanche_size:8.1f}{flag}")
        else:  # legacy CSV without rms columns — re-read the trees
            print(f"\nPer-point gain spread (per-event, from t_signals); cap={cap}:")
            for d in sorted(Path(args.csv).parent.glob("**/tgc_sim.root")):
                for pt, mean, rms in _point_rms(d):
                    flag = "  <-- near cap!" if mean > 0.8 * cap else ""
                    print(f"  {d.parent.name}/{pt:>18}  mean={mean:10.1f}  rms={rms:10.1f}{flag}")
    return 0


# ─────────────────── waveform shape test (position scan) ─────────────────────

def waveform(args) -> int:
    """Mean PEAK-ALIGNED, peak-normalized anode pulse shape (electron + ion) vs position.

    For every event (t_signals `anode` branch) under --out: find the current peak (the
    electron spike), normalize so the peak is +1, and TIME-ALIGN it so the spike sits at
    t'=0, then average per position. Aligning removes the per-event drift-time jitter
    (which otherwise smears the mean peak below 1) and the position-dependent drift
    offset, so the mean shape peaks at ~1 and positions overlay — differences then live
    in the ion tail. The spike time itself (= electron drift-time to the wire) is kept as
    a separate observable (waveform_drift.csv + map). Pools all batches (incremental);
    --tcut (ns of tail kept after the spike) is chosen here, no re-simulation. Needs an
    ion-drift run (config/scan_waveform.json) for the ion tail.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import uproot

    cfg = json.loads((TGC_DIR / args.config).read_text())
    dt = cfg["simulation"]["time_step_ns"]
    wires = _wire_positions_cm(cfg["geometry"])
    out_dir = Path(args.out)
    n_pre = int(round(5.0 / dt))                             # bins shown before the spike
    n_post = int(np.floor(args.tcut / dt - 0.5)) + 1        # bins of tail kept after it
    if n_post < 2:
        raise SystemExit(f"--tcut {args.tcut} too small for time_step {dt} ns")
    alen = n_pre + n_post
    search = int(round(50.0 / dt))                          # spike is within ~50 ns (max drift)

    asum, acnt = {}, {}          # per-bin aligned sum / count (edge-safe)
    dsum, dsq, nev = {}, {}, {}  # drift-time (spike-time) accumulators
    for root in sorted(out_dir.glob("**/tgc_sim.root")):
        try:
            with uproot.open(root) as f:
                for name in f.keys(cycle=False):
                    if name.rsplit("/", 1)[-1] != "t_signals":
                        continue
                    tag = name.rsplit("/", 1)[0]
                    if _parse_tag(tag)[0] is None:
                        continue
                    M = np.asarray(f[name]["anode"].array().tolist(), dtype=float)
                    if M.ndim != 2 or M.shape[0] == 0:
                        continue
                    nb = M.shape[1]
                    pk = np.argmax(np.abs(M[:, :min(search, nb)]), axis=1)  # spike bin / event
                    peak = M[np.arange(M.shape[0]), pk]
                    good = np.abs(peak) > 1e-12
                    if not good.any():
                        continue
                    a = asum.setdefault(tag, np.zeros(alen))
                    c = acnt.setdefault(tag, np.zeros(alen))
                    for e in np.nonzero(good)[0]:               # align each event to its spike
                        p = int(pk[e]); lo, hi = p - n_pre, p + n_post
                        sa, sb = max(0, lo), min(nb, hi)        # clip to the stored window
                        da = sa - lo
                        a[da:da + (sb - sa)] += M[e, sa:sb] / peak[e]
                        c[da:da + (sb - sa)] += 1.0
                    tpk = (pk[good] + 0.5) * dt                 # spike time = drift time to wire
                    dsum[tag] = dsum.get(tag, 0.0) + float(tpk.sum())
                    dsq[tag] = dsq.get(tag, 0.0) + float((tpk ** 2).sum())
                    nev[tag] = nev.get(tag, 0) + int(good.sum())
        except Exception as e:  # noqa: BLE001
            print(f"warning: skipping {root}: {e}", file=sys.stderr)
    if not nev:
        raise SystemExit(f"no t_signals/anode waveforms found under {out_dir}")

    t = (np.arange(alen) - n_pre) * dt                       # time since spike [ns] (peak at 0)
    means = {tag: np.divide(asum[tag], acnt[tag],
                            out=np.full(alen, np.nan), where=acnt[tag] > 0) for tag in asum}
    drift = {tag: (dsum[tag] / nev[tag],
                   float(np.sqrt(max(0.0, dsq[tag] / nev[tag] - (dsum[tag] / nev[tag]) ** 2))))
             for tag in nev}

    rows = []
    for tag, mean in means.items():
        d, x = _parse_tag(tag)
        for j in range(alen):
            rows.append({"source_distance_mm": d, "x_position_cm": x,
                         "time_since_peak_ns": t[j], "mean_norm_current": mean[j],
                         "n_events": nev[tag]})
    csv_path = out_dir / "waveform_shapes.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    drift_path = out_dir / "waveform_drift.csv"
    pd.DataFrame([{"source_distance_mm": _parse_tag(tag)[0], "x_position_cm": _parse_tag(tag)[1],
                   "mean_drift_ns": drift[tag][0], "rms_drift_ns": drift[tag][1],
                   "n_events": nev[tag]} for tag in drift]).to_csv(drift_path, index=False)

    def _dwire_mm(xc):  # distance from the nearest wire [mm], grid-agnostic
        return float(np.min(np.abs(xc - wires))) * 10.0

    def _map(ax, valfn, title, cbar, cmap):
        d = pd.DataFrame([{"x_position_cm": _parse_tag(tag)[1],
                           "source_distance_mm": _parse_tag(tag)[0],
                           "v": valfn(tag)} for tag in means])
        piv = d.pivot_table(index="source_distance_mm", columns="x_position_cm", values="v")
        im = ax.pcolormesh(piv.columns.to_numpy() * 10, piv.index.to_numpy(), piv.to_numpy(),
                           shading="nearest", cmap=cmap)
        for xw in wires:
            if piv.columns.min() <= xw <= piv.columns.max():
                ax.axvline(xw * 10, color="white", lw=0.8, ls="--", alpha=0.6)
        ax.set_xlabel("x [mm]  (dashed = wire)")
        ax.set_ylabel("source depth [mm]")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label=cbar)

    by_xy = {(round(_parse_tag(tag)[1], 4), round(_parse_tag(tag)[0], 4)): (means[tag], nev[tag])
             for tag in means}
    fig, (axo, axt, axd) = plt.subplots(1, 3, figsize=(18, 5))
    ys_all = sorted({round(_parse_tag(tag)[0], 4) for tag in means})
    xs_all = sorted({round(_parse_tag(tag)[1], 4) for tag in means})
    x_on, x_mid = min(xs_all, key=_dwire_mm), max(xs_all, key=_dwire_mm)
    for (xv, yv) in [(x_on, ys_all[0]), (x_on, ys_all[-1]),
                     (x_mid, ys_all[0]), (x_mid, ys_all[-1])]:
        hit = by_xy.get((round(xv, 4), round(yv, 4)))
        if hit:
            mean, n = hit
            axo.plot(t, mean, lw=1.3,
                     label=f"x={xv*10:g}mm (d_wire={_dwire_mm(xv):.2g}mm), y={yv:g}mm, N={n}")
    axo.axhline(0, color="k", lw=0.5)
    axo.axvline(0, color="0.6", lw=0.5, ls=":")
    axo.set_xlabel("time since spike [ns]")
    axo.set_ylabel("mean peak-normalized anode current")
    axo.set_title(f"Peak-aligned anode pulse shape  (T_cut={args.tcut:g} ns)")
    axo.legend(fontsize=8)

    tailmask = t > 5.0
    _map(axt, lambda tg: float(np.nanmean(means[tg][tailmask])),
         "Ion-tail height  (mean norm. current, >5 ns after spike)", "tail / peak", "magma")
    _map(axd, lambda tg: drift[tg][0],
         "Electron drift time to wire  (mean spike time)", "drift time [ns]", "viridis")

    fig.tight_layout()
    png = out_dir / "waveform_shapes.png"
    fig.savefig(png, dpi=130)
    print(f"pooled {len(nev)} positions, {sum(nev.values())} events (T_cut={args.tcut:g} ns)\n"
          f"  shapes -> {csv_path}\n  drift  -> {drift_path}\n  plot   -> {png}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    pr = sub.add_parser("run", help="launch the parallel scan locally; adds a batch and pools")
    pr.add_argument("--config", default="config/scan_gain.json")
    pr.add_argument("--out", default=str(DEFAULT_OUT))
    pr.add_argument("--label", default=None,
                    help="batch subdir name (default: batch_<timestamp>); re-runs never overwrite")
    pr.add_argument("--jobs", type=int, default=None, help="max concurrent processes")
    pr.set_defaults(func=run)

    pe = sub.add_parser("emit", help="write per-job configs + manifest for a scheduler")
    pe.add_argument("--config", default="config/scan_gain.json")
    pe.add_argument("--out", required=True, help="scan root directory")
    pe.add_argument("--label", default=None,
                    help="batch subdir name (default: batch_<timestamp>)")
    pe.add_argument("--xchunks", type=int, default=1,
                    help="split the x-list into K chunks per depth (more array tasks)")
    pe.set_defaults(func=emit)

    for name in ("accumulate", "merge"):  # merge kept as an alias
        pa = sub.add_parser(name, help="pool per-event gains of all batches under --out "
                                       "into gain_map.csv (incremental, no waste)")
        pa.add_argument("--out", required=True)
        pa.set_defaults(func=accumulate)

    pp = sub.add_parser("plot", help="render the 2D gain map from the pooled CSV")
    pp.add_argument("--csv", default=str(DEFAULT_OUT / "gain_map.csv"))
    pp.add_argument("--config", default="config/scan_gain.json")
    pp.add_argument("--rms", action="store_true",
                    help="also report per-point gain mean/rms/sem")
    pp.set_defaults(func=plot)

    pw = sub.add_parser("waveform",
                        help="mean peak-normalized anode current shape vs position "
                             "(ion-drift run; --tcut selectable)")
    pw.add_argument("--out", required=True, help="scan root (pools all batches under it)")
    pw.add_argument("--config", default="config/scan_waveform.json",
                    help="config used for the run (for time_step_ns + geometry)")
    pw.add_argument("--tcut", type=float, default=100.0,
                    help="analysis cutoff time [ns] (default 100; <= captured time_window_ns)")
    pw.set_defaults(func=waveform)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
