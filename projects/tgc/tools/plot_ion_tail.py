#!/usr/bin/env python3
"""
Overlay the analytic ion-tail model i(t) = i0 / (1 + t/t0) on the simulated
anode waveform, to show that the flat shoulder after the electron spike is the
expected ion-induced current (plateau for t < t0, then 1/t).

Physics
-------
Near a wire the field is radial, E ∝ 1/r.  The *induced* current scales as 1/r²
(drift velocity v = μE ∝ 1/r, times the wire Ramo weighting field E_w ∝ 1/r).
Solving the motion in the 1/r field gives r²(t) = r0²(1 + t/t0), hence

    i(t) = i0 / (1 + t/t0),     t0 = r0² / (2 μ k)

so the current is flat while the ion sits near its birth radius r0 (t < t0),
then decays as 1/t.  This is the differential form of Q(t) ∝ ln(1 + t/t0).

Fit
---
scipy is not assumed.  The model linearises as

    1/|i| = (1/i0) + (t - tc)/(i0 t0)

so a degree-1 numpy.polyfit of 1/|i| vs (t - tc) yields i0 and t0 directly.

Usage
-----
    python3 tools/plot_ion_tail.py [--root FILE] [--dist DIR] [--event N]

All arguments are optional; defaults pick the most recent results file.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import uproot

SCRIPT_DIR = Path(__file__).parent.resolve()
TGC_DIR = (SCRIPT_DIR / "..").resolve()


def find_latest_root() -> Path | None:
    """Return the most recently modified results/**/tgc_sim.root, if any."""
    candidates = list((TGC_DIR / "results").glob("**/tgc_sim.root"))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=None,
                    help="tgc_sim.root file (default: most recent under results/)")
    ap.add_argument("--dist", default=None,
                    help="distance sub-directory, e.g. dist_0p7mm_x0mm "
                         "(default: first dist_* found)")
    ap.add_argument("--event", type=int, default=0,
                    help="event index within the t_signals tree (default 0)")
    args = ap.parse_args()

    root_path = args.root or find_latest_root()
    if root_path is None or not Path(root_path).exists():
        print(f"No ROOT file found (looked for {root_path}).", file=sys.stderr)
        return 1

    with uproot.open(root_path) as f:
        dist_dirs = sorted({k.split("/")[0] for k in f.keys(cycle=False)
                            if k.split("/")[0].startswith("dist_")})
        if not dist_dirs:
            print("No dist_* directories in the ROOT file.", file=sys.stderr)
            return 1
        dist = args.dist or dist_dirs[0]
        if dist not in dist_dirs:
            print(f"Distance dir {dist!r} not found. Available: {dist_dirs}",
                  file=sys.stderr)
            return 1

        anode = f[f"{dist}/t_signals"]["anode"].array(library="np")
        if args.event >= len(anode):
            print(f"Event {args.event} out of range (n={len(anode)}).",
                  file=sys.stderr)
            return 1
        i = np.asarray(anode[args.event], dtype=float)        # fC/ns (scaled)
        t = f[f"{dist}/p_anode_signal"].axis(0).centers()      # ns

    a = np.abs(i)
    tc = t[int(np.argmax(a))]                                  # spike / ion-birth time

    # The i0/(1+t/t0) model only holds in the near-wire 1/r region.  Once ions
    # leave it (toward the uniform bulk) the current falls faster than 1/t, so
    # the fit window is restricted to t < tc + 200 ns.  Log-spaced bin averaging
    # suppresses single-event noise (1/|i| would otherwise amplify noisy
    # small-|i| samples into outliers that wreck the least-squares fit).
    lo, hi = tc + 1.5, tc + 200.0
    edges = np.logspace(np.log10(lo), np.log10(hi), 20)
    tb, ab = [], []
    for e0, e1 in zip(edges[:-1], edges[1:]):
        sel = (t >= e0) & (t < e1) & (a > 0.0)
        if sel.sum() == 0:
            continue
        tb.append(np.sqrt(e0 * e1))          # geometric-mean centre
        ab.append(float(a[sel].mean()))
    tb, ab = np.asarray(tb), np.asarray(ab)
    if len(tb) < 5:
        print("Not enough points in the fit window.", file=sys.stderr)
        return 1
    tpb = tb - tc
    slope, intercept = np.polyfit(tpb, 1.0 / ab, 1)            # 1/|i| = a + b*tp
    i0 = 1.0 / intercept
    t0 = intercept / slope
    print(f"ROOT : {root_path}")
    print(f"dist : {dist}   event : {args.event}")
    print(f"tc   = {tc:.2f} ns  (electron-spike / ion-birth reference)")
    print(f"i0   = {i0:.3f} fC/ns")
    print(f"t0   = {t0:.2f} ns")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Plot range: focus on the early plateau + 1/t decay.  Show the raw samples
    # faintly and the log-bin averages (what was fitted) as solid markers.
    pmask = (t >= 10.0) & (t <= 3000.0) & (a > 0.0)
    tt = t[pmask]
    mmask = tt > tc                       # model only valid for t > tc
    model = i0 / (1.0 + (tt[mmask] - tc) / t0)

    fig, (ax, axl) = plt.subplots(1, 2, figsize=(13, 5))

    # ── Left panel: log-log (power-law check) ────────────────────────────────
    ax.loglog(tt, a[pmask], ".", ms=1.5, color="#85B7EB", alpha=0.4,
              label="raw samples (single event)")
    ax.loglog(tb, ab, "o", ms=5, color="#185FA5",
              label="log-bin average (fitted)")
    ax.loglog(tt[mmask], model, "-", lw=2, color="#D85A30",
              label=r"$i_0/(1+t/t_0)$ model")
    ax.axvline(tc + t0, ls="--", lw=1, color="grey")
    ax.text(tc + t0, a[pmask].min() * 1.5, r"  $t_0$", color="grey", va="bottom")

    # Regime annotations.
    ax.annotate("plateau (t < t₀)", xy=(tc + 0.4 * t0, i0),
                xytext=(0, 16), textcoords="offset points",
                ha="center", color="#993C1D", fontsize=10)
    ax.annotate("1/t  (t > t₀)", xy=(800.0, i0 / (1.0 + (800.0 - tc) / t0)),
                xytext=(12, 10), textcoords="offset points",
                color="#993C1D", fontsize=10)

    ax.set_xlabel("time  t  [ns]")
    ax.set_ylabel("anode induced current  |i|  [fC/ns]")
    ax.set_title("log-log")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    # ── Right panel: linear (plateau shape, scope-trace view) ────────────────
    axl.plot(tt, a[pmask], ".", ms=2, color="#85B7EB", alpha=0.4,
             label="raw samples (single event)")
    axl.plot(tb, ab, "o", ms=5, color="#185FA5",
             label="log-bin average (fitted)")
    axl.plot(tt[mmask], model, "-", lw=2, color="#D85A30",
             label=r"$i_0/(1+t/t_0)$ model")
    axl.axvline(tc + t0, ls="--", lw=1, color="grey")
    axl.text(tc + t0, 0.04 * i0, r"  $t_0$", color="grey", va="bottom")
    axl.set_xlim(tc - 10.0, 500.0)
    axl.set_ylim(0.0, 1.4 * i0)
    axl.annotate("e⁻ spike off-scale", xy=(tc, 1.4 * i0),
                 xytext=(8, -14), textcoords="offset points",
                 color="grey", fontsize=9)
    axl.set_xlabel("time  t  [ns]")
    axl.set_ylabel("anode induced current  |i|  [fC/ns]")
    axl.set_title("linear")
    axl.grid(True, alpha=0.3)
    axl.legend()

    r0_um = np.sqrt(2 * 1.2 * 390 * t0 * 1e-9) * 1e4 if t0 > 0 else float("nan")
    fig.suptitle(f"Anode ion tail — i0 = {i0:.2f} fC/ns,  t0 = {t0:.1f} ns "
                 f"(implied r0 ≈ {r0_um:.0f} µm)")

    out = Path(root_path).parent / "ion_tail_fit.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"saved: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
