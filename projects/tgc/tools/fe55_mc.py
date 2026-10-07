#!/usr/bin/env python3
"""Light 55Fe Monte Carlo — synthesize realistic TGC pulses and the 55Fe charge
spectrum from the *measured* single-electron response, with no further Garfield runs.

Physics factorization (established by the position scan):
  * the peak-aligned, peak-normalized single-electron anode shape  S(t')  is
    POSITION-INDEPENDENT        -> one canonical shape (pooled over positions),
  * the drift delay             tau(x,y)   and its jitter rms_tau(x,y),  and
  * the avalanche-gain distribution  Gain(x,y)   (per-position Polya, incl. the
    zero-gain attachment fraction)
  are POSITION-DEPENDENT and already mapped.

A real event is a superposition of single-electron responses, so a full pulse is
    W(t) = sum_i  g_i * S(t - tau_i),      g_i ~ Gain(x_i,y_i),  tau_i ~ N(tau,rms)(x_i,y_i)
built by summing the measured shape over a realistic cloud of primary electrons.
`tgc_sim` instead fires ONE avalanche and scales by nPrimary=round(E/W); this MC is
the unfolding of that approximation (see src/tgc_sim.cc:1139).

55Fe model (Ka + Ar-escape; photoelectron track + Auger cloud):
  5.9 keV photon -> photoabsorb on Ar K (Eb=3.206 keV) -> photoelectron (2.69 keV)
  track + K-shell relaxation: Auger (~88%, deposits Eb locally -> main peak) or Ar
  Ka fluorescence (~12%, 2.96 keV photon; escapes the thin gap -> ESCAPE peak at
  2.94 keV, or reabsorbs -> satellite cluster -> main peak). Fano on each cluster.

Inputs (scan products; nothing re-simulated):
  --shape-dir  waveform scan -> S(t')  (waveform_shapes.csv, pooled)
                              + tau map (waveform_drift.csv: mean/rms drift)
  --gain-dir   gain scan     -> per-position empirical gain (avalanche_size in roots)

Outputs (--out):  fe55_spectrum.{png,csv}, fe55_waveforms.png, fe55_observables.csv
Reuses _parse_tag / _wire_positions_cm from gain_scan.py.  See config/fe55_mc.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from gain_scan import _parse_tag, _wire_positions_cm  # noqa: E402  (shared helpers)

TGC_DIR = (SCRIPT_DIR / "..").resolve()
DEFAULTS = {
    "geometry": {"wire_pitch_cm": 0.18, "wire_diameter_um": 50.0, "gap_cm": 0.14,
                 "n_wires": 10},
    "gas": {"w_value_eV": 26.0, "fano_factor": 0.20,
            "atten_len_5900_cm": 16.0, "atten_len_2960_cm": 2.0},
    "xray": {"lines": [{"name": "Ka", "energy_eV": 5895.0, "weight": 1.0}],
             "ar_k_binding_eV": 3206.0, "ar_k_fluor_yield": 0.12,
             "ar_ka_fluor_eV": 2957.0, "photoelectron_angular": "dipole"},
    "track": {"range_um_at_1keV": 30.0, "range_exponent": 1.7, "straggle_frac": 0.3},
    "diffusion": {"enable_transverse": True, "sigma_um_per_sqrt_cm": 120.0},
    "scan_points": [
        {"name": "on-wire shallow", "x_cm": 0.09, "depth_mm": 0.3},
        {"name": "on-wire deep", "x_cm": 0.09, "depth_mm": 1.3},
        {"name": "mid-gap shallow", "x_cm": 0.18, "depth_mm": 0.3},
        {"name": "mid-gap deep", "x_cm": 0.18, "depth_mm": 1.3},
    ],
    "run": {"n_photons": 20000, "n_photons_scan": 3000, "tcut_ns": 100.0,
            "out_time_step_ns": 0.5, "seed": 1},
}


# ───────────────────────────── input maps ──────────────────────────────────

def _merge(base: dict, over: dict) -> dict:
    out = json.loads(json.dumps(base))
    for k, v in (over or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def load_shape(shape_dir: Path, tcut: float):
    """Canonical single-electron shape S(t'): pool `mean_norm_current` over all
    positions (position-independent), restrict to [-5, tcut], renormalize peak to 1."""
    df = pd.read_csv(shape_dir / "waveform_shapes.csv")
    g = (df.groupby("time_since_peak_ns")["mean_norm_current"]
           .mean().reset_index().sort_values("time_since_peak_ns"))
    ts = g["time_since_peak_ns"].to_numpy(float)
    S = np.nan_to_num(g["mean_norm_current"].to_numpy(float), nan=0.0)
    keep = (ts >= -5.0 - 1e-9) & (ts <= tcut + 1e-9)
    ts, S = ts[keep], S[keep]
    pk = S.max()
    if pk > 0:
        S = S / pk
    ipk = int(np.argmin(np.abs(ts)))        # index of t'=0 (the spike)
    return ts, S, ipk


def _grid_from_drift(shape_dir: Path, wires: np.ndarray):
    """Build the regular (depth_mm, d_from_wire_mm) grid of mean/rms drift [ns].
    d = distance to the nearest wire (the physical variable; folds the full pitch)."""
    df = pd.read_csv(shape_dir / "waveform_drift.csv")
    df["d_mm"] = df["x_position_cm"].apply(
        lambda x: float(np.min(np.abs(x - wires))) * 10.0).round(4)
    df["depth_mm"] = df["source_distance_mm"].round(4)
    ys = np.array(sorted(df["depth_mm"].unique()))
    ds = np.array(sorted(df["d_mm"].unique()))
    tau = np.full((len(ys), len(ds)), np.nan)
    rms = np.full((len(ys), len(ds)), np.nan)
    for _, r in df.iterrows():
        iy = int(np.argmin(np.abs(ys - r["depth_mm"])))
        idd = int(np.argmin(np.abs(ds - r["d_mm"])))
        tau[iy, idd], rms[iy, idd] = r["mean_drift_ns"], r["rms_drift_ns"]
    return ys, ds, _fill_nan(tau), _fill_nan(rms)


def _fill_nan(Z: np.ndarray) -> np.ndarray:
    """Nearest-neighbour fill of any empty grid cell (missing scan point)."""
    if not np.isnan(Z).any():
        return Z
    iy, idx = np.where(~np.isnan(Z))
    jy, jdx = np.where(np.isnan(Z))
    for y, x in zip(jy, jdx):
        k = int(np.argmin((iy - y) ** 2 + (idx - x) ** 2))
        Z[y, x] = Z[iy[k], idx[k]]
    return Z


def load_gain_cells(gain_dir: Path, wires: np.ndarray, ys: np.ndarray, ds: np.ndarray):
    """Per grid-cell empirical gain samples (avalanche_size) + per-cell mean/rms,
    keyed onto the same (depth_mm, d_mm) grid as the drift map."""
    import uproot
    pools = defaultdict(list)
    for root in sorted(Path(gain_dir).glob("**/tgc_sim.root")):
        try:
            with uproot.open(root) as f:
                for name in f.keys(cycle=False):
                    if name.rsplit("/", 1)[-1] != "t_signals":
                        continue
                    dist, x = _parse_tag(name.rsplit("/", 1)[0])
                    if dist is None or x is None:
                        continue
                    arr = f[name]["avalanche_size"].array(library="np").astype(float)
                    if len(arr):
                        pools[(round(dist, 4), round(x, 4))].append(arr)
        except Exception as e:  # noqa: BLE001 — skip a partial/corrupt file
            print(f"warning: skipping {root}: {e}", file=sys.stderr)
    if not pools:
        raise SystemExit(f"no t_signals/avalanche_size under {gain_dir}")

    cells = np.empty((len(ys), len(ds)), dtype=object)
    gmean = np.full((len(ys), len(ds)), np.nan)
    grms = np.full((len(ys), len(ds)), np.nan)
    allg = []
    for (dist, x), arrs in pools.items():
        d_mm = round(float(np.min(np.abs(x - wires))) * 10.0, 4)
        iy = int(np.argmin(np.abs(ys - dist)))
        idd = int(np.argmin(np.abs(ds - d_mm)))
        g = np.concatenate(arrs)
        cells[iy, idd] = g
        gmean[iy, idd], grms[iy, idd] = g.mean(), g.std()
        allg.append(g)
    fallback = np.concatenate(allg)
    for iy in range(len(ys)):           # fill any empty cell with the global pool
        for idd in range(len(ds)):
            if cells[iy, idd] is None:
                cells[iy, idd] = fallback
                gmean[iy, idd], grms[iy, idd] = fallback.mean(), fallback.std()
    return cells, _fill_nan(gmean), _fill_nan(grms)


# ───────────────────────────── interpolation ───────────────────────────────

def _bilin(ys, ds, Z, yq, dq):
    """Bilinear interpolation of Z on the regular (ys, ds) grid, edge-clamped."""
    yq = np.clip(yq, ys[0], ys[-1])
    dq = np.clip(dq, ds[0], ds[-1])
    iy = np.clip(np.searchsorted(ys, yq) - 1, 0, len(ys) - 2)
    idd = np.clip(np.searchsorted(ds, dq) - 1, 0, len(ds) - 2)
    y0, y1 = ys[iy], ys[iy + 1]
    d0, d1 = ds[idd], ds[idd + 1]
    ty = np.where(y1 > y0, (yq - y0) / (y1 - y0), 0.0)
    td = np.where(d1 > d0, (dq - d0) / (d1 - d0), 0.0)
    return (Z[iy, idd] * (1 - ty) * (1 - td) + Z[iy, idd + 1] * (1 - ty) * td
            + Z[iy + 1, idd] * ty * (1 - td) + Z[iy + 1, idd + 1] * ty * td)


# ─────────────────────────── primary-pattern generation ────────────────────

def _perp(u):
    a = np.array([1.0, 0.0, 0.0]) if abs(u[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(u, a); e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    return e1, e2


def _iso(rng):
    c = rng.uniform(-1, 1); s = np.sqrt(1 - c * c); phi = rng.uniform(0, 2 * np.pi)
    return np.array([s * np.cos(phi), s * np.sin(phi), c])


def _photoelectron_dir(rng, p, beta, mode):
    if mode != "dipole":
        return _iso(rng)
    M = 1.0 + beta
    while True:                                   # reject onto (1-c^2)(1+beta c)
        c = rng.uniform(-1, 1)
        if rng.uniform(0, 1) < (1 - c * c) * (1 + beta * c) / M:
            break
    s = np.sqrt(1 - c * c); phi = rng.uniform(0, 2 * np.pi)
    e1, e2 = _perp(p)
    return c * p + s * (np.cos(phi) * e1 + np.sin(phi) * e2)


def _range_cm(E_keV, trk):
    return trk["range_um_at_1keV"] * max(E_keV, 1e-3) ** trk["range_exponent"] * 1e-4


def _fano_n(mean_n, F, rng):
    n = int(round(rng.normal(mean_n, np.sqrt(max(F * mean_n, 0.0)))))
    return max(n, 0)


def _track(origin, direction, E_eV, phys, rng):
    """Electrons of one energetic deposit: n=E/W (Fano) spread along a range-length
    segment with transverse straggling. Returns (n,3) cm positions."""
    n = _fano_n(E_eV / phys["w"], phys["F"], rng)
    if n == 0:
        return np.empty((0, 3))
    R = _range_cm(E_eV / 1000.0, phys["trk"])
    u = direction / np.linalg.norm(direction)
    e1, e2 = _perp(u)
    s = R * rng.random(n)
    st = phys["trk"]["straggle_frac"] * R
    pos = (origin + s[:, None] * u
           + st * (rng.standard_normal(n)[:, None] * e1
                   + rng.standard_normal(n)[:, None] * e2))
    return pos


def generate_electrons(conv, phys, rng):
    """One 55Fe photon at `conv`=(x,y,z) cm -> (positions (N,3) cm, truth label)."""
    xr = phys["xray"]
    E_photon = conv[3]                                  # chosen line energy [eV]
    E_pe = E_photon - xr["ar_k_binding_eV"]
    p = np.array([0.0, -1.0, 0.0])                      # photon travels into the gas
    beta = np.sqrt(1 - 1 / (1 + E_pe / 511000.0) ** 2)
    origin = np.array(conv[:3])
    emissions = [_track(origin, _photoelectron_dir(rng, p, beta, xr["photoelectron_angular"]),
                        E_pe, phys, rng)]
    truth = "main"
    if rng.random() < xr["ar_k_fluor_yield"]:          # K fluorescence
        E_L = xr["ar_k_binding_eV"] - xr["ar_ka_fluor_eV"]
        emissions.append(_track(origin, _iso(rng), E_L, phys, rng))   # L-shell Auger
        end = origin + rng.exponential(phys["lam2960"]) * _iso(rng)
        if abs(end[1]) >= phys["gap_cm"]:              # fluor photon leaves the gap
            truth = "escape"
        else:
            emissions.append(_track(end, _iso(rng), xr["ar_ka_fluor_eV"], phys, rng))
    else:                                              # Auger (full energy local)
        emissions.append(_track(origin, _iso(rng), xr["ar_k_binding_eV"], phys, rng))
    pos = np.vstack([e for e in emissions if len(e)]) if any(len(e) for e in emissions) \
        else np.empty((0, 3))
    return pos, truth


# ─────────────────────────────── transport + sum ───────────────────────────

def _sample_gain(iy, idd, cells, gmean, grms, mode, rng):
    """One gain draw per electron. Empirical: resample the cell's avalanche_size
    array (keeps the true Polya + zero-gain fraction). Polya: Gamma(mean,rms)."""
    n = len(iy)
    g = np.empty(n)
    cell = iy * cells.shape[1] + idd
    for c in np.unique(cell):
        m = cell == c
        ry, rd = divmod(int(c), cells.shape[1])
        k = int(m.sum())
        if mode == "polya":
            mu, sig = gmean[ry, rd], grms[ry, rd]
            theta = max(mu * mu / sig / sig - 1.0, 0.0) if sig > 0 else 0.0
            g[m] = rng.gamma(1.0 + theta, mu / (1.0 + theta)) if mu > 0 else 0.0
        else:
            g[m] = rng.choice(cells[ry, rd], size=k, replace=True)
    return g


def _risetime_ns(T, W):
    pk = int(np.argmax(W)); A = W[pk]
    if A <= 0 or pk == 0:
        return float("nan")

    def cross(frac):
        thr = frac * A
        seg = W[:pk + 1]
        below = np.where(seg < thr)[0]
        if len(below) == 0:
            return T[0]
        i = below[-1]
        if i + 1 > pk:
            return T[pk]
        w0, w1 = seg[i], seg[i + 1]
        f = 0.0 if w1 == w0 else (thr - w0) / (w1 - w0)
        return T[i] + f * (T[i + 1] - T[i])
    return cross(0.9) - cross(0.1)


def _fwhm_ns(t, y):
    """Full width at half maximum [ns] of a pulse y(t): interval between the two
    half-max (0.5*peak) crossings bracketing the peak, linearly interpolated. NaN if
    the peak is non-positive; clamps to the window edge if a side never drops below half."""
    y = np.asarray(y, float); t = np.asarray(t, float)
    if len(y) == 0 or not np.isfinite(np.nanmax(y)):
        return float("nan")
    pk = int(np.nanargmax(y)); A = y[pk]
    if A <= 0:
        return float("nan")
    half = 0.5 * A

    def edge(step):
        i = pk
        while 0 <= i + step < len(y) and y[i + step] >= half:
            i += step
        nxt = i + step
        if not (0 <= nxt < len(y)) or y[nxt] >= half:
            return t[i]                               # never crossed -> clamp to edge
        f = (y[i] - half) / (y[i] - y[nxt])           # y[i] >= half > y[nxt]
        return t[i] + f * (t[nxt] - t[i])
    return edge(1) - edge(-1)


def simulate(convs, phys, maps, S, ipk, T, dt, gain_mode, rng, keep_samples=0):
    """Transport+sum for a list of conversion points (each (x,y,z,E_eV)).
    Returns dict with sum_W, n, sample waveforms, and per-photon observable arrays."""
    ys, ds, tau, rms = maps["ys"], maps["ds"], maps["tau"], maps["rms"]
    cells, gmean, grms = maps["cells"], maps["gmean"], maps["grms"]
    wires, gap = maps["wires"], phys["gap_cm"]
    nT, nS = len(T), len(S)
    sum_W = np.zeros(nT)
    # Peak-aligned, peak-normalized accumulator (align each pulse on its own peak,
    # THEN average -> the intrinsic mean shape, free of the drift-time smearing that
    # an absolute-time average carries). Mirrors gain_scan.py waveform().
    n_pre = min(int(round(10.0 / dt)), nT - 1)          # bins shown before the peak
    n_post = nS - 1 - ipk                               # tail bins (= the S tail, ~tcut)
    alen = n_pre + n_post + 1
    t_rel = (np.arange(alen) - n_pre) * dt
    al_sum, al_wt = np.zeros(alen), np.zeros(alen)
    samples, Q, amp, tpeak, rise, truth = [], [], [], [], [], []
    for conv in convs:
        pos, tr = generate_electrons(conv, phys, rng)
        if len(pos) == 0:
            continue
        y = np.clip(pos[:, 1], -gap + 1e-4, gap - 1e-4)
        depth = np.clip(np.abs(y) * 10.0, ys[0], ys[-1])          # mm
        xx = pos[:, 0]
        if phys["diff"]["enable_transverse"]:                      # transverse diffusion
            sig = phys["diff"]["sigma_um_per_sqrt_cm"] * 1e-4 * np.sqrt(np.abs(y))
            xx = xx + sig * rng.standard_normal(len(xx))
        nearest = np.abs(xx[:, None] - wires).argmin(1)
        d_mm = np.clip(np.abs(xx - wires[nearest]) * 10.0, ds[0], ds[-1])
        tmean = _bilin(ys, ds, tau, depth, d_mm)
        trms = _bilin(ys, ds, rms, depth, d_mm)
        taus = np.clip(tmean + trms * rng.standard_normal(len(d_mm)), 0.0, None)
        iy = np.abs(depth[:, None] - ys).argmin(1)
        idd = np.abs(d_mm[:, None] - ds).argmin(1)
        g = _sample_gain(iy, idd, cells, gmean, grms, gain_mode, rng)

        bins = np.floor((taus - T[0]) / dt).astype(int)
        ok = (bins >= 0) & (bins < nT)
        A = np.bincount(bins[ok], weights=g[ok], minlength=nT)
        W = np.convolve(A, S)[ipk:ipk + nT]
        sum_W += W
        pk = int(np.argmax(W)); wp = W[pk]
        if wp > 0:                                      # align on this pulse's peak
            lo, hi = pk - n_pre, pk + n_post + 1
            sa, sb = max(0, lo), min(nT, hi)
            da = sa - lo
            al_sum[da:da + (sb - sa)] += W[sa:sb] / wp
            al_wt[da:da + (sb - sa)] += 1.0
        Q.append(float(g.sum())); amp.append(float(W.max()))
        tpeak.append(float(T[pk])); rise.append(_risetime_ns(T, W))
        truth.append(tr)
        if len(samples) < keep_samples:
            samples.append(W.copy())
    aligned_mean = np.divide(al_sum, al_wt, out=np.full(alen, np.nan), where=al_wt > 0)
    return {"sum_W": sum_W, "n": len(Q), "samples": samples,
            "aligned_mean": aligned_mean, "t_rel": t_rel,
            "Q": np.array(Q), "amp": np.array(amp), "tpeak": np.array(tpeak),
            "rise": np.array(rise), "truth": np.array(truth)}


# ─────────────────────────────── sampling conv points ──────────────────────

def _pick_line(xr, rng):
    lines = xr["lines"]
    w = np.array([ln["weight"] for ln in lines], float); w /= w.sum()
    return float(lines[int(rng.choice(len(lines), p=w))]["energy_eV"])


def realistic_convs(n, phys, rng):
    """Random conversions: depth from exponential X-ray absorption through the gap
    (truncated to the 2*gap crossing; lam>>gap -> ~uniform), x uniform over a pitch,
    z arbitrary. Vectorized via the inverse CDF of the gap-truncated exponential."""
    gap, pitch, lam = phys["gap_cm"], phys["pitch_cm"], phys["lam5900"]
    wire0 = phys["wires"][len(phys["wires"]) // 2]
    ell = -lam * np.log1p(-rng.random(n) * (1.0 - np.exp(-2.0 * gap / lam)))
    y = gap - ell                                     # in (-gap, gap)
    x = wire0 + rng.uniform(-pitch / 2, pitch / 2, n)
    lines = phys["xray"]["lines"]
    w = np.array([ln["weight"] for ln in lines], float); w /= w.sum()
    E = np.array([ln["energy_eV"] for ln in lines])[rng.choice(len(lines), size=n, p=w)]
    return list(zip(x, y, np.zeros(n), E))


# ─────────────────────────────────── plotting ──────────────────────────────

def _plot_spectrum(res, out, phys):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Q = res["Q"] / 1e3                                  # 10^3 gain-electrons
    tr = res["truth"]
    main, esc = Q[tr == "main"], Q[tr == "escape"]
    cen_main = float(np.mean(main)) if len(main) else float("nan")
    cen_esc = float(np.mean(esc)) if len(esc) else float("nan")
    res_main = float(np.std(main) / np.mean(main)) if len(main) else float("nan")
    keV_per_Q = 5.895 / cen_main if cen_main > 0 else float("nan")

    bins = np.linspace(0, np.percentile(Q, 99.5) * 1.1, 120)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(Q, bins=bins, color="0.6", label=f"all (N={len(Q)})")
    ax.hist(main, bins=bins, histtype="step", color="C0", lw=1.6,
            label=f"main 5.9 keV (centroid={cen_main:.0f}k, "
                  rf"$\sigma/E$={res_main*100:.1f}%)")
    if len(esc):
        ax.hist(esc, bins=bins, histtype="step", color="C3", lw=1.6,
                label=f"Ar escape (centroid={cen_esc:.0f}k, "
                      f"esc/main={cen_esc/cen_main:.2f})")
    ax.axvline(cen_main, color="C0", ls=":", lw=1)
    if len(esc):
        ax.axvline(cen_esc, color="C3", ls=":", lw=1)
    ax.set_xlabel(r"collected charge  $Q=\sum_i g_i$  [$10^3$ gain-electrons]")
    ax.set_ylabel("photons / bin")
    ax.set_title(f"55Fe charge spectrum  (Ka+escape; {len(Q)} photons)")
    sec = ax.secondary_xaxis("top", functions=(lambda q: q * keV_per_Q,
                                               lambda e: e / keV_per_Q))
    sec.set_xlabel("energy [keV]  (self-calibrated to the 5.9 keV peak)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "fe55_spectrum.png", dpi=130)
    pd.DataFrame({"Q_centre_k": 0.5 * (bins[1:] + bins[:-1]),
                  "counts_all": np.histogram(Q, bins)[0],
                  "counts_main": np.histogram(main, bins)[0],
                  "counts_escape": np.histogram(esc, bins)[0]}
                 ).to_csv(out / "fe55_spectrum.csv", index=False)
    return {"centroid_main_k": cen_main, "centroid_escape_k": cen_esc,
            "sigmaE_over_E_main": res_main,
            "escape_over_main": (cen_esc / cen_main) if cen_main > 0 else float("nan"),
            "escape_fraction": float(np.mean(tr == "escape"))}


def _plot_waveforms(scan_res, incl_res, T, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(13, 5))
    for name, r in scan_res.items():
        if r["n"] == 0:
            continue
        axa.plot(r["t_rel"], r["aligned_mean"], lw=1.5,
                 label=f"{name}  (peak-t={np.nanmean(r['tpeak']):.1f}ns, "
                       f"rise={np.nanmean(r['rise']):.1f}ns)")
    if incl_res and incl_res["n"]:                       # each sample aligned to its own peak
        for W in incl_res["samples"][:8]:
            pk = int(np.argmax(W))
            if W[pk] > 0:
                axb.plot(T - T[pk], W / W[pk], color="C0", lw=0.6, alpha=0.35)
        axb.plot(incl_res["t_rel"], incl_res["aligned_mean"], color="k", lw=1.8,
                 label=f"aligned mean (N={incl_res['n']})")
    for ax in (axa, axb):
        ax.set_xlim(-5, 80); ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="0.7", lw=0.5, ls=":")
        ax.set_xlabel("time since peak [ns]"); ax.set_ylabel("pulse / peak")
    axa.set_title("Mean 55Fe pulse vs conversion position (peak-aligned)")
    axb.set_title("Realistic-exposure pulses (peak-aligned, normalized)")
    for ax in (axa, axb):
        if ax.get_legend_handles_labels()[1]:
            ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out / "fe55_waveforms.png", dpi=130)


def _plot_grid(meta, scan_res, out, wires):
    """2D maps of the pulse observables across the full (depth, distance-from-wire) grid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for name, m in meta.items():
        r = scan_res[name]
        if not r["n"]:
            continue
        rows.append({"d_wire_mm": m["d_wire_mm"], "depth_mm": m["depth_mm"],
                     "peak_t": float(np.nanmean(r["tpeak"])), "rise": float(np.nanmean(r["rise"])),
                     "amp": float(r["amp"].mean()), "Q": float(r["Q"].mean()) / 1e3})
    df = pd.DataFrame(rows)
    panels = [("peak_t", "Peak time [ns]  (drift timing)", "viridis"),
              ("rise", "Rise time 10-90% [ns]  (pulse shape)", "magma"),
              ("amp", "Mean peak amplitude [arb.]", "cividis"),
              ("Q", r"Mean charge  $\sum_i g_i$  [$10^3$]", "plasma")]
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))
    for ax, (col, title, cmap) in zip(axes.flat, panels):
        piv = df.pivot_table(index="depth_mm", columns="d_wire_mm", values=col)
        im = ax.pcolormesh(piv.columns.to_numpy(), piv.index.to_numpy(), piv.to_numpy(),
                           shading="nearest", cmap=cmap)
        ax.set_xlabel("distance from wire [mm]  (0 = on wire)")
        ax.set_ylabel("conversion depth [mm]")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
    fig.suptitle("55Fe pulse observables across the full grid", y=1.0)
    fig.tight_layout()
    fig.savefig(out / "fe55_grid.png", dpi=130)


def _plot_fwhm(meta, scan_res, out):
    """Standalone 2D map of the pulse FWHM across the (depth, distance-from-wire) grid."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for name, m in meta.items():
        r = scan_res[name]
        if not r["n"]:
            continue
        rows.append({"d_wire_mm": m["d_wire_mm"], "depth_mm": m["depth_mm"],
                     "fwhm": _fwhm_ns(r["t_rel"], r["aligned_mean"])})
    piv = pd.DataFrame(rows).pivot_table(index="depth_mm", columns="d_wire_mm", values="fwhm")
    fig, ax = plt.subplots(figsize=(7.5, 6))
    im = ax.pcolormesh(piv.columns.to_numpy(), piv.index.to_numpy(), piv.to_numpy(),
                       shading="nearest", cmap="magma")
    ax.set_xlabel("distance from wire [mm]  (0 = on wire)")
    ax.set_ylabel("conversion depth [mm]")
    ax.set_title("55Fe pulse FWHM vs position")
    fig.colorbar(im, ax=ax, label="FWHM [ns]")
    fig.tight_layout()
    fig.savefig(out / "fe55_fwhm.png", dpi=130)


def _plot_sweep(meta, scan_res, out):
    """Aligned pulse shapes swept along one grid axis at a time: transverse (fixed
    depth) shows the broadening toward mid-gap; depth (fixed transverse) shows the
    shape is ~depth-independent. Uses the per-point aligned_mean from grid mode."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    ys_u = sorted({m["depth_mm"] for m in meta.values()})
    ds_u = sorted({m["d_wire_mm"] for m in meta.values()})
    depth_fixed = ys_u[len(ys_u) // 2]                   # transverse sweep at this depth
    d_fixed = ds_u[0]                                    # depth sweep on the wire

    def _find(depth, d):
        for name, m in meta.items():
            if abs(m["depth_mm"] - depth) < 1e-6 and abs(m["d_wire_mm"] - d) < 1e-6:
                return scan_res[name]
        return None

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (vals, fixed, cmap, norm_lo, norm_hi, cbar, title, finder) in [
        (axL, (ds_u, depth_fixed, "viridis", ds_u[0], ds_u[-1], "distance from wire [mm]",
               f"Transverse sweep at depth = {depth_fixed:g} mm", lambda v: _find(depth_fixed, v))),
        (axR, (ys_u, d_fixed, "plasma", ys_u[0], ys_u[-1], "conversion depth [mm]",
               f"Depth sweep at distance from wire = {d_fixed:g} mm", lambda v: _find(v, d_fixed))),
    ]:
        cm = plt.get_cmap(cmap); norm = Normalize(norm_lo, norm_hi)
        for v in vals:
            r = finder(v)
            if r and r["n"]:
                ax.plot(r["t_rel"], r["aligned_mean"], color=cm(norm(v)), lw=1.3)
        sm = plt.cm.ScalarMappable(norm=norm, cmap=cm); sm.set_array([])
        fig.colorbar(sm, ax=ax, label=cbar)
        ax.set_xlim(-5, 80); ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="0.6", lw=0.5, ls=":")
        ax.set_xlabel("time since peak [ns]"); ax.set_ylabel("mean pulse / peak")
        ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out / "fe55_sweep.png", dpi=130)


# ──────────────────────────────────── main ─────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description="Light 55Fe Monte Carlo from the single-electron response")
    ap.add_argument("--config", default=str(TGC_DIR / "config" / "fe55_mc.json"))
    ap.add_argument("--shape-dir", default=str(TGC_DIR / "results" / "scan_waveform_halfpitch"))
    ap.add_argument("--gain-dir", default=str(TGC_DIR / "results" / "gain_scan_halfpitch"))
    ap.add_argument("--out", default=str(TGC_DIR / "results" / "fe55_mc"))
    ap.add_argument("--n-photons", type=int, default=None, help="realistic-exposure photons")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--tcut", type=float, default=None, help="ns of single-electron tail used")
    ap.add_argument("--gain", choices=["empirical", "polya"], default="empirical")
    ap.add_argument("--mode", choices=["both", "scan", "realistic"], default="both")
    ap.add_argument("--scan-grid", action="store_true",
                    help="scan every (depth, x) point of the measured map grid (-> observable "
                         "maps fe55_grid.png) instead of the config's scan_points")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    cfg = _merge(DEFAULTS, json.loads(cfg_path.read_text()) if cfg_path.exists() else {})
    run = cfg["run"]
    seed = args.seed if args.seed is not None else run["seed"]
    tcut = args.tcut if args.tcut is not None else run["tcut_ns"]
    n_real = args.n_photons if args.n_photons is not None else run["n_photons"]
    dt = run["out_time_step_ns"]
    rng = np.random.default_rng(seed)

    shape_dir, gain_dir, out = Path(args.shape_dir), Path(args.gain_dir), Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    wires = _wire_positions_cm(cfg["geometry"])

    print(f"loading maps:  shape={shape_dir.name}  gain={gain_dir.name}")
    ts, S, ipk = load_shape(shape_dir, tcut)
    ys, ds, tau, rms = _grid_from_drift(shape_dir, wires)
    cells, gmean, grms = load_gain_cells(gain_dir, wires, ys, ds)
    print(f"  shape: {len(S)} bins, t'=[{ts[0]:g},{ts[-1]:g}]ns, peak@idx{ipk}; "
          f"grid: {len(ys)} depths x {len(ds)} d; <g>={np.nanmean(gmean):.0f}")

    phys = {
        "w": cfg["gas"]["w_value_eV"], "F": cfg["gas"]["fano_factor"],
        "gap_cm": cfg["geometry"]["gap_cm"], "pitch_cm": cfg["geometry"]["wire_pitch_cm"],
        "lam5900": cfg["gas"]["atten_len_5900_cm"], "lam2960": cfg["gas"]["atten_len_2960_cm"],
        "xray": cfg["xray"], "trk": cfg["track"], "diff": cfg["diffusion"], "wires": wires,
    }
    maps = {"ys": ys, "ds": ds, "tau": tau, "rms": rms,
            "cells": cells, "gmean": gmean, "grms": grms, "wires": wires}
    T = np.arange(-10.0, float(np.nanmax(tau)) + tcut + 15.0 + dt, dt)

    scan_res, scan_meta, incl_res, spec = {}, {}, None, None
    if args.mode in ("both", "scan"):
        if args.scan_grid:                               # every measured map grid point
            wire0 = wires[len(wires) // 2]
            specs = [(f"y{yv:g}_d{dv:g}", float(wire0 + dv / 10.0), float(yv), float(dv))
                     for yv in ys for dv in ds]
            print(f"  scan grid: {len(specs)} points x {run['n_photons_scan']} photons")
        else:
            specs = [(sp["name"], float(sp["x_cm"]), float(sp["depth_mm"]),
                      float(np.min(np.abs(sp["x_cm"] - wires)) * 10.0))
                     for sp in cfg["scan_points"]]
        for name, x_cm, depth_mm, d_mm in specs:
            E = _pick_line(cfg["xray"], rng)
            convs = [(x_cm, -depth_mm / 10.0, 0.0, E) for _ in range(run["n_photons_scan"])]
            scan_res[name] = simulate(convs, phys, maps, S, ipk, T, dt,
                                      args.gain, rng, keep_samples=0)
            scan_meta[name] = {"x_cm": x_cm, "depth_mm": depth_mm, "d_wire_mm": d_mm}
            r = scan_res[name]
            if not args.scan_grid:
                print(f"  scan '{name}': N={r['n']}  <Q>={r['Q'].mean():.0f}  "
                      f"peak-t={np.nanmean(r['tpeak']):.1f}ns  rise={np.nanmean(r['rise']):.1f}ns")
    if args.mode in ("both", "realistic"):
        convs = realistic_convs(n_real, phys, rng)
        incl_res = simulate(convs, phys, maps, S, ipk, T, dt, args.gain, rng, keep_samples=8)
        spec = _plot_spectrum(incl_res, out, phys)
        print(f"  realistic: N={incl_res['n']}  main={spec['centroid_main_k']:.0f}k  "
              f"escape={spec['centroid_escape_k']:.0f}k  esc/main={spec['escape_over_main']:.2f}  "
              f"sigmaE/E={spec['sigmaE_over_E_main']*100:.1f}%  "
              f"esc-frac={spec['escape_fraction']*100:.1f}%")

    outputs = ["fe55_spectrum.png/.csv"] if incl_res else []
    if args.scan_grid and scan_res:
        _plot_grid(scan_meta, scan_res, out, wires); outputs.append("fe55_grid.png")
        _plot_fwhm(scan_meta, scan_res, out); outputs.append("fe55_fwhm.png")
        _plot_sweep(scan_meta, scan_res, out); outputs.append("fe55_sweep.png")
        ys_u = sorted({m["depth_mm"] for m in scan_meta.values()})
        ds_u = sorted({m["d_wire_mm"] for m in scan_meta.values()})
        overlay = {}                                     # 4 corners -> readable shape overlay
        for yv, dv, lbl in [(ys_u[0], ds_u[0], "on-wire shallow"), (ys_u[-1], ds_u[0], "on-wire deep"),
                            (ys_u[0], ds_u[-1], "mid-gap shallow"), (ys_u[-1], ds_u[-1], "mid-gap deep")]:
            for name, m in scan_meta.items():
                if abs(m["depth_mm"] - yv) < 1e-6 and abs(m["d_wire_mm"] - dv) < 1e-6:
                    overlay[lbl] = scan_res[name]; break
        _plot_waveforms(overlay, incl_res, T, out); outputs.append("fe55_waveforms.png")
    elif scan_res or incl_res:
        _plot_waveforms(scan_res, incl_res, T, out); outputs.append("fe55_waveforms.png")

    rows = []
    for name, r in scan_res.items():
        if r["n"]:
            m = scan_meta.get(name, {})
            rows.append({"point": name, "x_cm": m.get("x_cm"), "depth_mm": m.get("depth_mm"),
                         "d_wire_mm": m.get("d_wire_mm"), "n": r["n"],
                         "mean_Q": float(r["Q"].mean()), "mean_amp": float(r["amp"].mean()),
                         "mean_peak_t_ns": float(np.nanmean(r["tpeak"])),
                         "mean_rise_ns": float(np.nanmean(r["rise"])),
                         "rms_rise_ns": float(np.nanstd(r["rise"])),
                         "fwhm_ns": float(_fwhm_ns(r["t_rel"], r["aligned_mean"]))})
    if rows:
        pd.DataFrame(rows).to_csv(out / "fe55_observables.csv", index=False)
        outputs.append("fe55_observables.csv")

    print(f"done -> {out}/  ({', '.join(outputs)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
