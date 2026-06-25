"""
bench_scaling.py — fixest-style THREE-WAY scaling benchmark:
  * R original  — GeoLift:::pvalueCalc loop, ns=1000 (exact Meta default)
  * R modified  — same engine, ns=100 (MARKETSEL_FAST preset: a PARAMETER, not a fork)
  * Python      — geolift_fast port, ns=1000 (full fidelity)

For each dataset size (number of locations L) all three run the SAME market-
selection workload (all C(L,2) treated pairs x 3 effect sizes). R timings come from
exploration/results/bench/bench_R_L*_ns*.json (bench_scaling_R.R, one fresh process
per point). This script times the Python port, joins per-cell results at the largest
size to prove the estimates agree (vs R original, same ns), writes
results/bench_scaling.json, and renders results/bench_scaling.png.

Run (after the R points exist):  python geolift_py/bench_scaling.py
"""
import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))   # find the geolift_fast package + siblings
from geolift_fast.market_selection import Panel, simulate_combo
from _compare_common import DATA, RES, BENCH, sig3, read_lower

SEED, TP, ALPHA = 42, 14, 0.1
ES = [-0.1, 0.0, 0.1]
NS_ORIG, NS_MOD = 1000, 100        # R original / R modified resample counts


def python_point(panel, L, ns):
    """Time Python over all C(L,2) pairs of the first L markets. Returns
    (time_s, cells_dataframe)."""
    markets = [f"market_{i:03d}" for i in range(1, L + 1)]
    combos = [list(c) for c in combinations(markets, 2)]
    rng = np.random.default_rng(SEED)
    t0 = time.perf_counter()
    rows = []
    for c in combos:
        rows += simulate_combo(panel, c, TP, ES, ns, rng, alpha=ALPHA)
    dt = time.perf_counter() - t0
    return dt, pd.DataFrame(rows), len(combos)


def main():
    full = Panel.from_long_csv(DATA / "realistic_panel.csv")
    r_orig = {int(p.stem.split("L")[1].split("_")[0]): json.loads(p.read_text())
              for p in BENCH.glob(f"bench_R_L*_ns{NS_ORIG}.json")}
    r_mod = {int(p.stem.split("L")[1].split("_")[0]): json.loads(p.read_text())
             for p in BENCH.glob(f"bench_R_L*_ns{NS_MOD}.json")}
    Ls = sorted(r_orig)

    rows, py_cells_by_L = [], {}
    for L in Ls:
        sub = Panel(full.locations[:L], full.times, full.Y[:L])
        py_t, py_df, n_combos = python_point(sub, L, NS_ORIG)   # Python at full ns
        py_cells_by_L[L] = py_df
        ro, rm = r_orig[L], r_mod[L]
        rows.append({
            "n_locations": L, "n_combos": n_combos, "n_cells": ro["n_cells"],
            "r_original_s": ro["time_s"], "r_modified_s": rm["time_s"],
            "python_s": round(py_t, 3),
            "speedup_py_vs_r_original": round(ro["time_s"] / py_t, 2),
            "speedup_py_vs_r_modified": round(rm["time_s"] / py_t, 2),
            "speedup_r_modified_vs_original": round(ro["time_s"] / rm["time_s"], 2)})
        print(f"L={L:2d}  R-orig={ro['time_s']:6.2f}s  R-mod(ns100)={rm['time_s']:5.2f}s  "
              f"Py={py_t:5.2f}s  ->  {rows[-1]['speedup_py_vs_r_original']}x vs orig")

    # --- results agreement at the largest size: R original vs Python (same ns) ---
    Lmax = Ls[-1]
    r_big = read_lower(BENCH / f"bench_R_L{Lmax}_ns{NS_ORIG}.csv")
    j = r_big.merge(py_cells_by_L[Lmax], on=["location", "EffectSize"],
                    suffixes=("_r", "_py"))
    agreement = {
        "at_n_locations": Lmax, "n_cells": len(j), "compared": "R original vs Python (both ns=1000)",
        "significance_agreement": round(float((j.power_r == j.power_py).mean()), 4),
        "att_max_abs_diff": sig3((j.att_estimator - j.AvgATT).abs().max()),
        "scaled_l2_max_abs_diff": sig3((j.ScaledL2Imbalance - j.AvgScaledL2Imbalance).abs().max()),
        "pvalue_mean_abs_diff": round(float((j.pvalue_r - j.pvalue_py).abs().mean()), 4)}

    out = {"setup": {"seed": SEED, "tp": TP, "effect_sizes": ES, "alpha": ALPHA,
                     "workload": "all C(L,2) treated pairs x 3 effect sizes",
                     "engines": {
                         "r_original": f"GeoLift:::pvalueCalc loop, ns={NS_ORIG} (Meta default)",
                         "r_modified": f"same engine, ns={NS_MOD} (MARKETSEL_FAST preset, parameter-only)",
                         "python": f"geolift_fast port, ns={NS_ORIG} (full fidelity)"}},
           "scaling": rows, "agreement": agreement}
    (RES / "bench_scaling.json").write_text(json.dumps(out, indent=2))

    # ------------------------------- plot ------------------------------------ #
    combos = [r["n_combos"] for r in rows]
    ro_t = [r["r_original_s"] for r in rows]
    rm_t = [r["r_modified_s"] for r in rows]
    py_t = [r["python_s"] for r in rows]
    sp_o = [r["speedup_py_vs_r_original"] for r in rows]
    sp_m = [r["speedup_py_vs_r_modified"] for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax1.plot(Ls, ro_t, "o-", color="#c0392b", lw=2, ms=7, label="R original (ns=1000)")
    ax1.plot(Ls, rm_t, "^--", color="#e67e22", lw=2, ms=7, label="R modified (ns=100 preset)")
    ax1.plot(Ls, py_t, "s-", color="#2471a3", lw=2, ms=7, label="Python port (ns=1000)")
    ax1.set_xlabel("number of locations  (workload = all C(L,2) market pairs)")
    ax1.set_ylabel("wall time (s)  ·  log scale")
    ax1.set_yscale("log")
    ax1.set_title("Market-selection inference: time vs dataset size")
    ax1.grid(True, which="both", ls=":", alpha=0.5)
    ax1.legend(frameon=False, loc="lower right", fontsize=9)
    ytop = max(ro_t) * 1.35
    for x, c in zip(Ls, combos):
        ax1.annotate(f"{c}", (x, ytop), fontsize=7, color="#999", ha="center", va="bottom")
    ax1.text(0.015, 0.97, "# market pairs evaluated", transform=ax1.transAxes,
             fontsize=7, color="#999", ha="left", va="top", style="italic")

    ax2.plot(Ls, sp_o, "D-", color="#1e8449", lw=2, ms=7, label="vs R original (same ns)")
    ax2.plot(Ls, sp_m, "v--", color="#7d3c98", lw=2, ms=7, label="vs R modified (ns=100)")
    ax2.axhline(np.mean(sp_o), ls=":", color="#1e8449", alpha=0.6)
    ax2.set_xlabel("number of locations")
    ax2.set_ylabel("Python speedup  (R time / Python time)")
    ax2.set_title("Python speedup vs dataset size")
    ax2.set_ylim(0, max(sp_o) * 1.25)
    ax2.grid(True, ls=":", alpha=0.5)
    ax2.legend(frameon=False, fontsize=9)

    fig.suptitle("GeoLift market selection — R original vs R modified vs Python port  "
                 f"({agreement['significance_agreement']*100:.0f}% significance match vs R, "
                 f"≤{agreement['att_max_abs_diff']:.0e} ATT diff)", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(RES / "bench_scaling.png", dpi=140, bbox_inches="tight")
    print(f"\nWrote {RES / 'bench_scaling.json'}")
    print(f"Wrote {RES / 'bench_scaling.png'}")
    a = agreement
    print(f"\nagreement @L={Lmax} (R-orig vs Py): sig={a['significance_agreement']*100:.0f}%  "
          f"ATT<{a['att_max_abs_diff']:.0e}  L2<{a['scaled_l2_max_abs_diff']:.0e}  "
          f"p|d|={a['pvalue_mean_abs_diff']}")
    print(f"mean speedup: vs R-original {np.mean(sp_o):.1f}x, vs R-modified {np.mean(sp_m):.1f}x")


if __name__ == "__main__":
    main()
