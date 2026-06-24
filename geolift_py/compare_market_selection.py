"""
compare_market_selection.py — head-to-head: the Python market-selection inner
loop vs R, on the EXACT same combos (results/ms_R_inner_powercurves.csv from
exploration/scripts/time_market_inner_R.R). Measures:
  * wall-clock speedup (Python sequential vs R sequential; + a multiprocessing
    figure to show the production multiplier R can't easily match)
  * decision agreement: per-cell significance, ATT and scaled-L2 closeness, and
    the selected-market RANKING (minimum detectable effect per combo).
Writes exploration/results/market_selection_python_compare.json.
Run:  python geolift_py/compare_market_selection.py
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from geolift_fast.market_selection import Panel, simulate_combo

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "exploration" / "results"
ALPHA, NS, TP, SEED = 0.1, 1000, 14, 42


def mde_table(df):
    """Minimum detectable effect per combo: smallest positive EffectSize whose
    cell is significant (power==1). NaN if none."""
    pos = df[df.EffectSize > 0]
    out = {}
    for loc, g in pos.groupby("location"):
        sig = g[g.power == 1]
        out[loc] = float(sig.EffectSize.min()) if len(sig) else np.nan
    return out


def run_python(panel, combos, es_list, parallel=False):
    rng = np.random.default_rng(SEED)
    t0 = time.perf_counter()
    if parallel:
        from concurrent.futures import ProcessPoolExecutor
        import os
        seeds = np.random.SeedSequence(SEED).spawn(len(combos))
        args = [(panel.locations, panel.times, panel.Y, c, TP, es_list, NS, s)
                for c, s in zip(combos, seeds)]
        with ProcessPoolExecutor(max_workers=min(8, os.cpu_count())) as ex:
            res = list(ex.map(_combo_worker, args))
        rows = [cell for r in res for cell in r]
    else:
        rows = []
        for c in combos:
            rows += simulate_combo(panel, c, TP, es_list, NS, rng, alpha=ALPHA)
    dt = time.perf_counter() - t0
    return pd.DataFrame(rows), dt


def _combo_worker(a):
    locs, times, Y, combo, tp, es_list, ns, seed = a
    panel = Panel(locs, times, Y)
    rng = np.random.default_rng(seed)
    return simulate_combo(panel, combo, tp, es_list, ns, rng, alpha=ALPHA)


def main():
    panel = Panel.from_long_csv(ROOT / "exploration" / "data" / "ms_subset_panel.csv")
    r = pd.read_csv(RES / "ms_R_inner_powercurves.csv")
    r["location"] = r["location"].str.lower()
    es_list = sorted(r.EffectSize.unique())
    combos = [s.split(", ") for s in pd.unique(r.location)]
    r_time = json.loads((RES / "ms_R_inner_time.json").read_text())["time_s"]

    py, py_t = run_python(panel, combos, es_list, parallel=False)
    # optional multiprocessing figure (skip silently if it fails on this platform)
    try:
        _, py_mp_t = run_python(panel, combos, es_list, parallel=True)
    except Exception:
        py_mp_t = None

    j = r.merge(py, on=["location", "EffectSize"], suffixes=("_r", "_py"))
    sig_agree = float((j.power_r == j.power_py).mean())
    att_diff = float((j.att_estimator - j.AvgATT).abs().max())
    l2_diff = float((j.ScaledL2Imbalance - j.AvgScaledL2Imbalance).abs().max())
    pval_mad = float((j.pvalue_r - j.pvalue_py).abs().mean())

    mde_r, mde_py = mde_table(r), mde_table(py)
    mde_match = float(np.mean([
        (np.isnan(mde_r[l]) and np.isnan(mde_py[l])) or (mde_r[l] == mde_py[l])
        for l in mde_r]))
    # ranking by MDE (best = smallest); NaN sorts last. Top-5 overlap.
    def rank(mde):
        return [l for l, _ in sorted(mde.items(),
                key=lambda kv: (np.isnan(kv[1]), kv[1], kv[0]))]
    k = min(5, len(combos))
    top5 = len(set(rank(mde_r)[:k]) & set(rank(mde_py)[:k])) / k

    out = {
        "setup": {"seed": SEED, "n_combos": len(combos), "n_cells": len(j),
                  "effect_sizes": es_list, "ns": NS, "tp": TP, "alpha": ALPHA,
                  "panel": f"{len(panel.locations)} markets x {panel.T} periods"},
        "timing": {"r_seconds": r_time, "python_seconds": round(py_t, 3),
                   "speedup_sequential": round(r_time / py_t, 2),
                   "python_mp_seconds": None if py_mp_t is None else round(py_mp_t, 3),
                   "speedup_mp": None if py_mp_t is None else round(r_time / py_mp_t, 2),
                   "r_ms_per_cell": round(1000 * r_time / len(j), 1),
                   "python_ms_per_cell": round(1000 * py_t / len(j), 1)},
        "agreement": {"significance_agreement": round(sig_agree, 4),
                      "att_max_abs_diff": float(f"{att_diff:.3g}"),
                      "scaled_l2_max_abs_diff": float(f"{l2_diff:.3g}"),
                      "pvalue_mean_abs_diff": round(pval_mad, 4),
                      "selected_mde_agreement": round(mde_match, 4),
                      "ranking_top5_overlap": round(top5, 3)},
    }
    (RES / "market_selection_python_compare.json").write_text(json.dumps(out, indent=2))

    t = out["timing"]; a = out["agreement"]
    print(f"\n--- Python vs R market-selection inner loop ({len(j)} cells) ---")
    print(f"  R   : {t['r_seconds']:.2f}s ({t['r_ms_per_cell']} ms/cell)")
    print(f"  Py  : {t['python_seconds']:.2f}s ({t['python_ms_per_cell']} ms/cell)"
          f"  -> {t['speedup_sequential']}x sequential")
    if t["python_mp_seconds"]:
        print(f"  Py+MP: {t['python_mp_seconds']:.2f}s -> {t['speedup_mp']}x")
    print(f"  significance agreement : {100*a['significance_agreement']:.0f}%")
    print(f"  ATT max abs diff       : {a['att_max_abs_diff']:.2e}")
    print(f"  scaled-L2 max abs diff : {a['scaled_l2_max_abs_diff']:.2e}")
    print(f"  p-value mean abs diff  : {a['pvalue_mean_abs_diff']:.3f}")
    print(f"  selected-MDE agreement : {100*a['selected_mde_agreement']:.0f}%")
    print(f"  ranking top-5 overlap  : {100*a['ranking_top5_overlap']:.0f}%")
    print(f"\nWrote {RES / 'market_selection_python_compare.json'}")


if __name__ == "__main__":
    main()
