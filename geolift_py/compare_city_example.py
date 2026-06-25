"""
compare_city_example.py — head-to-head on the CANONICAL GeoLift example data
(GeoLift_Test, real US city names). Runs the Python market-selection engine on the
EXACT combos R's GeoLiftMarketSelection chose (from build_city_example.R) and
measures agreement: per-cell significance / ATT / scaled-L2 / detected-lift, plus
the selected best-market set and per-combo MDE.

Inputs  (produced by exploration/scripts/build_city_example.R):
    exploration/data/geolift_test_panel.csv
    exploration/results/citylift_R_powercurves.csv
    exploration/results/citylift_R_bestmarkets.csv
Output:
    exploration/results/citylift_python_compare.json
Run:  python geolift_py/compare_city_example.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from geolift_fast import Panel, power_curves, best_markets

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "exploration" / "data"
RES = ROOT / "exploration" / "results"
ALPHA, NS, TP, SEED = 0.10, 1000, 14, 42


def main():
    panel = Panel.from_long_csv(DATA / "geolift_test_panel.csv")
    r_pc = pd.read_csv(RES / "citylift_R_powercurves.csv")
    r_pc["location"] = r_pc["location"].str.lower()
    r_bm = pd.read_csv(RES / "citylift_R_bestmarkets.csv")
    r_bm["location"] = r_bm["location"].str.lower()

    combos = [s.split(", ") for s in pd.unique(r_pc.location)]
    es = sorted(r_pc.EffectSize.unique())

    py_pc = power_curves(panel, combos, treatment_periods=TP, effect_sizes=es,
                         ns=NS, alpha=ALPHA, seed=SEED)

    j = r_pc.merge(py_pc, on=["location", "duration", "EffectSize"], suffixes=("_r", "_py"))
    per_cell = {
        "n_cells": int(len(j)),
        "significance_agreement": round(float((j.power_r == j.power_py).mean()), 4),
        "att_max_abs_diff": float(f"{(j.AvgATT_r - j.AvgATT_py).abs().max():.3g}"),
        "scaled_l2_max_abs_diff": float(f"{(j.AvgScaledL2Imbalance_r - j.AvgScaledL2Imbalance_py).abs().max():.3g}"),
        "detected_lift_max_abs_diff": float(f"{(j.AvgDetectedLift_r - j.AvgDetectedLift_py).abs().max():.3g}"),
        "pvalue_note": "R's $PowerCurves does not expose the p-value; only the significance flag is compared.",
    }

    py_bm = best_markets(py_pc)
    py_mde = dict(zip(py_bm.location, py_bm.MDE))
    r_mde = dict(zip(r_bm.location, r_bm.EffectSize))
    common = sorted(set(py_mde) & set(r_mde))
    same_set = set(py_bm.location) == set(r_bm.location)

    def topk_overlap(k):
        a = set(py_bm.sort_values("rank").location.head(k))
        b = set(r_bm.sort_values("rank").location.head(k))
        return round(len(a & b) / k, 3)

    best = {
        "n_markets_python": int(len(py_bm)),
        "n_markets_r": int(len(r_bm)),
        "selected_set_identical": bool(same_set),
        "mde_magnitude_agreement": round(float(np.mean([abs(py_mde[l]) == abs(r_mde[l]) for l in common])), 4),
        "mde_signed_agreement": round(float(np.mean([py_mde[l] == r_mde[l] for l in common])), 4),
        "top5_set_overlap": topk_overlap(5),
        "top10_set_overlap": topk_overlap(10),
        "note": ("Per-combo MDE magnitudes match; the few signed differences are direction-only "
                 "at equal magnitude (both ±es significant), and exact integer rank within a tie "
                 "group depends on a 3rd-decimal tiebreak — neither moves the selected market set."),
    }

    out = {
        "setup": {"data": "GeoLift_Test (40 US cities x 105 days)", "seed": SEED,
                  "combos": len(combos), "effect_sizes": es, "ns": NS, "tp": TP, "alpha": ALPHA},
        "per_cell": per_cell,
        "best_markets": best,
    }
    (RES / "citylift_python_compare.json").write_text(json.dumps(out, indent=2))

    print("--- Python vs R on GeoLift_Test (city names) ---")
    print(f"  cells: {per_cell['n_cells']}  significance agreement: "
          f"{100*per_cell['significance_agreement']:.1f}%")
    print(f"  ATT max abs diff {per_cell['att_max_abs_diff']:.2e}   "
          f"scaledL2 {per_cell['scaled_l2_max_abs_diff']:.2e}   "
          f"lift {per_cell['detected_lift_max_abs_diff']:.2e}")
    print(f"  best-market set identical: {best['selected_set_identical']}   "
          f"MDE |mag| agree: {100*best['mde_magnitude_agreement']:.0f}%   "
          f"top5/top10 overlap: {100*best['top5_set_overlap']:.0f}%/{100*best['top10_set_overlap']:.0f}%")
    print(f"\nWrote {RES / 'citylift_python_compare.json'}")


if __name__ == "__main__":
    main()
