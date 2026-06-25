"""
compare_posttest.py — head-to-head equivalence + timing for the POST-TEST stage
(R GeoLift()  vs  Python geolift_fast.geolift()) on the GeoLift_Test Walkthrough.

Produces exploration/results/posttest_compare.json for the REPORT. Equivalence is
read against the frozen R ground truth (posttest_truth.json, conformal_resids.json);
Python timing is measured here; R timing is merged from posttest_R_time.json if the
R timer has been run.

    Rscript exploration/scripts/time_posttest_R.R      # writes posttest_R_time.json
    python  exploration/scripts/compare_posttest.py    # writes posttest_compare.json
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "geolift_py"))
from geolift_fast.post_test import (                       # noqa: E402
    geolift, _conformal_resids_h0, _fit, _permute_pval, _stat_total)
from geolift_fast.market_selection import Panel            # noqa: E402

RES = ROOT / "exploration/results"
truth = json.loads((RES / "posttest_truth.json").read_text())
cresid = json.loads((RES / "conformal_resids.json").read_text())
panel = Panel.from_long_csv(str(ROOT / "exploration/data/geolift_test_panel.csv"))
LOCS = ["chicago", "portland"]
trt = np.array([1 if l in LOCS else 0 for l in panel.locations])
M = panel.Y[:, :105].astype(float)
t0 = 90


def equivalence(model):
    t = truth[f"{model}_jackknife"]
    r = geolift(panel, LOCS, 91, 105, model=model, confidence_intervals=True,
                method="jackknife+", ns=1000, seed=42)
    exp_w = np.array([t["weights"][loc] for loc in r.donor_locations])
    # deterministic conformal residual vector -> limiting (infinite-ns) p-value
    lam = r.lam if model == "ridge" else None
    py_res = _conformal_resids_h0(M, trt, t0, 0.0, model, True, lam)
    r_res = np.array(cresid[model]["resids"])
    p_lim_py = _permute_pval(py_res, t0, 500_000, np.random.default_rng(0), _stat_total)
    p_lim_r = _permute_pval(r_res, t0, 500_000, np.random.default_rng(0), _stat_total)
    out = {
        "att_py": r.att, "att_r": t["att_estimator"],
        "att_rel_diff": abs(r.att - t["att_estimator"]) / abs(t["att_estimator"]),
        "percent_lift_py": r.percent_lift, "percent_lift_r": t["percent_lift"],
        "incremental_py": r.incremental, "incremental_r": t["incremental"],
        "scaled_l2_py": r.scaled_l2_imbalance, "scaled_l2_r": t["scaled_l2_imbalance"],
        "weights_max_abs_diff": float(np.max(np.abs(r.weights - exp_w))),
        "yhat_max_abs_diff": float(np.max(np.abs(r.y_hat - np.array(t["y_hat"])))),
        "jackknife_ci_py": [r.lower_conf_int, r.upper_conf_int],
        "jackknife_ci_r": [t["lower_conf_int"], t["upper_conf_int"]],
        "jackknife_ci_max_abs_diff": float(max(
            abs(r.lower_conf_int - t["lower_conf_int"]),
            abs(r.upper_conf_int - t["upper_conf_int"]))),
        "conformal_resid_max_abs_diff": float(np.max(np.abs(py_res - r_res))),
        "pvalue_limit_py": p_lim_py, "pvalue_limit_r_from_R_residuals": p_lim_r,
        "pvalue_r_ns1000_draw": t["pvalue"],
    }
    if model == "ridge":
        out["lambda_py"] = r.lam
        out["lambda_r"] = t["lambda"]
        out["lambda_rel_diff"] = abs(r.lam - t["lambda"]) / abs(t["lambda"])
    return out


def time_python(model, ci, method, reps=5):
    best = float("inf")
    for _ in range(reps):
        s = time.perf_counter()
        geolift(panel, LOCS, 91, 105, model=model, confidence_intervals=ci,
                method=method, grid_size=250, ns=1000, seed=42)
        best = min(best, time.perf_counter() - s)
    return best


py_timing = {
    "point_pval_s": round(time_python("none", False, "conformal"), 4),
    "with_ci_conformal_s": round(time_python("none", True, "conformal"), 4),
    "with_ci_jackknife_s": round(time_python("none", True, "jackknife+"), 4),
}

r_time_path = RES / "posttest_R_time.json"
r_timing = json.loads(r_time_path.read_text()) if r_time_path.exists() else None

timing = {"python": py_timing, "r": r_timing}
if r_timing:
    timing["speedup"] = {
        "point_pval": round(r_timing["point_pval_s"] / py_timing["point_pval_s"], 1),
        "with_ci_conformal": round(
            r_timing["with_ci_conformal_s"] / py_timing["with_ci_conformal_s"], 1),
    }

out = {
    "setup": {
        "example": "GeoLift Walkthrough: chicago+portland, t=91..105",
        "panel": "GeoLift_Test (40 locations x 105 periods)",
        "alpha": 0.1, "ns": 1000, "grid_size": 250, "seed": 42,
        "augsynth_commit": "06415b4 (pre-PR#88)",
        "note": ("Deterministic quantities match R to solver tolerance. The conformal "
                 "p-value is Monte-Carlo: its residual vector is identical to R, so the "
                 "limiting p-value is identical; R's ns=1000 value is one MC draw."),
    },
    "equivalence": {"none": equivalence("none"), "ridge": equivalence("ridge")},
    "timing": timing,
}
(RES / "posttest_compare.json").write_text(json.dumps(out, indent=2))
print("Wrote posttest_compare.json")
print(json.dumps(out["timing"], indent=2))
for m in ("none", "ridge"):
    e = out["equivalence"][m]
    print(f"{m:6s} ATT rel-diff {e['att_rel_diff']:.1e}  weights-d {e['weights_max_abs_diff']:.1e}  "
          f"jackCI-d {e['jackknife_ci_max_abs_diff']:.1e}  resid-d {e['conformal_resid_max_abs_diff']:.1e}  "
          f"p_lim py/R {e['pvalue_limit_py']:.4f}/{e['pvalue_limit_r_from_R_residuals']:.4f}")
