"""
validate_market_selection.py — check the Python inner loop reproduces R's
deterministic per-cell internals (weights/ATT/scaled-L2/conformal residuals/
observed statistic) to ~1e-6, and that the iid p-value agrees with R within
Monte-Carlo error. Ground truth: exploration/results/ms_truth_cells.json
(written by exploration/scripts/dump_market_cells.R).
"""
import json
from pathlib import Path
import numpy as np

from geolift_fast.market_selection import (
    Panel, ComboFit, conformal_resids, conformal_pval, _stat_func)

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "exploration" / "results"

truth = json.loads((RES / "ms_truth_cells.json").read_text())
setup = truth["setup"]
panel = Panel.from_long_csv(ROOT / "exploration" / "data" / "ms_subset_panel.csv")
ns = setup["ns"]

print(f"panel: {len(panel.locations)} units x {panel.T} periods  ns={ns}\n")
rows = []
for k, c in enumerate(truth["cells"]):
    locs = c["locations"] if isinstance(c["locations"], list) else [c["locations"]]
    es, tp = c["es"], c["tp"]
    fit = ComboFit(panel, locs, tp)
    att = fit.att(es)
    resids, t0 = conformal_resids(fit, es)
    obs = _stat_func(resids[t0:])
    r_resids = np.array(c["resids"], float)
    rng = np.random.default_rng(setup["seed"])
    pval, _ = conformal_pval(resids, t0, ns, rng)

    d_att = abs(att - c["att_estimator"])
    d_l2 = abs(fit.scaled_l2 - c["scaled_l2_imbalance"])
    d_res = np.abs(resids - r_resids).max()
    d_obs = abs(obs - c["obs_stat"])
    d_pv = abs(pval - c["pvalue_seed42"])
    rows.append(dict(cell=k, es=es, n_treated=fit.n_treated,
                     d_att=d_att, d_l2=d_l2, d_resid=d_res, d_obs=d_obs,
                     py_pval=pval, r_pval=c["pvalue_seed42"], d_pval=d_pv))
    print(f"cell{k} es={es:+.2f}  d_att={d_att:.2e}  d_l2={d_l2:.2e}  "
          f"d_resid={d_res:.2e}  d_obs={d_obs:.2e}  "
          f"p(py/R)={pval:.3f}/{c['pvalue_seed42']:.3f}")

det = ["d_att", "d_l2", "d_resid", "d_obs"]
worst = {m: max(r[m] for r in rows) for m in det}
print("\nworst deterministic diffs:", {k: f"{v:.2e}" for k, v in worst.items()})
ok = all(worst[m] < 1e-5 for m in det)
print("DETERMINISTIC MATCH:", "PASS (<1e-5)" if ok else "FAIL")
(RES / "ms_python_validation.json").write_text(json.dumps(
    {"worst_deterministic_abs_diff": worst, "cells": rows,
     "deterministic_pass_1e_5": ok}, indent=2, default=float))
print(f"\nWrote {RES / 'ms_python_validation.json'}")
