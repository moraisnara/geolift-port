"""
Fidelity test: reproduce augsynth's SCM fit (progfunc=None, scm=TRUE, fixedeff=TRUE)
for the GeoLift Walkthrough example and compare against the R benchmark.

We compare two strategies for the synthetic-control weights:
  (A) Direct simplex-constrained QP   (cvxpy + OSQP)   -- we control this
  (B) pysyncon.Synth                                   -- the off-the-shelf lib

Benchmark target: exploration/ground_truth/benchmark.json (weights, y_hat, ATT).
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

GT = Path(__file__).resolve().parents[1] / "exploration" / "ground_truth"

# ---- load augsynth's exact input matrices --------------------------------
Z0 = pd.read_csv(GT / "Z0_donors_pre.csv").to_numpy(float)      # 90 pre x 38 donors
Z1 = pd.read_csv(GT / "Z1_treated_pre.csv").to_numpy(float).ravel()   # 90
Y0 = pd.read_csv(GT / "Y0_donors_full.csv").to_numpy(float)     # 105 x 38
Y1 = pd.read_csv(GT / "Y1_treated_full.csv").to_numpy(float).ravel()  # 105
bench = json.loads((GT / "benchmark.json").read_text())["None"]
w_r = np.array(bench["weights"], float)
yhat_r = np.array(bench["y_hat"], float)
att_r = bench["ATT"]
t_int = int(bench["treatment_start"])      # 91 (1-based); pre = periods 1..90

n_pre, n_donor = Z0.shape
print(f"Z0 {Z0.shape}  Z1 {Z1.shape}  Y0 {Y0.shape}  donors={n_donor}  t_int={t_int}")

# ---- fixed effects: demean by pre-treatment unit means -------------------
mu_trt = Z1.mean()
mu_don = Z0.mean(axis=0)                    # per-donor pre mean
Z1_dm = Z1 - mu_trt
Z0_dm = Z0 - mu_don                         # broadcast over rows


def counterfactual(w):
    """Re-centered counterfactual over all 105 periods."""
    return (Y0 - mu_don) @ w + mu_trt


RESULTS = []


def metrics(w, label, fit_time, ok=True):
    yhat = counterfactual(w)
    post = slice(t_int - 1, Y1.shape[0])    # periods 91..105 (0-based 90..104)
    att = (Y1[post] - yhat[post]).mean()
    print(f"\n=== {label} ===")
    print(f"  fit time           : {fit_time*1000:8.2f} ms")
    print(f"  weights max abs diff vs R: {np.abs(w - w_r).max():.3e}")
    print(f"  y_hat  max abs diff vs R : {np.abs(yhat - yhat_r).max():.3e}")
    print(f"  ATT  python / R     : {att:.4f} / {att_r:.4f}   (diff={att-att_r:+.4f})")
    RESULTS.append({
        "strategy": label,
        "weights_max_abs_diff": float(np.abs(w - w_r).max()),
        "yhat_max_abs_diff": float(np.abs(yhat - yhat_r).max()),
        "ATT_python": round(float(att), 4),
        "ATT_R": round(float(att_r), 4),
        "fit_time_ms": round(fit_time * 1000, 2),
    })
    return att


# ---- (A) direct simplex QP -----------------------------------------------
import cvxpy as cp

w = cp.Variable(n_donor, nonneg=True)
obj = cp.Minimize(cp.sum_squares(Z1_dm - Z0_dm @ w))
prob = cp.Problem(obj, [cp.sum(w) == 1])
t0 = time.perf_counter()
prob.solve(solver=cp.CLARABEL)
t_qp = time.perf_counter() - t0
print(f"  solver status: {prob.status}")
w_qp = np.clip(np.asarray(w.value, float).ravel(), 0, None)
w_qp = w_qp / w_qp.sum()
metrics(w_qp, "(A) Direct simplex QP (cvxpy/CLARABEL)", t_qp)

# ---- (B) pysyncon --------------------------------------------------------
try:
    from pysyncon import Dataprep, Synth

    # Build long panel from the demeaned matrices is awkward; pysyncon wants raw
    # panel data. Reconstruct a tidy frame from the FULL series (donors + treated).
    don_names = [f"d{i}" for i in range(n_donor)]
    rows = []
    for t in range(Y0.shape[0]):
        for i, nm in enumerate(don_names):
            rows.append((nm, t + 1, Y0[t, i]))
        rows.append(("treated", t + 1, Y1[t]))
    df = pd.DataFrame(rows, columns=["unit", "time", "Y"])
    dp = Dataprep(
        foo=df, predictors=["Y"], predictors_op="mean",
        dependent="Y", unit_variable="unit", time_variable="time",
        treatment_identifier="treated", controls_identifier=don_names,
        time_predictors_prior=range(1, t_int),
        time_optimize_ssr=range(1, t_int),
    )
    syn = Synth()
    t0 = time.perf_counter()
    syn.fit(dataprep=dp)
    t_ps = time.perf_counter() - t0
    w_ps = np.array([syn.W[nm] for nm in don_names], float)
    metrics(w_ps, "(B) pysyncon.Synth", t_ps)
except Exception as e:
    print(f"\n=== (B) pysyncon.Synth ===\n  FAILED: {type(e).__name__}: {e}")
    RESULTS.append({"strategy": "(B) pysyncon.Synth", "status": "failed",
                    "error": f"{type(e).__name__}: {e}"})

# Persist for the report (numbers are rendered from this file, not hardcoded).
res_dir = GT.parent / "results"
res_dir.mkdir(exist_ok=True)
(res_dir / "fidelity.json").write_text(json.dumps(RESULTS, indent=2))
print(f"\nWrote {res_dir / 'fidelity.json'}")
print("Done.")
