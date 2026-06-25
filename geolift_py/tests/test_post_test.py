"""Equivalence tests for the post-test port (geolift_fast.post_test.geolift)
against frozen R ground truth.

Deterministic quantities (weights, counterfactual, ATT, lift, incremental,
scaled-L2, ridge lambda, jackknife+ CI) must match augsynth-R to solver
tolerance. The conformal p-value is Monte-Carlo, so we instead assert the
DETERMINISTIC fact it rests on: the full-series-refit residual vector matches R,
which forces the limiting permutation p-value to be identical.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from geolift_fast.market_selection import Panel
from geolift_fast.post_test import (
    geolift, _conformal_resids_h0, _fit, _permute_pval, _stat_total)

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "exploration/results"
PANEL = Panel.from_long_csv(str(ROOT / "exploration/data/geolift_test_panel.csv"))
LOCS = ["chicago", "portland"]


def _trt_M():
    trt = np.array([1 if l in LOCS else 0 for l in PANEL.locations])
    return trt, PANEL.Y[:, :105].astype(float), 90


@pytest.fixture(scope="module")
def truth():
    return json.loads((RES / "posttest_truth.json").read_text())


@pytest.mark.parametrize("model", ["none", "ridge"])
def test_point_estimates_match_R(truth, model):
    t = truth[f"{model}_jackknife"]
    r = geolift(PANEL, LOCS, 91, 105, model=model,
                confidence_intervals=True, method="jackknife+", ns=200, seed=42)
    rel = lambda a, b: abs(a - b) / (abs(b) + 1e-12)
    assert rel(r.att, t["att_estimator"]) < 1e-5
    assert rel(r.incremental, t["incremental"]) < 1e-5
    assert rel(r.scaled_l2_imbalance, t["scaled_l2_imbalance"]) < 1e-5
    assert r.percent_lift == pytest.approx(t["percent_lift"], abs=1e-9)
    if model == "ridge":
        assert rel(r.lam, t["lambda"]) < 1e-5
    # donor weights
    exp_w = np.array([t["weights"][loc] for loc in r.donor_locations])
    assert np.max(np.abs(r.weights - exp_w)) < 1e-5
    # full counterfactual
    assert np.max(np.abs(r.y_hat - np.array(t["y_hat"]))) < 1e-3
    # deterministic jackknife+ CI (the binding CI on this data)
    assert rel(r.lower_conf_int, t["lower_conf_int"]) < 1e-5
    assert rel(r.upper_conf_int, t["upper_conf_int"]) < 1e-5


@pytest.mark.skipif(not (RES / "conformal_resids.json").exists(),
                    reason="conformal_resids.json ground truth not dumped")
@pytest.mark.parametrize("model", ["none", "ridge"])
def test_conformal_residuals_match_R(model):
    """The p-value is MC, but its residual vector is deterministic. Matching it
    proves the limiting p-value is identical to R's."""
    cr = json.loads((RES / "conformal_resids.json").read_text())[model]
    trt, M, t0 = _trt_M()
    lam = _fit(M, trt, t0, "ridge", True)[2] if model == "ridge" else None
    res = _conformal_resids_h0(M, trt, t0, 0.0, model, True, lam)
    assert np.max(np.abs(res - np.array(cr["resids"]))) < 1e-3
    # limiting permutation p-value identical whether fed py or R residuals
    p_py = _permute_pval(res, t0, 50_000, np.random.default_rng(0), _stat_total)
    p_r = _permute_pval(np.array(cr["resids"]), t0, 50_000,
                        np.random.default_rng(0), _stat_total)
    assert abs(p_py - p_r) < 5e-4


def test_conformal_falls_back_to_jackknife():
    """On this data the conformal grid retains no point, so R (and the port)
    fall back to jackknife+. The reported CI must equal the jackknife+ CI."""
    r = geolift(PANEL, LOCS, 91, 105, model="none", confidence_intervals=True,
                method="conformal", grid_size=30, ns=300, seed=1)
    assert "jackknife+" in r.ci_method_used
    assert r.lower_conf_int == pytest.approx(-81.69, abs=0.5)
    assert r.upper_conf_int == pytest.approx(378.33, abs=0.5)
