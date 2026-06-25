"""
Self-contained unit tests for geolift_fast.market_selection.

These need NO R artifacts and NO generated data — they build tiny synthetic
panels in-memory and assert the mathematical contracts the rest of the project
relies on. Run from the repo root or from geolift_py/:

    pip install -e ".[test]"
    pytest geolift_py/tests -q
"""
import numpy as np
import pandas as pd
import pytest

from geolift_fast.market_selection import (
    Panel, ComboFit, scm_weights, conformal_resids, conformal_pval,
    power_curves, best_markets,
)

POWER_COLS = {"location", "duration", "EffectSize", "pvalue", "power",
              "AvgATT", "AvgDetectedLift", "AvgScaledL2Imbalance"}
BEST_COLS = ["rank", "location", "duration", "MDE", "AvgDetectedLift", "abs_lift_in_zero"]


def _toy_panel(seed=0, n_units=5, T=14):
    """A small wide panel with a clear donor structure."""
    rng = np.random.default_rng(seed)
    base = rng.normal(100, 5, size=(n_units, 1))
    trend = np.linspace(0, 3, T)[None, :]
    noise = rng.normal(0, 1, size=(n_units, T))
    Y = base + trend + noise
    locs = [f"m{i:02d}" for i in range(n_units)]
    return Panel(locs, np.arange(1, T + 1), Y)


# --------------------------------------------------------------------------- #
# scm_weights: simplex constraints + exact recovery                            #
# --------------------------------------------------------------------------- #
def test_scm_weights_simplex_and_recovery():
    rng = np.random.default_rng(1)
    t0, n_c = 12, 3
    Z0 = rng.normal(size=(t0, n_c))
    w_true = np.array([0.2, 0.5, 0.3])
    Z1 = Z0 @ w_true                       # Z1 is an exact convex combo of donors
    w = scm_weights(Z0, Z1)
    assert np.all(w >= -1e-9)              # non-negative
    assert w.sum() == pytest.approx(1.0, abs=1e-6)   # sums to one
    assert np.allclose(w, w_true, atol=1e-4)         # recovers the true weights


def test_scm_weights_degenerate_fallback():
    # zero donors -> uniform fallback, still a valid simplex point
    Z0 = np.zeros((5, 4))
    Z1 = np.zeros(5)
    w = scm_weights(Z0, Z1)
    assert w.sum() == pytest.approx(1.0, abs=1e-9)
    assert np.all(w >= 0)


# --------------------------------------------------------------------------- #
# ComboFit.att_lift: ATT is affine in the effect size (lever-2 invariant)       #
# --------------------------------------------------------------------------- #
def test_att_is_affine_in_effect_size():
    panel = _toy_panel(seed=2)
    fit = ComboFit(panel, ["m00"], tp=3)
    es = [0.0, 0.05, 0.10, 0.15, 0.20]
    att = np.array([fit.att(e) for e in es])
    # affine => constant first differences
    diffs = np.diff(att)
    assert np.allclose(diffs, diffs[0], atol=1e-9)
    # the pre-period fit (weights) is reused, never re-solved per es
    assert fit.w.sum() == pytest.approx(1.0, abs=1e-6)


def test_att_lift_zero_effect_matches_att():
    panel = _toy_panel(seed=3)
    fit = ComboFit(panel, ["m01"], tp=4)
    att, lift = fit.att_lift(0.07)
    assert att == pytest.approx(fit.att(0.07))
    assert np.isfinite(lift)


# --------------------------------------------------------------------------- #
# conformal_pval: determinism, bounds, degenerate residuals                     #
# --------------------------------------------------------------------------- #
def test_conformal_pval_deterministic_and_bounded():
    panel = _toy_panel(seed=4)
    fit = ComboFit(panel, ["m02"], tp=3)
    resids, t0 = conformal_resids(fit, 0.1)
    p1, obs1 = conformal_pval(resids, t0, ns=500, rng=np.random.default_rng(42))
    p2, obs2 = conformal_pval(resids, t0, ns=500, rng=np.random.default_rng(42))
    assert p1 == p2 and obs1 == obs2            # same seed -> identical
    assert 0.0 <= p1 <= 1.0


def test_conformal_pval_zero_residuals_is_one():
    resids = np.zeros(14)
    p, obs = conformal_pval(resids, t0=11, ns=100, rng=np.random.default_rng(0))
    assert obs == 0.0
    assert p == pytest.approx(1.0)              # obs <= all (zero) stats


# --------------------------------------------------------------------------- #
# power_curves: column contract + reproducibility                              #
# --------------------------------------------------------------------------- #
def test_power_curves_columns_and_seed_reproducible():
    panel = _toy_panel(seed=5, n_units=6)
    combos = [["m00", "m01"], ["m02", "m03"]]
    es = [-0.1, 0.0, 0.1]
    pc1 = power_curves(panel, combos, treatment_periods=3, effect_sizes=es, ns=200, seed=7)
    pc2 = power_curves(panel, combos, treatment_periods=3, effect_sizes=es, ns=200, seed=7)
    assert set(pc1.columns) == POWER_COLS
    assert len(pc1) == len(combos) * len(es)
    pd.testing.assert_frame_equal(pc1, pc2)     # deterministic given seed


# --------------------------------------------------------------------------- #
# best_markets: 6-column contract, negative-side MDE, es==0-only drop           #
# --------------------------------------------------------------------------- #
def _pc_row(loc, es, power, lift):
    return dict(location=loc, duration=3, EffectSize=es, pvalue=0.0 if power else 1.0,
                power=power, AvgATT=0.0, AvgDetectedLift=lift, AvgScaledL2Imbalance=0.0)


def test_best_markets_negative_side_and_contract():
    # combo A: significant on the negative side only -> MDE must be negative
    # combo B: significant only at es==0 -> dropped
    # combo C: significant on the small positive side -> positive MDE
    rows = []
    for es, p, lift in [(-0.10, 1, -0.10), (-0.05, 1, -0.05), (0.0, 0, 0.0),
                        (0.05, 0, 0.05), (0.10, 0, 0.10)]:
        rows.append(_pc_row("a", es, p, lift))
    for es, p, lift in [(-0.10, 0, -0.10), (0.0, 1, 0.0), (0.10, 0, 0.10)]:
        rows.append(_pc_row("b", es, p, lift))
    for es, p, lift in [(-0.10, 0, -0.10), (0.0, 0, 0.0), (0.05, 1, 0.05),
                        (0.10, 1, 0.10)]:
        rows.append(_pc_row("c", es, p, lift))
    pc = pd.DataFrame(rows)

    bm = best_markets(pc)
    assert list(bm.columns) == BEST_COLS
    locs = set(bm.location)
    assert "b" not in locs                       # es==0-only combo dropped
    assert {"a", "c"} <= locs
    mde = dict(zip(bm.location, bm.MDE))
    assert mde["a"] == pytest.approx(-0.05)      # closest-to-zero significant negative
    assert mde["c"] == pytest.approx(0.05)       # closest-to-zero significant positive


def test_best_markets_empty_when_nothing_significant():
    rows = [_pc_row("a", es, 0, es) for es in (-0.1, 0.0, 0.1)]
    bm = best_markets(pd.DataFrame(rows))
    assert len(bm) == 0


# --------------------------------------------------------------------------- #
# Panel.from_long_df: lowercasing + time ordering round-trip                    #
# --------------------------------------------------------------------------- #
def test_panel_from_long_df_roundtrip():
    df = pd.DataFrame({
        "location": ["Chicago", "Chicago", "Portland", "Portland"],
        "date": ["2021-01-02", "2021-01-01", "2021-01-02", "2021-01-01"],
        "Y": [11.0, 10.0, 21.0, 20.0],
    })
    panel = Panel.from_long_df(df)
    assert panel.locations == ["chicago", "portland"]     # lowercased + sorted
    # date order (not row order) drives the columns: 01-01 then 01-02
    assert panel.Y[panel.idx["chicago"]].tolist() == [10.0, 11.0]
    assert panel.Y[panel.idx["portland"]].tolist() == [20.0, 21.0]
    sub = panel.subset(["Portland"])
    assert sub.locations == ["portland"]
    assert sub.Y.shape == (1, 2)
