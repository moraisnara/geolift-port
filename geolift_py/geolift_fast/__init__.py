"""
geolift_fast — a fast Python port of GeoLift, both halves of the workflow.

Pre-test (market selection / power): reproduces GeoLift's per-cell decisions (ATT,
scaled-L2, conformal iid p-value, selected-market ranking) to solver tolerance, while
running ~10-14x faster than the R inner loop single-threaded — by vectorising the
permutation draw and reusing the effect-size-invariant SCM fit across effect sizes.

Post-test (single-test measurement): `geolift()` mirrors R's `GeoLift()` — augsynth-faithful
fixed-effects / ridge SCM, conformal p-value, conformal-CI grid inversion with jackknife+
fallback. Deterministic quantities match R to solver tolerance; the conformal p-value matches
in the limit (identical residual vector). See post_test.py.

Quick start
-----------
    import pandas as pd
    from geolift_fast import power_curves, best_markets, all_pairs, Panel

    df = pd.read_csv("panel.csv")            # columns: location, date, Y
    panel = Panel.from_long_df(df)
    combos = all_pairs(panel, size=2)        # or your own list of treated-market lists
    pc = power_curves(panel, combos, treatment_periods=14,
                      effect_sizes=[-0.1, -0.05, 0.0, 0.05, 0.1], ns=1000, seed=42)
    ranking = best_markets(pc)               # markets ranked by minimum detectable effect
"""
from .market_selection import (
    Panel,
    ComboFit,
    simulate_combo,
    conformal_resids,
    conformal_pval,
    scm_weights,
    power_curves,
    best_markets,
    all_pairs,
)
from .post_test import geolift, GeoLiftResult

__all__ = [
    "Panel", "ComboFit", "simulate_combo", "conformal_resids", "conformal_pval",
    "scm_weights", "power_curves", "best_markets", "all_pairs",
    "geolift", "GeoLiftResult",
]
__version__ = "0.1.0"
