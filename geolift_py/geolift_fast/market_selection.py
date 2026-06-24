"""
geolift_fast.market_selection — fast Python port of GeoLift's market-selection /
power INNER LOOP.

This reproduces, in numpy, exactly what `GeoLift::pvalueCalc` does per simulation
cell (see GeoLift's pre_test_power.R and augsynth's conformal inference), for the
default market-selection config:
    fixed_effects = TRUE, model = "none", conformal_type = "iid",
    side_of_test = "two_sided"  (=> stat_func = sum(|x|)), alpha = 0.1.

Per cell (combo of treated locations, duration tp, effect size es):
  1. fixed-effects SCM fit on the PRE period (simplex QP, same as augsynth/osqp)
     -> weights, scaled_l2_imbalance, ATT estimator.   [es-invariant fit]
  2. conformal "iid" p-value: refit SCM over the FULL series, take the
     treated-minus-synthetic residual vector, then draw `ns` random permutations
     of it and compare the post-period statistic.       [vectorized over ns]

The structural win vs R: R draws the ns permutations with an interpreted
`sapply(1:ns, ...)`; numpy draws all ns at once as a matrix. And because the
PRE-period fit is invariant to es (es only inflates treated POST values), we fit
it ONCE per combo and reuse it across every effect size (GeoLift's "lever 2").

Public API (see also geolift_fast.__init__):
    Panel, ComboFit, simulate_combo        — building blocks
    power_curves, best_markets             — high-level, DataFrame in/out
"""
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import sparse
import osqp


# --------------------------------------------------------------------------- #
# Simplex-constrained SCM weights:  min ||Z1 - Z0 w||^2  s.t. w>=0, sum w = 1   #
# Z0: (t x n_c) donor matching matrix, Z1: (t,) treated target. Uses OSQP, the  #
# same solver augsynth calls, so weights match to solver tolerance.            #
# --------------------------------------------------------------------------- #
def scm_weights(Z0, Z1):
    n_c = Z0.shape[1]
    P = sparse.csc_matrix(2.0 * (Z0.T @ Z0))
    q = -2.0 * (Z0.T @ Z1)
    # constraints: sum(w)=1 ; w>=0
    A = sparse.vstack([sparse.csc_matrix(np.ones((1, n_c))), sparse.eye(n_c)]).tocsc()
    lo = np.concatenate([[1.0], np.zeros(n_c)])
    up = np.concatenate([[1.0], np.full(n_c, np.inf)])
    m = osqp.OSQP()
    m.setup(P=P, q=q, A=A, l=lo, u=up, verbose=False,
            eps_abs=1e-9, eps_rel=1e-9, max_iter=20000, polish=True)
    res = m.solve()
    w = np.clip(np.asarray(res.x, float), 0, None)
    s = w.sum()
    return w / s if s > 0 else np.full(n_c, 1.0 / n_c)


def _scaled_l2(Z0, Z1, w):
    n_c = Z0.shape[1]
    l2 = np.sqrt(np.sum((Z0 @ w - Z1) ** 2))
    uni = np.full(n_c, 1.0 / n_c)
    unif = np.sqrt(np.sum((Z0 @ uni - Z1) ** 2))
    return l2 / unif if unif > 0 else np.nan


def _stat_func(x):
    # two_sided test statistic used by GeoLift's market selection
    return np.abs(x).sum()


class Panel:
    """Wide outcome matrix Y (n_units x T) plus a location index."""
    def __init__(self, locations, times, Y):
        self.locations = list(locations)
        self.times = np.asarray(times)
        self.Y = np.asarray(Y, float)          # (n_units, T)
        self.idx = {loc: i for i, loc in enumerate(self.locations)}
        self.T = self.Y.shape[1]

    @classmethod
    def from_long_df(cls, df, loc_col="location", time_col="date", y_col="Y"):
        """Build a Panel from a long/tidy DataFrame (one row per location-period).
        Locations are lowercased and time is ordered by sorted unique value, matching
        GeoLift's GeoDataRead (which assigns time = 1..T in date order)."""
        df = df.copy()
        df[loc_col] = df[loc_col].astype(str).str.lower()
        tvals = sorted(df[time_col].unique())
        tpos = {t: i for i, t in enumerate(tvals)}
        locs = sorted(df[loc_col].unique())
        lpos = {l: i for i, l in enumerate(locs)}
        Y = np.zeros((len(locs), len(tvals)))
        Y[df[loc_col].map(lpos).to_numpy(), df[time_col].map(tpos).to_numpy()] = \
            df[y_col].to_numpy(float)
        return cls(locs, np.arange(1, len(tvals) + 1), Y)

    @classmethod
    def from_long_csv(cls, path, loc_col="location", time_col="date", y_col="Y"):
        return cls.from_long_df(pd.read_csv(path), loc_col, time_col, y_col)

    def subset(self, locations):
        """A new Panel restricted to `locations` (lowercased), order preserved."""
        wanted = [l.lower() for l in locations]
        rows = [self.idx[l] for l in wanted]
        return Panel(wanted, self.times, self.Y[rows])


class ComboFit:
    """The es-invariant PRE-period fit for one combo of treated locations,
    cached so every effect size reuses it (lever 2)."""
    def __init__(self, panel: Panel, treated, tp):
        self.tp = tp
        self.T = panel.T
        self.t0 = panel.T - tp                          # pre length (post = last tp)
        ti = [panel.idx[l.lower()] for l in treated]
        self.trt = np.zeros(panel.Y.shape[0], bool); self.trt[ti] = True
        self.Y = panel.Y
        self.don = ~self.trt
        # --- pre-period fixed-effects SCM fit (es-invariant) ---
        Xpre = panel.Y[:, :self.t0]                     # (units, t0)
        mu = Xpre.mean(axis=1)                           # per-unit pre mean
        self.mu = mu
        Z0 = (Xpre[self.don] - mu[self.don, None]).T     # (t0, n_c)
        Z1 = (Xpre[self.trt] - mu[self.trt, None]).mean(axis=0)  # (t0,)
        self.w = scm_weights(Z0, Z1)
        self.scaled_l2 = _scaled_l2(Z0, Z1, self.w)
        self.n_treated = int(self.trt.sum())

    def att(self, es):
        """ATT estimator for this combo at effect size es (pvalueCalc formula)."""
        tp, t0 = self.tp, self.t0
        Ypost = self.Y[:, t0:].copy()                    # (units, tp)
        Ypost[self.trt] *= (1.0 + es)                    # inflate treated post
        ave_treatment = Ypost[self.trt].sum() / self.n_treated
        # synthetic control over post periods (fixed-effects predict)
        mu_d = self.mu[self.don]
        m1 = self.mu[self.trt].mean()
        y0_post = m1 + (Ypost[self.don] - mu_d[:, None]).T @ self.w   # (tp,)
        ave_incremental = ave_treatment - y0_post.sum()
        return ave_incremental / tp


def conformal_resids(fit: ComboFit, es):
    """Residual vector for the iid conformal test: refit SCM over the FULL
    series (treated post inflated by es) and return treated-minus-synthetic
    over all T periods, plus t0 (pre length)."""
    Y = fit.Y.copy()
    Y[fit.trt, fit.t0:] *= (1.0 + es)                    # inflate treated post
    mu = Y.mean(axis=1)                                  # full-series unit means
    Z0 = (Y[fit.don] - mu[fit.don, None]).T              # (T, n_c)
    Z1 = (Y[fit.trt] - mu[fit.trt, None]).mean(axis=0)   # (T,)
    w = scm_weights(Z0, Z1)
    m1 = mu[fit.trt].mean()
    y0 = m1 + (Y[fit.don] - mu[fit.don, None]).T @ w     # (T,)
    treated_mean = Y[fit.trt].mean(axis=0)               # (T,)
    return treated_mean - y0, fit.t0


def conformal_pval(resids, t0, ns, rng):
    """iid permutation p-value: shuffle the full residual vector ns times,
    take the post slice, compare the statistic. Fully vectorized."""
    post = resids[t0:]
    obs = _stat_func(post)
    n_post = post.size
    T = resids.size
    # ns independent permutations of the T residuals, post slice = first n_post
    perm = np.argsort(rng.random((ns, T)), axis=1)[:, :n_post]
    stats = np.abs(resids[perm]).sum(axis=1)
    # R: mean(obs <= test_stats)
    return float(np.mean(obs <= stats)), float(obs)


def simulate_combo(panel: Panel, treated, tp, effect_sizes, ns, rng, alpha=0.1):
    """One combo across all effect sizes. Returns a list of cell dicts.
    The PRE fit is computed once and reused (lever 2)."""
    fit = ComboFit(panel, treated, tp)
    out = []
    for es in effect_sizes:
        resids, t0 = conformal_resids(fit, es)
        pval, obs = conformal_pval(resids, t0, ns, rng)
        out.append(dict(
            location=", ".join(l.lower() for l in treated),
            duration=tp, EffectSize=es,
            pvalue=pval, power=int(pval < alpha),
            AvgATT=fit.att(es), AvgScaledL2Imbalance=fit.scaled_l2))
    return out


# --------------------------------------------------------------------------- #
# High-level, DataFrame-in / DataFrame-out API (the "library" surface).        #
# --------------------------------------------------------------------------- #
def power_curves(data, locations, treatment_periods, effect_sizes,
                 ns=1000, alpha=0.1, seed=42):
    """Run the market-selection power simulation and return a tidy PowerCurves
    DataFrame — the Python analogue of `GeoLift::GeoLiftMarketSelection(...)$PowerCurves`.

    Parameters
    ----------
    data : pandas.DataFrame | str | Panel
        Long/tidy panel (columns location, date, Y), a path to such a CSV, or a Panel.
    locations : list[list[str]]
        Candidate treated-market combinations to evaluate, e.g. [["chicago","portland"], ...].
    treatment_periods : int | list[int]
        Test-window length(s) tp (number of post-treatment periods).
    effect_sizes : list[float]
        Simulated lifts, e.g. [-0.1, -0.05, 0.0, 0.05, 0.1].
    ns : int
        Conformal permutation count (Monte-Carlo p-value precision). Default 1000 (Meta default).
    alpha : float
        Significance level for the `power` flag. Default 0.1.
    seed : int
        RNG seed for the permutation draws (statistics deterministic given seed).

    Returns
    -------
    pandas.DataFrame with columns:
        location, duration, EffectSize, pvalue, power, AvgATT, AvgScaledL2Imbalance
    """
    if isinstance(data, Panel):
        panel = data
    elif isinstance(data, (str,)) or hasattr(data, "__fspath__"):
        panel = Panel.from_long_csv(data)
    else:
        panel = Panel.from_long_df(data)
    tps = [treatment_periods] if np.isscalar(treatment_periods) else list(treatment_periods)
    rng = np.random.default_rng(seed)
    rows = []
    for combo in locations:
        for tp in tps:
            rows += simulate_combo(panel, combo, tp, effect_sizes, ns, rng, alpha=alpha)
    return pd.DataFrame(rows)


def all_pairs(panel, size=2, locations=None):
    """All location combinations of `size` (default pairs) — a convenience for
    building the `locations` argument of power_curves. `locations` restricts the pool."""
    pool = [l.lower() for l in (locations if locations is not None else panel.locations)]
    return [list(c) for c in combinations(pool, size)]


def best_markets(pc, alpha=0.1):
    """Rank candidate markets by minimum detectable effect (MDE) — the smallest
    positive EffectSize that is significant — best (smallest MDE) first. Returns a
    DataFrame: location, duration, MDE (NaN if never significant), rank."""
    pos = pc[pc.EffectSize > 0]
    recs = []
    for (loc, dur), g in pos.groupby(["location", "duration"]):
        sig = g[g.pvalue < alpha]
        recs.append({"location": loc, "duration": dur,
                     "MDE": float(sig.EffectSize.min()) if len(sig) else np.nan})
    out = pd.DataFrame(recs).sort_values(
        ["MDE", "location"], na_position="last").reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    return out
