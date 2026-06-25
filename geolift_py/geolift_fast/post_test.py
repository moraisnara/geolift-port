"""
post_test.py — Python port of GeoLift's POST-TEST estimation/inference.

R's `GeoLift()` (post_test_analysis.R) is a thin orchestration layer that
delegates the synthetic-control fit to the **augsynth** package and adds
GeoLift's conformal permutation inference on top. This module mirrors that
architecture faithfully: it re-implements augsynth's fixed-effects (de-meaned)
SCM and ridge-augmented SCM exactly, then layers GeoLift's conformal p-value,
conformal confidence-interval grid inversion, and jackknife+ confidence
intervals.

Equivalence to R (the design goal):
  * Deterministic quantities — SCM weights, counterfactual (`y_hat`), ATT, the
    ATT series, Perc.Lift, Incremental, scaled-L2 imbalance, ridge lambda, and
    the **jackknife+** CI — match augsynth-R to solver tolerance (~1e-6).
  * Monte-Carlo quantities — the conformal "iid" p-value and the conformal CI
    bounds — match R only in expectation (independent R/numpy RNG); they differ
    by Monte-Carlo error, exactly as two R runs with different seeds would.

Faithful to augsynth commit 06415b4 (the version GeoLift 2.7.5 pins). Mapped
source: fit_augsynth_internal, fit_ridgeaug_formatted/_inner, demean_data,
predict.augsynth, conformal_inf, compute_permute_pval/_test_stats,
time_jackknife_plus, drop_time_t; GeoLift ConfIntervals + type_of_test.

NOT ported: progfunc="GSYN" — augsynth delegates that to the separate `gsynth`
factor-model package (fit_prog_gsynth just calls gsynth::gsynth), so reproducing
it bit-for-bit would mean re-porting a third package. `model="best"` is supported
only insofar as it chooses among {none, ridge} by scaled-L2 (GSYN excluded).
"""
from dataclasses import dataclass, field

import numpy as np

from .market_selection import Panel, scm_weights


# --------------------------------------------------------------------------- #
# Test statistics — GeoLift::type_of_test (pre_test_power.R).                  #
# --------------------------------------------------------------------------- #
def _stat_total(x):       # two_sided: sum(abs(x))
    return float(np.abs(x).sum())


def _stat_negative(x):    # one_sided, H0: ES>=0  ->  -sum(x)
    return float(-np.sum(x))


def _stat_positive(x):    # one_sided, H0: ES<=0  ->  sum(x)
    return float(np.sum(x))


_STAT = {"total": _stat_total, "negative": _stat_negative, "positive": _stat_positive}


# --------------------------------------------------------------------------- #
# Faithful augsynth fit (progfunc in {"none","ridge"}, fixedeff, scm=TRUE).    #
# Operates on a wide outcome block split into `t0` "pre" columns (matched) and #
# the remaining "post" columns (predicted). Returns weights + the de-mean      #
# offsets, from which the counterfactual / ATT path is reconstructed exactly   #
# as predict.augsynth does.                                                    #
# --------------------------------------------------------------------------- #
def _fit(M, trt, t0, progfunc="none", fixedeff=True, lam=None,
         lambda_min_ratio=1e-8, n_lambda=20):
    """Fit the (augmented) SCM on wide matrix M (n_units x ttot), matching the
    first `t0` columns. Returns (weights, means, lam_used).

    `weights` are the donor weights (n_c,), `means` the per-unit de-mean offset
    (n,) added back by predict (0 if not fixedeff). For ridge, `weights` are the
    SCM weights plus augsynth's closed-form ridge correction."""
    don = trt == 0
    Xpre = M[:, :t0]
    if fixedeff:
        means = Xpre.mean(axis=1)                      # demean_data: rowMeans(X)
    else:
        means = np.zeros(M.shape[0])
    Xd = Xpre - means[:, None]
    # fit_ridgeaug_formatted: per-column control-mean centering (cancels under
    # sum-to-1 weights, but kept so the ridge correction matches augsynth exactly)
    col_c = Xd[don].mean(axis=0)
    X_c = Xd[don] - col_c[None, :]                     # (n_c, t0)
    X_1 = (Xd[trt == 1] - col_c[None, :]).mean(axis=0) # (t0,)
    syn = scm_weights(X_c.T, X_1)                      # simplex SCM weights (n_c,)
    weights = syn.copy()
    lam_used = None
    if progfunc == "ridge":
        if lam is None:
            lam = _cv_lambda(X_c, X_1, lambda_min_ratio, n_lambda)
        lam_used = float(lam)
        # fit_ridgeaug_inner: ridge_w = (X_1 - X_c' syn)' (X_c'X_c + lam I)^-1 X_c'
        imb = X_1 - X_c.T @ syn                        # (t0,)
        G = X_c.T @ X_c + lam * np.eye(t0)             # (t0, t0)
        ridge_w = imb @ np.linalg.solve(G, X_c.T)      # (n_c,)
        weights = syn + ridge_w
    return weights, means, lam_used


def _counterfactual(M, trt, means, weights):
    """predict.augsynth(att=FALSE): full-timeline counterfactual for the treated
    aggregate.  y0[t] = mean(means[trt]) + sum_j w_j (M[control_j, t] - means_j)."""
    don = trt == 0
    m1 = means[trt == 1].mean()
    resid = M[don] - means[don, None]                  # (n_c, ttot)
    return m1 + resid.T @ weights                      # (ttot,)


def _att_path(M, trt, means, weights):
    """predict.augsynth(att=TRUE): treated-aggregate minus counterfactual."""
    return M[trt == 1].mean(axis=0) - _counterfactual(M, trt, means, weights)


# --------------------------------------------------------------------------- #
# Ridge lambda cross-validation (augsynth cv_lambda, scm=TRUE, min_1se=TRUE,   #
# holdout_length=1). Deterministic; reproduces augsynth's choice exactly.      #
# --------------------------------------------------------------------------- #
def _cv_lambda(X_c, X_t, lambda_min_ratio, n_lambda):
    """X_c is (n_c, t0) centered controls, X_t is (t0,) centered treated mean.
    Mirrors augsynth cv_lambda exactly: leave-one-PRE-PERIOD-out holdout MSE of
    the ridge-augmented fit, over the first t0-1 periods (holdout_length=1)."""
    # get_lambda_max: largest singular value^2 of X_c
    lambda_max = float(np.linalg.svd(X_c, compute_uv=False)[0] ** 2)
    # create_lambda_list: lambda_max * scaler^(0..n_lambda), scaler = ratio^(1/n_lambda)
    scaler = lambda_min_ratio ** (1.0 / n_lambda)
    lambdas = lambda_max * scaler ** np.arange(n_lambda + 1)   # length n_lambda+1
    n_c, t0 = X_c.shape
    n_hold = t0 - 1                                            # ncol(X_c) - holdout_length
    errors = np.zeros((n_hold, len(lambdas)))
    for i in range(n_hold):                                    # drop pre-period i
        keep = np.ones(t0, bool); keep[i] = False
        X0 = X_c[:, keep]                                      # (n_c, t0-1)
        X1 = X_t[keep]                                         # (t0-1,)
        X0v, X1v = X_c[:, i], X_t[i]                           # held-out column / scalar
        syn = scm_weights(X0.T, X1)                            # refit synth (scm=TRUE)
        imb = X1 - X0.T @ syn
        G = X0.T @ X0
        for j, lam in enumerate(lambdas):
            ridge_w = imb @ np.linalg.solve(G + lam * np.eye(t0 - 1), X0.T)
            aug = syn + ridge_w
            errors[i, j] = float((X1v - X0v @ aug) ** 2)
    mean_err = errors.mean(axis=0)
    se_err = errors.std(axis=0, ddof=1) / np.sqrt(n_hold)
    # choose_lambda, min_1se=TRUE: largest lambda within 1 SE of the min-error lambda
    k = int(np.argmin(mean_err))
    cand = lambdas[mean_err <= mean_err[k] + se_err[k]]
    return float(cand.max())


# --------------------------------------------------------------------------- #
# Conformal permutation p-value (augsynth compute_permute_pval), generalized   #
# to an arbitrary null h0. Refits the SCM over the FULL series with h0 removed  #
# from the treated post block, then permutes residuals (iid).                  #
# --------------------------------------------------------------------------- #
def _conformal_resids_h0(M, trt, t0, h0, progfunc, fixedeff, lam):
    """Residual path under H0: ATT==h0 (constant). Subtract h0 from treated post,
    refit over the whole series, return full-length att path."""
    Madj = M.copy()
    Madj[trt == 1, t0:] -= h0
    ttot = Madj.shape[1]
    w, means, _ = _fit(Madj, trt, ttot, progfunc, fixedeff, lam)  # match ALL periods
    return _att_path(Madj, trt, means, w)


def _permute_pval(resids, t0, ns, rng, stat_func):
    """iid permutation p-value: mean(stat(obs_post) <= stat(shuffled_post))."""
    post = resids[t0:]
    obs = stat_func(post)
    T, n_post = resids.size, post.size
    perm = np.argsort(rng.random((ns, T)), axis=1)[:, :n_post]
    stats = np.array([stat_func(row) for row in resids[perm]]) if stat_func is not _stat_total \
        else np.abs(resids[perm]).sum(axis=1)
    return float(np.mean(obs <= stats))


# --------------------------------------------------------------------------- #
# Confidence intervals — GeoLift::ConfIntervals.                              #
# --------------------------------------------------------------------------- #
def _conformal_ci(M, trt, t0, alpha, grid_size, ns, rng, stat_func,
                  progfunc, fixedeff, lam):
    """Grid-inversion conformal CI on the average post ATT (ConfIntervals,
    method='conformal'). Returns (lo, hi); (inf, -inf) if no grid point is
    retained (the caller then falls back to jackknife+, as R does)."""
    att = _att_path(M, trt, *_fit_mw(M, trt, t0, progfunc, fixedeff, lam))
    post_att = att[t0:]
    post_sd = np.sqrt(np.mean(post_att ** 2))
    mean_att = post_att.mean()
    grid = np.linspace(mean_att - 6 * post_sd, mean_att + 6 * post_sd, grid_size)
    grid = np.concatenate([[0.0], grid])
    keep = []
    for null in grid:
        resids = _conformal_resids_h0(M, trt, t0, null, progfunc, fixedeff, lam)
        p = _permute_pval(resids, t0, ns, rng, stat_func)
        if p >= alpha:
            keep.append(null)
    if not keep:
        return np.inf, -np.inf
    return float(min(keep)), float(max(keep))


def _fit_mw(M, trt, t0, progfunc, fixedeff, lam):
    w, means, _ = _fit(M, trt, t0, progfunc, fixedeff, lam)
    return means, w


def _jackknife_plus_ci(M, trt, t0, alpha, progfunc, fixedeff, lam):
    """Deterministic jackknife+ CI on the average post ATT
    (augsynth time_jackknife_plus, non-conservative). Matches R exactly."""
    n, ttot = M.shape
    tpost = ttot - t0
    w, means, _ = _fit(M, trt, t0, progfunc, fixedeff, lam)
    att = _att_path(M, trt, means, w)
    # leave-one-pre-period-out refits
    est_plus, est_minus = [], []        # each row: post estimates + their mean (len tpost+1)
    for td in range(t0):
        cols = [c for c in range(ttot) if c != td]      # drop pre col td; it becomes held-out
        order = cols[:t0 - 1] + [td] + cols[t0 - 1:]    # X_new (t0-1) | [dropped, post...]
        Mr = M[:, order]
        wr, mr, _ = _fit(Mr, trt, t0 - 1, progfunc, fixedeff, lam)
        y0_full = _counterfactual(Mr, trt, mr, wr)      # length ttot
        est_post = y0_full[t0:]                          # counterfactual over true post
        err = float(M[trt == 1, td].mean() - y0_full[t0 - 1])  # held-out error at dropped col
        est = np.append(est_post, est_post.mean())       # len tpost+1
        est_plus.append(est + abs(err))
        est_minus.append(est - abs(err))
    est_plus = np.array(est_plus)        # (t0, tpost+1)
    est_minus = np.array(est_minus)
    lb_raw = np.quantile(est_minus, alpha / 2, axis=0, method="linear")
    ub_raw = np.quantile(est_plus, 1 - alpha / 2, axis=0, method="linear")
    # shift onto the ATT scale: y1 = observed treated; CI = y1 - [ub_raw, lb_raw]
    y0 = _counterfactual(M, trt, means, w)
    y1 = M[trt == 1].mean(axis=0)                        # = y0 + att, length ttot
    y1_agg = y1[t0:].mean()
    # aggregate (average post ATT) bound is the last element (index tpost)
    lo = float(y1_agg - ub_raw[-1])
    hi = float(y1_agg - lb_raw[-1])
    return lo, hi


# --------------------------------------------------------------------------- #
# Public API — geolift(): the post-test analogue of R's GeoLift().             #
# --------------------------------------------------------------------------- #
@dataclass
class GeoLiftResult:
    att: float                      # average post-treatment ATT (per period)
    percent_lift: float
    pvalue: float
    incremental: float
    lower_conf_int: float = None
    upper_conf_int: float = None
    l2_imbalance: float = None
    scaled_l2_imbalance: float = None
    model: str = "none"
    lam: float = None
    ci_method_used: str = None
    treatment_start: int = None
    treatment_end: int = None
    weights: np.ndarray = field(default=None, repr=False)
    donor_locations: list = field(default=None, repr=False)
    y_hat: np.ndarray = field(default=None, repr=False)      # counterfactual, full timeline
    att_series: np.ndarray = field(default=None, repr=False)

    def summary(self):
        ci = (f"[{self.lower_conf_int:.2f}, {self.upper_conf_int:.2f}] "
              f"({self.ci_method_used})" if self.lower_conf_int is not None else "(not requested)")
        return (f"GeoLift post-test [{self.model}]  t={self.treatment_start}..{self.treatment_end}\n"
                f"  ATT (avg/period) : {self.att:.4f}\n"
                f"  Percent Lift     : {self.percent_lift:.2f}%\n"
                f"  Incremental      : {self.incremental:.2f}\n"
                f"  p-value          : {self.pvalue:.4f}\n"
                f"  Scaled L2        : {self.scaled_l2_imbalance:.5f}\n"
                f"  Confidence Int.  : {ci}")


def geolift(data, locations, treatment_start_time, treatment_end_time,
            alpha=0.1, model="none", fixed_effects=True,
            confidence_intervals=False, method="conformal", grid_size=250,
            stat_test="Total", ns=1000, seed=42):
    """Post-test geo-lift inference — the Python analogue of `GeoLift::GeoLift()`.

    Parameters mirror R. `data` is a long/tidy DataFrame, CSV path, or Panel
    (columns location/date/Y). `locations` is the list of treated locations.
    Times are 1-based indices (pre = 1..start-1, post = start..end). `model` in
    {"none","ridge","best"}; `stat_test` in {"Total","Negative","Positive"};
    `method` in {"conformal","jackknife+"}.

    Returns a GeoLiftResult. Deterministic fields match augsynth-R to ~1e-6; the
    p-value and (conformal) CI differ from R only by Monte-Carlo error.
    """
    panel = data if isinstance(data, Panel) else (
        Panel.from_long_csv(data) if isinstance(data, str) or hasattr(data, "__fspath__")
        else Panel.from_long_df(data))
    locs = [l.lower() for l in locations]
    ts, te = int(treatment_start_time), int(treatment_end_time)
    t0 = ts - 1                                   # pre length
    ttot = te                                     # R filters time <= treatment_end_time
    trt = np.array([1 if l in locs else 0 for l in panel.locations])
    M = panel.Y[:, :ttot].astype(float)
    stat_func = _STAT[stat_test.lower()]

    # model="best": pick none vs ridge by scaled-L2 (GSYN excluded — see module docstring)
    if model == "best":
        cands = {}
        for m in ("none", "ridge"):
            w, means, lam = _fit(M, trt, t0, m, fixed_effects)
            cands[m] = (_scaled_l2(M, trt, t0, means, w), w, means, lam)
        model = min(cands, key=lambda m: cands[m][0])
        scaled_l2, weights, means, lam = cands[model]
    else:
        weights, means, lam = _fit(M, trt, t0, model, fixed_effects)
        scaled_l2 = _scaled_l2(M, trt, t0, means, weights)

    y_hat = _counterfactual(M, trt, means, weights)              # full timeline
    att_series = M[trt == 1].mean(axis=0) - y_hat
    post_att = att_series[t0:]
    att = float(post_att.mean())

    # Lift / Incremental — GeoLift formulas (multi-location uses treated colMeans)
    n_treated = int(trt.sum())
    treated_post_mean = M[trt == 1][:, t0:].mean(axis=0)         # colMeans over treated
    pred_post = y_hat[t0:]
    lift = (treated_post_mean.sum() - pred_post.sum()) / abs(pred_post.sum())
    treated_post_sum = M[trt == 1][:, t0:].sum()
    incremental = float(treated_post_sum - pred_post.sum() * n_treated)

    # p-value — conformal iid at h0=0 over the post block (Monte-Carlo)
    rng = np.random.default_rng(seed)
    resids0 = _conformal_resids_h0(M, trt, t0, 0.0, model, fixed_effects, lam)
    pvalue = _permute_pval(resids0, t0, ns, rng, stat_func)

    lo = hi = None
    ci_used = None
    if confidence_intervals:
        meth = method.lower()
        if stat_test.lower() != "total":
            meth = "conformal"                                  # jackknife+ is Total-only
        if meth == "jackknife+":
            lo, hi = _jackknife_plus_ci(M, trt, t0, alpha, model, fixed_effects, lam)
            ci_used = "jackknife+"
        else:
            lo, hi = _conformal_ci(M, trt, t0, alpha, grid_size, ns, rng,
                                   stat_func, model, fixed_effects, lam)
            ci_used = "conformal"
            if lo == np.inf or hi == -np.inf:                   # R's fallback
                lo, hi = _jackknife_plus_ci(M, trt, t0, alpha, model, fixed_effects, lam)
                ci_used = "jackknife+ (conformal failed)"

    donor_locs = [panel.locations[i] for i in range(len(trt)) if trt[i] == 0]
    return GeoLiftResult(
        att=att, percent_lift=float(100 * round(lift, 3)), pvalue=pvalue,
        incremental=incremental, lower_conf_int=lo, upper_conf_int=hi,
        l2_imbalance=float(_l2(M, trt, t0, means, weights)),
        scaled_l2_imbalance=float(scaled_l2), model=model, lam=lam,
        ci_method_used=ci_used, treatment_start=ts, treatment_end=te,
        weights=weights, donor_locations=donor_locs,
        y_hat=y_hat, att_series=att_series)


# --- imbalance helpers (augsynth l2 over the matched pre block) ------------- #
def _l2(M, trt, t0, means, weights):
    don = trt == 0
    Xd = M[:, :t0] - means[:, None]
    Z0 = Xd[don].T                                  # (t0, n_c)
    Z1 = (M[trt == 1, :t0] - means[trt == 1, None]).mean(axis=0)
    return np.sqrt(np.sum((Z0 @ weights - Z1) ** 2))


def _scaled_l2(M, trt, t0, means, weights):
    don = trt == 0
    Xd = M[:, :t0] - means[:, None]
    Z0 = Xd[don].T
    Z1 = (M[trt == 1, :t0] - means[trt == 1, None]).mean(axis=0)
    l2 = np.sqrt(np.sum((Z0 @ weights - Z1) ** 2))
    n_c = Z0.shape[1]
    uni = np.full(n_c, 1.0 / n_c)
    unif = np.sqrt(np.sum((Z0 @ uni - Z1) ** 2))
    return l2 / unif if unif > 0 else np.nan
