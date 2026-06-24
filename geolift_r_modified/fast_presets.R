# =============================================================================
# fast_presets.R  —  faster PARAMETER configurations for the UNMODIFIED GeoLift
#
# Design principle: our speedups should be *parameters you pass to the original
# Meta functions*, not forked functions. So you can run the EXACT default Meta
# package, or our tuned configuration, by changing arguments only — same engine.
#
#   library(GeoLift)
#   source("geolift_r_modified/fast_presets.R")
#
#   # exact default Meta:
#   GeoLift(...)                                   # or do.call(GeoLift, c(args, GEOLIFT_DEFAULT))
#   # our tuned config (same function, different params):
#   do.call(GeoLift, c(args, GEOLIFT_FAST))
#
# WHY these parameters (verified — see exploration/scripts/ and REPORT.md):
#   * Point estimates (ATT, Lift, ScaledL2Imbalance, selected markets' ranking
#     basis) do NOT depend on `ns`. `ns` only sets the resample count for the
#     conformal p-value, whose dominant cost is exactly those `ns` permutations.
#   * In the power / market-selection path that p-value runs ONCE PER SIMULATION
#     (#effect_sizes x #durations x #combinations x #N), so lowering `ns` is a
#     *multiplied* saving. We validate the selected markets / power curves still
#     match the default (exploration/scripts/run_market_selection.R).
#   * For single-call CIs, `method = "jackknife+"` skips the conformal grid that
#     often fails on geo data and falls back to jackknife+ anyway.
#
# These are the EXPOSED-PARAMETER improvements. Structural improvements that are
# NOT yet parameters of the Meta functions (reuse the es-invariant fit across
# effect sizes; drop duplicated AggYperLoc / per-sim normalization) would need a
# code change to run_simulations()/pvalueCalc() and are tracked separately as a
# future modified engine — kept out of here on purpose.
# =============================================================================

# ---- single post-test call: GeoLift() ---------------------------------------
GEOLIFT_DEFAULT <- list(ns = 1000, method = "conformal")       # exact Meta defaults
GEOLIFT_FAST    <- list(ns = 200)                              # tuned: fewer resamples
# For confidence intervals specifically, add method = "jackknife+":
GEOLIFT_FAST_CI <- list(ns = 200, ConfidenceIntervals = TRUE, method = "jackknife+")

# ---- power / market selection: GeoLiftMarketSelection(), GeoLiftPower(), ... --
MARKETSEL_DEFAULT <- list(ns = 1000)                           # exact Meta default
MARKETSEL_FAST    <- list(ns = 100)                            # tuned: fewer resamples/sim

# Convenience: run any GeoLift function with a preset merged over your arguments.
# Caller-supplied args win over the preset (so you can still override ns, etc.).
run_with_preset <- function(fn, preset = list(), ...) {
  user <- list(...)
  merged <- modifyList(preset, user)          # user args override preset
  do.call(fn, merged)
}
