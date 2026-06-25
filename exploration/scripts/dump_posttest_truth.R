# =============================================================================
# Post-test ground truth for the Python geolift() port.
# Runs the documented Walkthrough example (GeoLift_Test, chicago+portland,
# t=91..105) through R's GeoLift() post-test with confidence intervals, for
# model = "none" and "ridge", and dumps every numeric result the Python port
# must reproduce: ATT, Lift, Incremental, p-value, conformal CI, jackknife+ CI,
# weights, counterfactual, plus the ridge lambda.
#
# Statistics are deterministic given SEED *except* the conformal p-value / CI,
# which carry Monte-Carlo noise (independent of Python's RNG). Run:
#   Rscript exploration/scripts/dump_posttest_truth.R
# Writes exploration/results/posttest_truth.json  (a stats artifact, not timing).
# =============================================================================
suppressMessages({library(GeoLift); library(jsonlite)})
set.seed(42)

data(GeoLift_Test)
geo <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                   Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)
LOCS <- c("chicago", "portland"); TS <- 91; TE <- 105; ALPHA <- 0.1

run_one <- function(model, method) {
  set.seed(42)
  g <- GeoLift(Y_id = "Y", data = geo, locations = LOCS,
               treatment_start_time = TS, treatment_end_time = TE,
               alpha = ALPHA, model = model,
               ConfidenceIntervals = TRUE, method = method)
  inf <- g$inference
  w <- g$df_weights
  list(
    model = model, method_requested = method,
    ATT = g$ATT_se,                         # NA for conformal (placeholder)
    att_estimator = inf$ATT,
    percent_lift = inf$Perc.Lift,
    pvalue = inf$pvalue,
    lower_conf_int = inf$Lower.Conf.Int,
    upper_conf_int = inf$Upper.Conf.Int,
    lower_bound = g$lower_bound,            # CI x periods x n_treated (aggregate)
    upper_bound = g$upper_bound,
    incremental = g$incremental,
    l2_imbalance = g$summary$l2_imbalance,
    scaled_l2_imbalance = g$summary$scaled_l2_imbalance,
    lambda = if (!is.null(g$results$lambda)) g$results$lambda else NA,
    weights = setNames(as.list(w$weight), w$location),
    y_hat = as.numeric(g$y_hat),
    att_series = as.numeric(g$ATT)
  )
}

out <- list(
  meta = list(
    example = "GeoLift Walkthrough: chicago+portland, t=91..105",
    seed = 42, alpha = ALPHA, augsynth_commit = "06415b4 (pre-PR#88)",
    note = paste("Deterministic stats match to solver tolerance; conformal",
                 "p-value and conformal CI carry Monte-Carlo noise.")),
  none_conformal   = run_one("none",  "conformal"),
  none_jackknife   = run_one("none",  "jackknife+"),
  ridge_conformal  = run_one("ridge", "conformal"),
  ridge_jackknife  = run_one("ridge", "jackknife+")
)

resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = 12, na = "null"),
           file.path(resdir, "posttest_truth.json"))
cat("Wrote", file.path(resdir, "posttest_truth.json"), "\n")
for (nm in c("none_conformal","none_jackknife","ridge_conformal","ridge_jackknife")) {
  r <- out[[nm]]
  cat(sprintf("%-16s ATT=%.4f lift=%.2f p=%.3f CI=[%s, %s] lambda=%s\n",
              nm, r$att_estimator, r$percent_lift, r$pvalue,
              format(r$lower_conf_int), format(r$upper_conf_int), format(r$lambda)))
}
