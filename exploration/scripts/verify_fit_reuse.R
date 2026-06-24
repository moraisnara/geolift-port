# =============================================================================
# verify_fit_reuse.R  —  test the "reuse work across effect sizes" hypothesis
# for the power/market-selection inner loop (pvalueCalc).
#
# Question: for a FIXED (market combination, duration), as we sweep effect size
# `es`, what is reusable?
#   H1: the SCM weights / counterfactual are es-invariant  (es only inflates the
#       treated unit's POST values; SCM fits on the PRE period).
#   H2: the permutation NULL distribution is es-invariant  (=> could compute once
#       and reuse across es for a big speedup).
#
# We replicate pvalueCalc's exact construction and, for each es, record:
#   * the SCM weights and the pre-period fit (to test H1)
#   * the att estimator and the conformal p-value with a FIXED seed (to test H2)
#
# Writes exploration/results/fit_reuse_verify.json (numbers, not hardcoded).
# Run from Geolift/ root:  Rscript exploration/scripts/verify_fit_reuse.R
# =============================================================================
suppressMessages(library(GeoLift))
suppressMessages(library(augsynth))
suppressMessages(library(jsonlite))

SEED <- 42
resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

data(GeoLift_Test)
d <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

LOCS <- c("chicago", "portland")
tp   <- 15                 # treatment duration (periods)
max_time <- max(d$time)
sim <- 1
es_grid <- c(0.00, 0.05, 0.10, 0.15)

fn_treatment <- getFromNamespace("fn_treatment", "GeoLift")
treatment_start_time <- max_time - tp - sim + 2
treatment_end_time   <- treatment_start_time + tp - 1

# ---- replicate pvalueCalc's fit + p-value for a single es -------------------
one_es <- function(es) {
  data_aux <- fn_treatment(d, locations = LOCS, treatment_start_time, treatment_end_time)
  data_aux$Y_inc <- data_aux$Y
  data_aux$Y_inc[data_aux$D == 1] <- data_aux$Y_inc[data_aux$D == 1] * (1 + es)

  ascm <- augsynth::augsynth(Y_inc ~ D, unit = location, time = time, data = data_aux,
                             t_int = treatment_start_time, progfunc = "none",
                             scm = TRUE, fixedeff = TRUE)

  # counterfactual on the post window (es-invariant under H1)
  cf_post <- as.numeric(predict(ascm)[treatment_start_time:treatment_end_time])
  ave_trt <- sum(ascm$data$y[which(ascm$data$trt == 1), ]) / sum(ascm$data$trt == 1)
  att <- (ave_trt - sum(cf_post)) / tp

  # conformal p-value, exactly as pvalueCalc builds it, with a FIXED seed
  wide <- ascm$data
  nw <- wide; nw$X <- cbind(wide$X, wide$y); nw$y <- matrix(1, nrow = nrow(wide$X), ncol = 1)
  set.seed(SEED)
  pv <- augsynth:::compute_permute_pval(wide_data = nw, ascm = ascm, h0 = 0,
                                        post_length = ncol(wide$y), type = "iid",
                                        q = 1, ns = 1000, stat_func = NULL)

  list(es = es,
       weights = round(as.numeric(ascm$weights), 8),
       counterfactual_post = round(cf_post, 6),
       att = round(att, 4),
       pvalue = round(pv, 4))
}

runs <- lapply(es_grid, one_es)

# ---- compare against the es=0 baseline -------------------------------------
base <- runs[[1]]
cmp <- lapply(runs, function(r) {
  list(es = r$es,
       weights_max_abs_diff_vs_es0 = signif(max(abs(r$weights - base$weights)), 3),
       cf_max_abs_diff_vs_es0      = signif(max(abs(r$counterfactual_post - base$counterfactual_post)), 3),
       att = r$att,
       pvalue = r$pvalue)
})

out <- list(
  setup = list(seed = SEED, locations = paste(LOCS, collapse = "+"), duration_tp = tp,
               ns = 1000, es_grid = es_grid),
  hypothesis = list(
    H1_weights_es_invariant = "expect weights_max_abs_diff ~ 0 across es",
    H2_null_es_invariant    = "if pvalue changes with es, the null is NOT reusable across es"),
  per_es = cmp)
write_json(out, file.path(resdir, "fit_reuse_verify.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8)

cat("\n--- fit-reuse verification (seed", SEED, ") ---\n")
cat(sprintf("%-6s %22s %18s %10s %9s\n", "es", "weightΔ vs es0", "cfΔ vs es0", "att", "pvalue"))
for (r in cmp) cat(sprintf("%-6.2f %22g %18g %10s %9s\n",
                           r$es, r$weights_max_abs_diff_vs_es0, r$cf_max_abs_diff_vs_es0,
                           format(r$att), format(r$pvalue)))
cat("\nWrote", file.path(resdir, "fit_reuse_verify.json"), "\n")
