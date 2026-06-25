# =============================================================================
# Dump the DETERMINISTIC conformal residual vector that feeds GeoLift's headline
# p-value, so the Python port can be checked against it exactly (instead of only
# comparing Monte-Carlo p-values).
#
# The headline conformal p-value (augsynth conformal_inf -> compute_permute_pval)
# refits the SCM over the FULL series with H0: ATT==0 subtracted (a no-op at 0),
# then permutes the resulting ATT-path residuals. That residual vector is fully
# deterministic. If R and Python produce the same vector, their permutation
# p-values converge to the same limit -- equivalence is then a fact about the
# residuals, not an argument about RNG.
#
#   Rscript exploration/scripts/dump_conformal_resids.R
# Writes exploration/results/conformal_resids.json
# =============================================================================
suppressMessages({library(GeoLift); library(augsynth); library(jsonlite)})
set.seed(42)

data(GeoLift_Test)
geo <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                   Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)
LOCS <- c("chicago", "portland"); TS <- 91; TE <- 105

stat_total <- function(x) sum(abs(x))   # GeoLift stat_test = "Total" (two-sided)

dump_one <- function(model) {
  set.seed(42)
  g <- GeoLift(Y_id = "Y", data = geo, locations = LOCS,
               treatment_start_time = TS, treatment_end_time = TE,
               alpha = 0.1, model = model, ConfidenceIntervals = FALSE)
  ascm <- g$results                       # the augsynth fit object
  wide <- ascm$data
  post_length <- ncol(wide$y)
  # replicate conformal_inf's headline-pval data: X <- [X | y], y <- dummy
  nw <- wide
  nw$X <- cbind(wide$X, wide$y)
  nw$y <- matrix(1, nrow = nrow(wide$X), ncol = 1)
  out <- augsynth:::compute_permute_test_stats(nw, ascm, 0, post_length,
                                               "iid", 1, 1, stat_total)
  t0 <- ncol(nw$X) - post_length
  tpost <- t0 + post_length
  list(model = model,
       resids = as.numeric(out$resids),                 # full-series-refit ATT path
       obs_stat = stat_total(out$resids[(t0 + 1):tpost]),
       t0 = t0, post_length = post_length)
}

out <- list(none = dump_one("none"), ridge = dump_one("ridge"),
            note = "Deterministic full-series-refit ATT residuals; H0 ATT=0.")
resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = 12),
           file.path(resdir, "conformal_resids.json"))
cat("Wrote conformal_resids.json; none obs_stat =", out$none$obs_stat,
    " ridge obs_stat =", out$ridge$obs_stat, "\n")
