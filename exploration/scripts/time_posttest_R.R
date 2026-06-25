# =============================================================================
# Time R's GeoLift() POST-TEST on the Walkthrough, to pair against the Python
# geolift_fast.geolift() timing. Two operations:
#   point_pval        : GeoLift(..., ConfidenceIntervals = FALSE)  (ATT + p-value)
#   with_ci_conformal : GeoLift(..., ConfidenceIntervals = TRUE)   (+ CI grid)
# Single-threaded; best-of-N to mirror the Python timer.
#
#   Rscript exploration/scripts/time_posttest_R.R
# Writes exploration/results/posttest_R_time.json  (a TIMING artifact).
# =============================================================================
suppressMessages({library(GeoLift); library(jsonlite)})

data(GeoLift_Test)
geo <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                   Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)
LOCS <- c("chicago", "portland"); TS <- 91; TE <- 105

best_time <- function(fn, reps = 5) {
  best <- Inf
  for (i in seq_len(reps)) {
    set.seed(42)
    t <- system.time(fn())[["elapsed"]]
    best <- min(best, t)
  }
  best
}

t_point <- best_time(function() {
  GeoLift(Y_id = "Y", data = geo, locations = LOCS,
          treatment_start_time = TS, treatment_end_time = TE,
          alpha = 0.1, model = "none", ConfidenceIntervals = FALSE)
})

t_ci <- best_time(function() {
  GeoLift(Y_id = "Y", data = geo, locations = LOCS,
          treatment_start_time = TS, treatment_end_time = TE,
          alpha = 0.1, model = "none", ConfidenceIntervals = TRUE,
          method = "conformal")
}, reps = 3)

out <- list(point_pval_s = round(t_point, 4),
            with_ci_conformal_s = round(t_ci, 4),
            ns = 1000, grid_size = 250, reps = "best-of-5 / best-of-3",
            note = "Single-threaded R GeoLift() wall time, best-of-N.")
resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)
writeLines(toJSON(out, auto_unbox = TRUE, pretty = TRUE, digits = 6),
           file.path(resdir, "posttest_R_time.json"))
cat(sprintf("R post-test: point+pval=%.3fs  with-CI(conformal)=%.3fs\n", t_point, t_ci))
