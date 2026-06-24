# =============================================================================
# run_experiments.R  —  reproducible source of every number in REPORT.md
#
# Writes machine-readable results to exploration/results/*.json.
# Statistics (ATT, lift, p-value) are made deterministic with a fixed SEED so
# anyone re-running gets identical values. Timings are wall-clock and will vary
# by machine/run; the RELATIVE pattern (and the % breakdown) is what matters.
#
# Run from the Geolift/ root:  Rscript exploration/scripts/run_experiments.R
# =============================================================================
suppressMessages(library(GeoLift))
suppressMessages(library(jsonlite))

SEED   <- 42
resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

data(GeoLift_Test)
d <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

LOC <- c("chicago", "portland"); T0 <- 91; T1 <- 105
fit_geolift <- function(...) GeoLift(Y_id = "Y", data = d, locations = LOC,
                                     treatment_start_time = T0, treatment_end_time = T1, ...)

# -----------------------------------------------------------------------------
# Experiment 1 — where the runtime goes (function-level profile, bucketed)
# -----------------------------------------------------------------------------
cat("[1/4] profiling a default GeoLift() call ...\n")
prof <- tempfile(fileext = ".out")
set.seed(SEED)
Rprof(prof, interval = 0.01, line.profiling = FALSE)
invisible(fit_geolift(ns = 1000))
Rprof(NULL)
self <- summaryRprof(prof)$by.self
buckets <- list(
  "Resampling RNG (permutation draws)" = c("sample", "sample.int"),
  "Input re-validation (stopifnot)"    = c("stopifnot"),
  "SCM model fit (quadratic program)"  = c("synth_qp", "osqpSolve", "osqpSetup",
                                           "fit_ridgeaug_inner", "fit_synth_formatted",
                                           "solve_osqp", "single_augsynth"),
  "Matrix / mean ops"                  = c("%*%", "mean", "apply", "t", "abs", "sum"))
total_self <- sum(self$self.time)
nm <- gsub('"', "", rownames(self))
used <- c()
prof_cat <- lapply(names(buckets), function(b) {
  idx <- nm %in% buckets[[b]]; used <<- c(used, nm[idx])
  list(category = b, self_time_s = round(sum(self$self.time[idx]), 2),
       pct = round(100 * sum(self$self.time[idx]) / total_self, 1))
})
other <- total_self - sum(self$self.time[nm %in% used])
prof_cat[[length(prof_cat) + 1]] <- list(category = "Other (glue/dispatch)",
  self_time_s = round(other, 2), pct = round(100 * other / total_self, 1))
write_json(list(total_sampled_s = round(total_self, 2), categories = prof_cat),
           file.path(resdir, "runtime_profile.json"), auto_unbox = TRUE, pretty = TRUE)

# -----------------------------------------------------------------------------
# Experiment 2 — Tier 1: cutting `ns` (the resample count) for the point p-value
# Shows time drops ~linearly while ATT / Lift are unchanged.
# -----------------------------------------------------------------------------
cat("[2/4] ns sweep (point p-value) ...\n")
ns_grid <- c(1000, 500, 200, 100)
sweep <- do.call(rbind, lapply(ns_grid, function(nsv) {
  set.seed(SEED)
  t0 <- proc.time()[["elapsed"]]
  g  <- fit_geolift(ns = nsv)
  dt <- proc.time()[["elapsed"]] - t0
  data.frame(ns = nsv, time_s = round(dt, 2),
             ATT = round(g$inference$ATT, 3),
             lift_pct = round(g$inference$Perc.Lift, 3),
             pvalue = round(g$inference$pvalue, 4),
             significant_at_0.10 = g$inference$pvalue < 0.10)
}))
write_json(sweep, file.path(resdir, "tier1_ns_sweep.json"),
           dataframe = "rows", auto_unbox = TRUE, pretty = TRUE)

# -----------------------------------------------------------------------------
# Experiment 3 — Tier 1: CI via conformal grid vs. jackknife+ directly
# (the conformal grid often fails on this data and falls back to jackknife+,
#  so calling jackknife+ directly skips wasted work)
# -----------------------------------------------------------------------------
cat("[3/4] CI method comparison ...\n")
ci_cfg <- list(
  list(label = "conformal grid=250, ns=1000", method = "conformal", grid_size = 250, ns = 1000),
  list(label = "jackknife+ (direct)",          method = "jackknife+", grid_size = 250, ns = 1000))
ci_res <- do.call(rbind, lapply(ci_cfg, function(cfg) {
  set.seed(SEED)
  t0 <- proc.time()[["elapsed"]]
  g  <- tryCatch(fit_geolift(ConfidenceIntervals = TRUE, method = cfg$method,
                             grid_size = cfg$grid_size, ns = cfg$ns),
                 error = function(e) NULL)
  dt <- proc.time()[["elapsed"]] - t0
  if (is.null(g)) return(data.frame(method = cfg$label, time_s = NA, lower = NA, upper = NA))
  data.frame(method = cfg$label, time_s = round(dt, 2),
             lower = round(g$lower, 1), upper = round(g$upper, 1))
}))
write_json(ci_res, file.path(resdir, "tier1_ci_methods.json"),
           dataframe = "rows", auto_unbox = TRUE, pretty = TRUE)

# -----------------------------------------------------------------------------
# Experiment 4 — same UNMODIFIED GeoLift(), two parameter configs: the exact Meta
# default (ns=1000) vs our tuned preset (ns=200). The "improved version" is a set
# of PARAMETERS, not a forked function. Confirms identical estimate, faster runtime.
# -----------------------------------------------------------------------------
cat("[4/4] default params vs fast preset (same GeoLift engine) ...\n")
source("geolift_r_modified/fast_presets.R")
time_it <- function(expr) { t0 <- proc.time()[["elapsed"]]; g <- force(expr)
  list(g = g, dt = proc.time()[["elapsed"]] - t0) }

# Both calls hit the SAME GeoLift::GeoLift; only the parameter list differs.
set.seed(SEED); orig <- time_it(do.call(fit_geolift, GEOLIFT_DEFAULT))
set.seed(SEED); fast <- time_it(do.call(fit_geolift, GEOLIFT_FAST))
compare <- data.frame(
  version  = c("default params (ns=1000)", "fast preset (ns=200)"),
  time_s   = c(round(orig$dt, 2), round(fast$dt, 2)),
  ATT      = c(round(orig$g$inference$ATT, 3),  round(fast$g$inference$ATT, 3)),
  lift_pct = c(round(orig$g$inference$Perc.Lift, 3), round(fast$g$inference$Perc.Lift, 3)),
  pvalue   = c(round(orig$g$inference$pvalue, 4), round(fast$g$inference$pvalue, 4)))
compare$speedup <- round(compare$time_s[1] / compare$time_s, 2)
write_json(compare, file.path(resdir, "tier1_modified_compare.json"),
           dataframe = "rows", auto_unbox = TRUE, pretty = TRUE)

# Provenance so the report can state how numbers were produced.
write_json(list(seed = SEED, r_version = R.version.string,
                example = "GeoLift Walkthrough: chicago+portland, t=91..105",
                augsynth_commit = "06415b4 (pre-PR#88)",
                note = "Statistics are deterministic given SEED; timings vary by machine/run."),
           file.path(resdir, "provenance.json"), auto_unbox = TRUE, pretty = TRUE)

cat("Done. Results in", resdir, "\n")
