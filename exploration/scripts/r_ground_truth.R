# =============================================================================
# GeoLift R ground-truth + profiling
# Reproduces the documented Walkthrough examples, captures every numeric result
# as a benchmark for the Python port, and times each step to locate bottlenecks.
# =============================================================================
suppressMessages(library(GeoLift))

outdir <- "exploration/ground_truth"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)   # make the permutation p-values stored in the benchmark reproducible

have_json <- requireNamespace("jsonlite", quietly = TRUE)
timings <- list()
timed <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- force(expr)
  dt <- proc.time()[["elapsed"]] - t0
  timings[[label]] <<- dt
  cat(sprintf("[time] %-28s %8.3f s\n", label, dt))
  val
}

# ---- 1. Data read (data wrangling path) -------------------------------------
data(GeoLift_PreTest)
data(GeoLift_Test)

GeoTestData_PreTest <- timed("GeoDataRead_PreTest", GeoDataRead(
  data = GeoLift_PreTest, date_id = "date", location_id = "location",
  Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE))

GeoTestData_Test <- timed("GeoDataRead_Test", GeoDataRead(
  data = GeoLift_Test, date_id = "date", location_id = "location",
  Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE))

cat(sprintf("PreTest: %d rows, %d locations, %d periods\n",
            nrow(GeoTestData_PreTest), length(unique(GeoTestData_PreTest$location)),
            length(unique(GeoTestData_PreTest$time))))

# ---- 2. Core GeoLift inference (the headline documented example) ------------
GeoTest <- timed("GeoLift_None", GeoLift(
  Y_id = "Y", data = GeoTestData_Test,
  locations = c("chicago", "portland"),
  treatment_start_time = 91, treatment_end_time = 105))

GeoTestBest <- timed("GeoLift_best", GeoLift(
  Y_id = "Y", data = GeoTestData_Test,
  locations = c("chicago", "portland"),
  treatment_start_time = 91, treatment_end_time = 105, model = "best"))

# ---- 3. Extract scalar results defensively ----------------------------------
grab <- function(obj) {
  s <- summary(obj)
  pick <- function(x, nm) if (!is.null(x[[nm]])) x[[nm]] else NA
  list(
    incremental = pick(obj, "incremental"),
    att_estimator = pick(obj, "att_estimator"),
    percent_lift = pick(obj, "ATT_se"),  # placeholder, overwritten below
    obj_names = names(obj),
    summary_names = names(s)
  )
}
# Dump full structure so we can see exact field names, then extract precisely.
sink(file.path(outdir, "GeoTest_str.txt")); str(GeoTest, max.level = 2); sink()

# Known GeoLift object fields (from post_test_analysis.R)
extract <- function(obj) {
  f <- function(nm) { v <- tryCatch(obj[[nm]], error = function(e) NA); if (is.null(v)) NA else v }
  list(
    incremental            = f("incremental"),
    att_estimator          = f("att_estimator"),
    percent_lift           = f("percent_lift"),
    pvalue                 = f("p_value"),
    lower_bound            = f("lower_bound"),
    upper_bound            = f("upper_bound"),
    l2_imbalance           = f("l2_imbalance"),
    scaled_l2_imbalance    = f("scaled_l2_imbalance"),
    treatment_start_time   = f("TreatmentStart"),
    ci_widths              = f("ci_width")
  )
}
res_none <- extract(GeoTest)
res_best <- extract(GeoTestBest)
cat("\n--- GeoLift(None) extracted ---\n"); print(res_none)

# Synthetic-control weights (the modelling output we must match in Python)
w <- tryCatch(GetWeights(
  data = GeoTestData_Test, locations = c("chicago", "portland"),
  treatment_start_time = 91, treatment_end_time = 105), error = function(e) NULL)
if (!is.null(w)) write.csv(w, file.path(outdir, "weights_None.csv"), row.names = FALSE)

# ---- 4. Save data matrices Python will load (exact same inputs) -------------
write.csv(GeoTestData_PreTest, file.path(outdir, "GeoTestData_PreTest.csv"), row.names = FALSE)
write.csv(GeoTestData_Test,    file.path(outdir, "GeoTestData_Test.csv"),    row.names = FALSE)

# ---- 5. Profile a small power run (likely the real bottleneck) --------------
pwr <- timed("GeoLiftPower_small", GeoLiftPower(
  data = GeoTestData_PreTest,
  locations = c("chicago", "portland"),
  effect_size = seq(0, 0.15, 0.05),
  treatment_periods = c(15),
  cpic = 7.5, side_of_test = "one_sided", parallel = FALSE))

# ---- 6. Persist everything --------------------------------------------------
saveRDS(list(none = GeoTest, best = GeoTestBest), file.path(outdir, "geolift_objects.rds"))
bundle <- list(
  results_None = res_none, results_best = res_best,
  timings_seconds = timings,
  doc_expected = list(percent_lift = 5.4, incremental = 4667,
                      att_per_day = 155.556, pvalue_approx = 0.01)
)
if (have_json) {
  writeLines(jsonlite::toJSON(bundle, auto_unbox = TRUE, pretty = TRUE, na = "null"),
             file.path(outdir, "ground_truth.json"))
} else {
  dput(bundle, file.path(outdir, "ground_truth.dput"))
}
cat("\n=== DONE. Wrote ground truth to", outdir, "===\n")
