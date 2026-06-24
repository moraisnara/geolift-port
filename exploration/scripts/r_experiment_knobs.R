# Quantify how inference knobs (ns, grid_size, ConfidenceIntervals) drive runtime,
# and whether the headline numbers stay stable when we cut them.
suppressMessages(library(GeoLift))
set.seed(1)
data(GeoLift_Test)
d <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

run <- function(label, ...) {
  t0 <- proc.time()[["elapsed"]]
  g  <- GeoLift(Y_id = "Y", data = d, locations = c("chicago", "portland"),
                treatment_start_time = 91, treatment_end_time = 105, ...)
  dt <- proc.time()[["elapsed"]] - t0
  inf <- g$inference
  cat(sprintf("%-34s  %6.2fs | ATT=%.2f  lift=%.2f%%  p=%.4f  CI=(%s, %s)\n",
      label, dt, inf$ATT, inf$Perc.Lift, inf$pvalue,
      ifelse(is.na(g$lower), "NA", round(g$lower,1)),
      ifelse(is.na(g$upper), "NA", round(g$upper,1))))
  invisible(dt)
}

cat("== POINT P-VALUE ONLY (ConfidenceIntervals = FALSE, default) ==\n")
run("default  ns=1000",            ns = 1000)
run("modified ns=500",             ns = 500)
run("modified ns=200",             ns = 200)
run("modified ns=100",             ns = 100)

cat("\n== WITH CONFIDENCE INTERVALS (grid_size x ns inversion) ==\n")
run("default  CI grid=250 ns=1000", ConfidenceIntervals = TRUE, grid_size = 250, ns = 1000)
run("modified CI grid=100 ns=500",  ConfidenceIntervals = TRUE, grid_size = 100, ns = 500)
run("modified CI grid=40  ns=250",  ConfidenceIntervals = TRUE, grid_size = 40,  ns = 250)
