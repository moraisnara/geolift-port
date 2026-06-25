# =============================================================================
# build_city_example.R — run GeoLift's market-selection example on the CANONICAL
# GeoLift_Test data (real US city names) and dump the artifacts the Python
# example notebook compares against, so we can show R and Python produce the
# SAME results on the SAME example.
#
# Writes:
#   exploration/data/geolift_test_panel.csv          raw long panel (location,date,Y)
#   exploration/results/citylift_R_powercurves.csv   R $PowerCurves (the combos R chose)
#   exploration/results/citylift_R_bestmarkets.csv   R $BestMarkets (ranked, with MDE)
#
# Run from Geolift/ root:  Rscript exploration/scripts/build_city_example.R
# =============================================================================
suppressMessages(library(GeoLift))
suppressMessages(library(dplyr))

SEED <- 42
datadir <- "exploration/data"
resdir  <- "exploration/results"
dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
dir.create(resdir,  showWarnings = FALSE, recursive = TRUE)

# ---- canonical GeoLift example data (40 US cities x 105 days) ---------------
data(GeoLift_Test)
write.csv(GeoLift_Test, file.path(datadir, "geolift_test_panel.csv"), row.names = FALSE)

d <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

# ---- market selection: same knobs the Python notebook uses ------------------
set.seed(SEED)
out <- suppressMessages(GeoLift::GeoLiftMarketSelection(
  data = d,
  treatment_periods = c(14),
  N = c(2),
  Y_id = "Y", location_id = "location", time_id = "time",
  effect_size = seq(-0.10, 0.10, 0.05),
  ns = 1000,
  alpha = 0.10,
  fixed_effects = TRUE,
  Correlations = TRUE,
  side_of_test = "two_sided",
  parallel = FALSE,
  print = FALSE
))

pc <- out$PowerCurves %>% dplyr::arrange(location, duration, EffectSize)
bm <- out$BestMarkets

write.csv(pc, file.path(resdir, "citylift_R_powercurves.csv"), row.names = FALSE)
write.csv(bm, file.path(resdir, "citylift_R_bestmarkets.csv"), row.names = FALSE)

cat("\n--- R market selection on GeoLift_Test (city names) ---\n")
cat(sprintf("  panel        : %d locations x %d periods\n",
            length(unique(GeoLift_Test$location)), length(unique(GeoLift_Test$time))))
cat(sprintf("  combos tested: %d   (cells: %d)\n",
            length(unique(pc$location)), nrow(pc)))
cat("  R BestMarkets (top 5):\n")
print(utils::head(bm, 5))
cat("\nWrote:\n  ", file.path(datadir, "geolift_test_panel.csv"),
    "\n  ", file.path(resdir, "citylift_R_powercurves.csv"),
    "\n  ", file.path(resdir, "citylift_R_bestmarkets.csv"), "\n")
