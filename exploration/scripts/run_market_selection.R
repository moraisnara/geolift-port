# =============================================================================
# run_market_selection.R  —  power/market-selection path: default vs fast PARAMS
#
# Demonstrates the principle: the SAME unmodified GeoLift::GeoLiftMarketSelection
# is run repeatedly, differing ONLY in the `ns` parameter (1000 = Meta default,
# down to 100). We then check that the selected markets and the power curves do
# NOT move, while wall-clock drops. Small, reproducible config (sequential +
# fixed seed). Scale-up is a parameter change (bigger N / ES) — see the scaled
# variant. Shared sweep/scoring logic lives in ms_sweep_common.R.
#
# Run from Geolift/ root:  Rscript exploration/scripts/run_market_selection.R
# =============================================================================
suppressMessages(library(jsonlite))
source("geolift_r_modified/fast_presets.R")
source("exploration/scripts/ms_sweep_common.R")

SEED <- 42
resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

# ---- small representative config -------------------------------------------
# Keep a real donor pool but limit the search so it runs in ~1-2 min sequential.
data(GeoLift_Test)
keep_locs <- sort(unique(GeoLift_Test$location))[1:20]      # deterministic subset
sub <- GeoLift_Test[GeoLift_Test$location %in% keep_locs, ]
d <- GeoDataRead(data = sub, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

CFG <- list(
  data = d, treatment_periods = c(15), N = c(2),
  effect_size = seq(-0.15, 0.15, 0.05),
  Y_id = "Y", location_id = "location", time_id = "time",
  parallel = FALSE,                 # sequential => set.seed reproducible
  print = FALSE
)

sweep <- run_ns_sweep(CFG, c(1000, 500, 250, 100), SEED)
rows <- sweep$rows

summary <- list(
  setup = list(seed = SEED, locations_in_pool = length(keep_locs),
               treatment_periods = CFG$treatment_periods, N = CFG$N,
               effect_sizes = CFG$effect_size, parallel = FALSE,
               n_simulations_per_run = rows[[1]]$n_cells,
               reference = "ns=1000 (Meta default); lower ns scored against it"),
  sweep = rows)
write_json(summary, file.path(resdir, "market_selection_compare.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8)

cat("\n--- market selection: ns sweep vs ns=1000 reference ---\n")
cat(sprintf("  sims/run: %d\n", rows[[1]]$n_cells))
print_sweep(rows)
cat("\nWrote", file.path(resdir, "market_selection_compare.json"), "\n")
