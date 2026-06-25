# =============================================================================
# run_market_selection_scaled.R  —  the same ns sweep as run_market_selection.R,
# but on the LARGER realistic panel (exploration/data/realistic_panel.csv) to test
# whether the small-config finding holds at scale: does the ns speedup grow with
# data size, do point estimates stay invariant, does the ranking stay stable?
#
# Same principle: one UNMODIFIED GeoLift::GeoLiftMarketSelection, only `ns` differs.
# Sequential + fixed seed => reproducible and isolates the ns effect (parallelism
# is an orthogonal, ~linear multiplier on top). Shared sweep/scoring logic lives
# in ms_sweep_common.R.
#
# Run from Geolift/ root:  Rscript exploration/scripts/run_market_selection_scaled.R
# =============================================================================
suppressMessages(library(jsonlite))
source("geolift_r_modified/fast_presets.R")
source("exploration/scripts/ms_sweep_common.R")

SEED <- 42
resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

info <- fromJSON("exploration/data/realistic_panel_info.json")
panel <- read.csv("exploration/data/realistic_panel.csv", stringsAsFactors = FALSE)
d <- GeoDataRead(data = panel, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

CFG <- list(
  data = d, treatment_periods = c(14), N = c(3),
  effect_size = seq(-0.1, 0.1, 0.05),
  Y_id = "Y", location_id = "location", time_id = "time",
  parallel = FALSE, print = FALSE)

sweep <- run_ns_sweep(CFG, c(1000, 500, 250), SEED)
rows <- sweep$rows

summary <- list(
  dataset = list(n_locations = info$n_locations, n_periods = info$n_periods,
                 n_categories = info$n_categories, avg_pairwise_correlation = info$avg_pairwise_correlation,
                 compared_to = info$compared_to_GeoLift_Test),
  setup = list(seed = SEED, treatment_periods = CFG$treatment_periods, N = CFG$N,
               effect_sizes = CFG$effect_size, parallel = FALSE,
               n_simulations_per_run = rows[[1]]$n_cells,
               reference = "ns=1000 (Meta default); lower ns scored against it"),
  sweep = rows)
write_json(summary, file.path(resdir, "market_selection_scaled_compare.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8)

cat("\n--- SCALED market selection: ns sweep vs ns=1000 reference ---\n")
cat(sprintf("  dataset: %d locations x %d periods, %d sims/run\n",
            info$n_locations, info$n_periods, rows[[1]]$n_cells))
print_sweep(rows)
cat("\nWrote", file.path(resdir, "market_selection_scaled_compare.json"), "\n")
