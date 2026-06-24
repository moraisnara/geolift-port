# =============================================================================
# run_market_selection.R  —  power/market-selection path: default vs fast PARAMS
#
# Demonstrates the principle: the SAME unmodified GeoLift::GeoLiftMarketSelection
# is run twice, differing ONLY in the `ns` parameter:
#     * default (Meta): ns = 1000
#     * fast preset    : ns = 100
# We then check that the selected markets and the power curves do NOT move, while
# wall-clock drops. Small, reproducible config (sequential + fixed seed) so anyone
# re-running gets the same numbers. Scale-up is a parameter change (bigger N / ES).
#
# Run from Geolift/ root:  Rscript exploration/scripts/run_market_selection.R
# =============================================================================
suppressMessages(library(GeoLift))
suppressMessages(library(jsonlite))
suppressMessages(library(dplyr))
source("geolift_r_modified/fast_presets.R")

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

run_ms <- function(ns_val) {
  set.seed(SEED)
  t0 <- proc.time()[["elapsed"]]
  out <- suppressMessages(do.call(GeoLift::GeoLiftMarketSelection, c(CFG, list(ns = ns_val))))
  list(out = out, dt = proc.time()[["elapsed"]] - t0)
}

# Sweep ns. ns=1000 (Meta default) is the REFERENCE; every lower ns is scored
# against it, so we can see where the selected markets / power curves stabilise.
ns_grid <- c(1000, 500, 250, 100)
runs <- list()
for (i in seq_along(ns_grid)) {
  cat(sprintf("[%d/%d] ns=%d ...\n", i, length(ns_grid), ns_grid[i]))
  runs[[i]] <- run_ms(ns_grid[i])
}
ref <- runs[[1]]                                  # ns = 1000
pc_ref <- ref$out$PowerCurves %>% dplyr::arrange(location, duration, EffectSize)
bm_ref <- ref$out$BestMarkets

score_vs_ref <- function(r) {
  pc <- r$out$PowerCurves %>% dplyr::arrange(location, duration, EffectSize)
  j <- dplyr::inner_join(pc_ref, pc, by = c("location", "duration", "EffectSize"),
                         suffix = c("_ref", "_x"))
  bm <- r$out$BestMarkets
  top_k <- min(5, nrow(bm_ref), nrow(bm))
  top_overlap <- if (top_k > 0)
    length(intersect(head(bm_ref$location, top_k), head(bm$location, top_k))) / top_k else NA
  mde_join <- dplyr::inner_join(
    bm_ref[, c("location", "EffectSize")], bm[, c("location", "EffectSize")],
    by = "location", suffix = c("_ref", "_x"))
  mde_same <- if (nrow(mde_join) > 0) mean(mde_join$EffectSize_ref == mde_join$EffectSize_x) else NA
  list(n_cells = nrow(j),
       power_significance_agreement = round(mean(j$power_ref == j$power_x), 4),
       avg_att_max_abs_diff = signif(max(abs(j$AvgATT_ref - j$AvgATT_x)), 3),
       scaled_l2_imbalance_max_abs_diff = signif(max(abs(j$AvgScaledL2Imbalance_ref - j$AvgScaledL2Imbalance_x)), 3),
       bestmarkets_top5_overlap = round(top_overlap, 3),
       bestmarkets_selected_mde_agreement = round(mde_same, 4),
       n_markets = nrow(bm))
}

rows <- lapply(seq_along(ns_grid), function(i) {
  s <- score_vs_ref(runs[[i]])
  c(list(ns = ns_grid[i], time_s = round(runs[[i]]$dt, 2),
         speedup = round(ref$dt / runs[[i]]$dt, 2)), s)
})

summary <- list(
  setup = list(seed = SEED, locations_in_pool = length(keep_locs),
               treatment_periods = CFG$treatment_periods, N = CFG$N,
               effect_sizes = CFG$effect_size, parallel = FALSE,
               n_simulations_per_run = rows[[1]]$n_cells,
               reference = "ns=1000 (Meta default); lower ns scored against it"),
  sweep = rows)
write_json(summary, file.path(resdir, "market_selection_compare.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8)

# Save the reference ranked table for inspection.
write.csv(bm_ref, file.path(resdir, "market_selection_bestmarkets_default.csv"), row.names = FALSE)

cat("\n--- market selection: ns sweep vs ns=1000 reference ---\n")
cat(sprintf("  sims/run: %d\n", rows[[1]]$n_cells))
cat(sprintf("%-6s %8s %8s %12s %12s %10s\n", "ns", "time_s", "speedup",
            "sig.agree", "top5overlap", "MDEagree"))
for (r in rows) cat(sprintf("%-6d %8.1f %8s %11.0f%% %11.0f%% %9.0f%%\n",
    r$ns, r$time_s, paste0(r$speedup, "x"),
    100 * r$power_significance_agreement, 100 * r$bestmarkets_top5_overlap,
    100 * r$bestmarkets_selected_mde_agreement))
cat("\nWrote", file.path(resdir, "market_selection_compare.json"), "\n")
