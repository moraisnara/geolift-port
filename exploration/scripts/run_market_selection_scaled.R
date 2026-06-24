# =============================================================================
# run_market_selection_scaled.R  —  the same ns sweep as run_market_selection.R,
# but on the LARGER realistic panel (exploration/data/realistic_panel.csv) to test
# whether the small-config finding holds at scale: does the ns speedup grow with
# data size, do point estimates stay invariant, does the ranking stay stable?
#
# Same principle: one UNMODIFIED GeoLift::GeoLiftMarketSelection, only `ns` differs.
# Sequential + fixed seed => reproducible and isolates the ns effect (parallelism
# is an orthogonal, ~linear multiplier on top).
#
# Run from Geolift/ root:  Rscript exploration/scripts/run_market_selection_scaled.R
# =============================================================================
suppressMessages(library(GeoLift))
suppressMessages(library(jsonlite))
suppressMessages(library(dplyr))
source("geolift_r_modified/fast_presets.R")

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

run_ms <- function(ns_val) {
  set.seed(SEED)
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(
    suppressMessages(do.call(GeoLift::GeoLiftMarketSelection, c(CFG, list(ns = ns_val)))),
    error = function(e) { message("  ERROR at ns=", ns_val, ": ", conditionMessage(e)); NULL })
  list(out = out, dt = proc.time()[["elapsed"]] - t0)
}

ns_grid <- c(1000, 500, 250)
runs <- list()
for (i in seq_along(ns_grid)) {
  cat(sprintf("[%d/%d] ns=%d ...\n", i, length(ns_grid), ns_grid[i])); flush.console()
  runs[[i]] <- run_ms(ns_grid[i])
}

ref <- runs[[1]]
pc_ref <- ref$out$PowerCurves %>% dplyr::arrange(location, duration, EffectSize)
bm_ref <- ref$out$BestMarkets

score_vs_ref <- function(r) {
  if (is.null(r$out)) return(list(n_cells = NA, power_significance_agreement = NA,
    avg_att_max_abs_diff = NA, scaled_l2_imbalance_max_abs_diff = NA,
    bestmarkets_top5_overlap = NA, bestmarkets_selected_mde_agreement = NA, n_markets = NA))
  pc <- r$out$PowerCurves %>% dplyr::arrange(location, duration, EffectSize)
  j <- dplyr::inner_join(pc_ref, pc, by = c("location", "duration", "EffectSize"),
                         suffix = c("_ref", "_x"))
  bm <- r$out$BestMarkets
  top_k <- min(5, nrow(bm_ref), nrow(bm))
  top_overlap <- if (top_k > 0)
    length(intersect(head(bm_ref$location, top_k), head(bm$location, top_k))) / top_k else NA
  mde_join <- dplyr::inner_join(bm_ref[, c("location", "EffectSize")],
    bm[, c("location", "EffectSize")], by = "location", suffix = c("_ref", "_x"))
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
write.csv(bm_ref, file.path(resdir, "market_selection_scaled_bestmarkets.csv"), row.names = FALSE)

cat("\n--- SCALED market selection: ns sweep vs ns=1000 reference ---\n")
cat(sprintf("  dataset: %d locations x %d periods, %d sims/run\n",
            info$n_locations, info$n_periods, rows[[1]]$n_cells))
cat(sprintf("%-6s %8s %8s %12s %12s %10s\n", "ns", "time_s", "speedup",
            "sig.agree", "top5overlap", "MDEagree"))
for (r in rows) cat(sprintf("%-6d %8.1f %8s %11s %11s %9s\n",
    r$ns, r$time_s, paste0(r$speedup, "x"),
    ifelse(is.na(r$power_significance_agreement), "NA", paste0(round(100*r$power_significance_agreement), "%")),
    ifelse(is.na(r$bestmarkets_top5_overlap), "NA", paste0(round(100*r$bestmarkets_top5_overlap), "%")),
    ifelse(is.na(r$bestmarkets_selected_mde_agreement), "NA", paste0(round(100*r$bestmarkets_selected_mde_agreement), "%"))))
cat("\nWrote", file.path(resdir, "market_selection_scaled_compare.json"), "\n")
