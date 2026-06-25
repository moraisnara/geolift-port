# =============================================================================
# ms_sweep_common.R — shared ns-sweep + scoring used by run_market_selection.R
# and run_market_selection_scaled.R. Both run the SAME unmodified
# GeoLift::GeoLiftMarketSelection, sequential + fixed seed, varying only `ns`;
# ns=1000 (Meta default) is the reference and every lower ns is scored against it.
#
# Not a standalone script — it is sourced. The caller supplies CFG (the
# GeoLiftMarketSelection args, sans `ns`), the ns grid and SEED.
# =============================================================================
suppressMessages({library(GeoLift); library(dplyr)})

# One GeoLiftMarketSelection run at a given ns (tryCatch so one bad ns doesn't
# abort the whole sweep — the scaled panel can occasionally segfault a cell).
run_ms <- function(CFG, ns_val, SEED) {
  set.seed(SEED)
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(
    suppressMessages(do.call(GeoLift::GeoLiftMarketSelection, c(CFG, list(ns = ns_val)))),
    error = function(e) { message("  ERROR at ns=", ns_val, ": ", conditionMessage(e)); NULL })
  list(out = out, dt = proc.time()[["elapsed"]] - t0)
}

# Score one run's PowerCurves / BestMarkets against the reference run.
score_vs_ref <- function(r, pc_ref, bm_ref) {
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

# Run the full sweep over ns_grid. Returns list(rows = per-ns scored rows,
# ref = the ns_grid[1] reference run) — ns_grid[1] must be the reference (1000).
run_ns_sweep <- function(CFG, ns_grid, SEED) {
  runs <- list()
  for (i in seq_along(ns_grid)) {
    cat(sprintf("[%d/%d] ns=%d ...\n", i, length(ns_grid), ns_grid[i])); flush.console()
    runs[[i]] <- run_ms(CFG, ns_grid[i], SEED)
  }
  ref <- runs[[1]]
  pc_ref <- ref$out$PowerCurves %>% dplyr::arrange(location, duration, EffectSize)
  bm_ref <- ref$out$BestMarkets
  rows <- lapply(seq_along(ns_grid), function(i) {
    s <- score_vs_ref(runs[[i]], pc_ref, bm_ref)
    c(list(ns = ns_grid[i], time_s = round(runs[[i]]$dt, 2),
           speedup = round(ref$dt / runs[[i]]$dt, 2)), s)
  })
  list(rows = rows, ref = ref)
}

# NA-safe table print of the scored rows.
print_sweep <- function(rows) {
  pct <- function(x) if (is.na(x)) "NA" else paste0(round(100 * x), "%")
  cat(sprintf("%-6s %8s %8s %12s %12s %10s\n", "ns", "time_s", "speedup",
              "sig.agree", "top5overlap", "MDEagree"))
  for (r in rows) cat(sprintf("%-6d %8.1f %8s %11s %11s %9s\n",
      r$ns, r$time_s, paste0(r$speedup, "x"),
      pct(r$power_significance_agreement), pct(r$bestmarkets_top5_overlap),
      pct(r$bestmarkets_selected_mde_agreement)))
}
