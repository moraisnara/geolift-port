# =============================================================================
# bench_scaling_R.R  —  one point of the R-vs-Python scaling benchmark.
#
# Usage:  Rscript exploration/scripts/bench_scaling_R.R <n_locations> [ns]
#   ns = 1000 (default)  -> "R original" series (exact Meta default)
#   ns = 100             -> "R modified" series (MARKETSEL_FAST preset; a PARAMETER
#                           passed to the unmodified engine, not a fork)
# Times R's market-selection inference loop (GeoLift:::pvalueCalc) over ALL pairs
# of the first <n_locations> markets x 3 effect sizes, and writes:
#   results/bench/bench_R_L<n>_ns<ns>.json  {n_locations, n_combos, n_cells, time_s}
#   results/bench/bench_R_L<n>_ns<ns>.csv   per-cell pvalue/att/l2/power (for agreement)
#
# Run as a FRESH process per point (this is why it takes args): standalone repeated
# augsynth/osqp calls occasionally segfault on Windows R, and a fresh process both
# isolates any crash to one point and avoids cumulative state. Timing per cell is
# independent of the effect-size VALUE, so 3 es suffice to characterise scaling
# while keeping each R point well within the safe regime.
# =============================================================================
suppressMessages({library(GeoLift); library(jsonlite)})

args <- commandArgs(trailingOnly = TRUE)
L <- as.integer(args[1])
NS <- if (length(args) >= 2) as.integer(args[2]) else 1000L
B <- "c:/Users/naral/Desktop/Nara/DP6/Geolift"
SEED <- 42; TP <- 14; ES <- c(-0.1, 0.0, 0.1); ALPHA <- 0.1
outdir <- file.path(B, "exploration/results/bench")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

full <- read.csv(file.path(B, "exploration/data/realistic_panel.csv"), stringsAsFactors = FALSE)
markets <- sprintf("market_%03d", seq_len(L))
sub <- full[full$location %in% markets, c("location", "date", "Y")]
d <- GeoDataRead(data = sub, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

combos <- combn(markets, 2, simplify = FALSE)
pvc <- getFromNamespace("pvalueCalc", "GeoLift"); sf <- function(x) sum(abs(x))
mt <- max(d$time)
grid <- expand.grid(ci = seq_along(combos), es = ES); n <- nrow(grid)
pv <- numeric(n); att <- numeric(n); l2 <- numeric(n); loc <- character(n)

set.seed(SEED)
t0 <- proc.time()[["elapsed"]]
for (i in seq_len(n)) {
  cb <- combos[[grid$ci[i]]]
  r <- suppressMessages(pvc(data = d, sim = 1, max_time = mt, tp = TP, es = grid$es[i],
    locations = as.list(cb), cpic = 1, X = c(), type = "pValue", normalize = FALSE,
    fixed_effects = TRUE, stat_func = sf, model = "none", conformal_type = "iid", ns = NS))
  pv[i] <- r$pvalue; att[i] <- r$att_estimator; l2[i] <- r$scaled_l2_imbalance
  loc[i] <- paste(cb, collapse = ", ")
}
dt <- proc.time()[["elapsed"]] - t0

out <- data.frame(location = loc, EffectSize = grid$es, pvalue = pv,
                  att_estimator = att, ScaledL2Imbalance = l2,
                  power = as.integer(pv < ALPHA))
write.csv(out, file.path(outdir, sprintf("bench_R_L%d_ns%d.csv", L, NS)), row.names = FALSE)
write_json(list(n_locations = L, n_combos = length(combos), n_cells = n,
                tp = TP, ns = NS, effect_sizes = ES, seed = SEED,
                time_s = round(dt, 3), ms_per_cell = round(1000 * dt / n, 2)),
           file.path(outdir, sprintf("bench_R_L%d_ns%d.json", L, NS)),
           auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("L=%d ns=%d: %d combos, %d cells in %.2fs (%.1f ms/cell)\n",
            L, NS, length(combos), n, dt, 1000 * dt / n))
