# =============================================================================
# time_market_inner_R.R  —  time R's market-selection inference loop on the EXACT
# combos the Python port runs (apples-to-apples). Mirrors the Python port's plain
# sequential loop over the inference cell (GeoLift:::pvalueCalc), no foreach.
#
#   results/ms_R_inner_powercurves.csv   (per-cell pvalue/att/l2/power)
#   results/ms_R_inner_time.json         (wall time + workload)
# Run from Geolift/ root:  Rscript exploration/scripts/time_market_inner_R.R
#
# NOTE: standalone repeated augsynth/osqp calls occasionally segfault on Windows
# R (a BLAS/osqp interaction GeoLift's own foreach path avoids). It is benign and
# intermittent — just re-run; results are deterministic given the seed.
# =============================================================================
suppressMessages({library(GeoLift); library(jsonlite); library(dplyr)})

B <- getwd()   # run from the Geolift/ repo root (see header)
SEED <- 42; TP <- 14; ES <- seq(-0.1, 0.1, 0.05); NS <- 1000; ALPHA <- 0.1
resdir <- file.path(B, "exploration/results")

sub <- read.csv(file.path(B, "exploration/data/ms_subset_panel.csv"), stringsAsFactors = FALSE)
d <- GeoDataRead(data = sub, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)
pc <- read.csv(file.path(resdir, "ms_truth_powercurves.csv"), stringsAsFactors = FALSE)
cm <- do.call(rbind, strsplit(unique(pc$location), ", "))
pvc <- getFromNamespace("pvalueCalc", "GeoLift"); sf <- function(x) sum(abs(x))
mt <- max(d$time)
cells <- expand.grid(combo = seq_len(nrow(cm)), es = ES); n <- nrow(cells)
pv <- numeric(n); att <- numeric(n); l2 <- numeric(n)

set.seed(SEED)
t0 <- proc.time()[["elapsed"]]
for (i in seq_len(n)) {
  r <- suppressMessages(pvc(data = d, sim = 1, max_time = mt, tp = TP, es = cells$es[i],
    locations = as.list(cm[cells$combo[i], ]), cpic = 1, X = c(), type = "pValue",
    normalize = FALSE, fixed_effects = TRUE, stat_func = sf, model = "none",
    conformal_type = "iid", ns = NS))
  pv[i] <- r$pvalue; att[i] <- r$att_estimator; l2[i] <- r$scaled_l2_imbalance
}
dt <- proc.time()[["elapsed"]] - t0

out <- data.frame(location = apply(cm[cells$combo, ], 1, paste, collapse = ", "),
                  duration = TP, EffectSize = cells$es, pvalue = pv,
                  att_estimator = att, ScaledL2Imbalance = l2,
                  power = as.integer(pv < ALPHA))
write.csv(out, file.path(resdir, "ms_R_inner_powercurves.csv"), row.names = FALSE)
write_json(list(seed = SEED, n_combos = nrow(cm), N = ncol(cm), effect_sizes = ES,
                ns = NS, tp = TP, alpha = ALPHA, n_cells = n,
                engine = "R pvalueCalc loop (sequential)", time_s = round(dt, 3)),
           file.path(resdir, "ms_R_inner_time.json"), auto_unbox = TRUE, pretty = TRUE)
cat(sprintf("R inner: %d cells in %.2fs (%.1f ms/cell)\n", n, dt, 1000 * dt / n))
