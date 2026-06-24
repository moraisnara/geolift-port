# =============================================================================
# dump_market_cells.R  —  ground truth for the Python market-selection port.
#
# Runs the UNMODIFIED GeoLift::GeoLiftMarketSelection on a small subset of the
# realistic panel (fixed seed) and dumps everything the Python version must
# reproduce:
#   * the exact panel subset            -> data/ms_subset_panel.csv
#   * the full PowerCurves R produced    -> results/ms_truth_powercurves.csv
#       (location-combo, duration, EffectSize, pvalue, att_estimator,
#        ScaledL2Imbalance, ...)  -> the cells Python must match
#   * per-cell DETERMINISTIC internals for a handful of cells
#     (weights, att_estimator, scaled_l2_imbalance, the conformal residual
#      vector and the observed post-period statistic) -> results/ms_truth_cells.json
#     These are seed-independent, so they pin the Python fit exactly (~1e-6);
#     the p-value itself is Monte-Carlo and only needs to agree in distribution.
#
# Run from Geolift/ root:  Rscript exploration/scripts/dump_market_cells.R
# =============================================================================
suppressMessages({library(GeoLift); library(jsonlite); library(dplyr)})

SEED <- 42
datadir <- "exploration/data"; resdir <- "exploration/results"
dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

# ---- small, fixed workload (keep R fast; big enough to be meaningful) --------
N_MARKETS <- 24
TP        <- 14
ES        <- seq(-0.1, 0.1, 0.05)
NS        <- 1000
NVAL      <- c(2)          # combination size
ALPHA     <- 0.1

full <- read.csv(file.path(datadir, "realistic_panel.csv"), stringsAsFactors = FALSE)
keep <- sprintf("market_%03d", seq_len(N_MARKETS))
sub  <- full[full$location %in% keep, c("location", "date", "Y")]
write.csv(sub, file.path(datadir, "ms_subset_panel.csv"), row.names = FALSE)

d <- GeoDataRead(data = sub, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

# ---- run the real market selection -> PowerCurves ground truth ---------------
set.seed(SEED)
ms <- suppressMessages(GeoLift::GeoLiftMarketSelection(
  data = d, treatment_periods = c(TP), N = NVAL, effect_size = ES,
  alpha = ALPHA, parallel = FALSE, print = FALSE, ns = NS))

pc <- ms$PowerCurves
write.csv(pc, file.path(resdir, "ms_truth_powercurves.csv"), row.names = FALSE)
cat(sprintf("PowerCurves: %d cells, %d distinct combos\n",
            nrow(pc), length(unique(pc$location))))
cat("PowerCurves columns:", paste(names(pc), collapse = ", "), "\n")

# ---- per-cell deterministic internals (pin the Python fit) -------------------
# Replicate pvalueCalc's fit + conformal residual construction for a few cells,
# extracting the seed-INDEPENDENT pieces so Python can match to ~1e-6.
fn_treatment <- getFromNamespace("fn_treatment", "GeoLift")
stat_func <- function(x) sum(abs(x))            # two_sided (MarketSelection default)
max_time <- max(d$time)

dump_cell <- function(locs, tp, es, seed = SEED) {
  treatment_start_time <- max_time - tp - 1 + 2     # sim = lookback = 1
  treatment_end_time   <- treatment_start_time + tp - 1
  da <- fn_treatment(d, locations = as.list(locs), treatment_start_time, treatment_end_time)
  da$Y_inc <- da$Y
  da$Y_inc[da$D == 1] <- da$Y_inc[da$D == 1] * (1 + es)
  ascm <- augsynth::augsynth(Y_inc ~ D, unit = location, time = time, data = da,
                             t_int = treatment_start_time, progfunc = "none",
                             scm = TRUE, fixedeff = TRUE)
  ave_treatment_convs <- sum(ascm$data$y[which(ascm$data$trt == 1), ]) / sum(ascm$data$trt == 1)
  ave_pred_control    <- predict(ascm)[treatment_start_time:treatment_end_time]
  ave_incremental     <- ave_treatment_convs - sum(ave_pred_control)
  att <- ave_incremental / tp

  # conformal residual vector (deterministic) + observed post statistic
  wide <- ascm$data
  nw <- wide; nw$X <- cbind(wide$X, wide$y); nw$y <- matrix(1, nrow(wide$X), 1)
  post_length <- ncol(wide$y)
  cts <- augsynth:::compute_permute_test_stats(nw, ascm, 0, post_length,
                                               "iid", 1, NS, stat_func)
  t0 <- ncol(nw$X) - post_length; tpost <- t0 + post_length
  obs_stat <- stat_func(cts$resids[(t0 + 1):tpost])
  set.seed(seed)
  pval <- augsynth:::compute_permute_pval(nw, ascm, 0, post_length, "iid", 1, NS, stat_func)

  list(locations = unlist(locs), tp = tp, es = es,
       treatment_start = treatment_start_time,
       weights = as.numeric(ascm$weights),
       donor_units = rownames(ascm$data$X)[ascm$data$trt == 0],
       att_estimator = att, scaled_l2_imbalance = ascm$scaled_l2_imbalance,
       t0 = t0, tpost = tpost, post_length = post_length,
       resids = as.numeric(cts$resids), obs_stat = obs_stat,
       pvalue_seed42 = pval)
}

# pick the first tested combo across a few effect sizes + a second combo at es=0.1
combos <- strsplit(unique(pc$location), ", ")
sel <- list(
  list(locs = combos[[1]], es = 0.0),
  list(locs = combos[[1]], es = 0.1),
  list(locs = combos[[1]], es = -0.05),
  list(locs = combos[[min(2, length(combos))]], es = 0.1))
cells <- lapply(sel, function(s) dump_cell(s$locs, TP, s$es))

write_json(list(
  setup = list(seed = SEED, n_markets = N_MARKETS, tp = TP, effect_sizes = ES,
               ns = NS, N = NVAL, alpha = ALPHA, side_of_test = "two_sided",
               conformal_type = "iid", model = "none", fixed_effects = TRUE,
               max_time = max_time),
  cells = cells),
  file.path(resdir, "ms_truth_cells.json"), auto_unbox = TRUE, pretty = TRUE, digits = 10)

cat("Wrote ms_subset_panel.csv, ms_truth_powercurves.csv, ms_truth_cells.json\n")
cat(sprintf("Sample cell[1]: att=%.4f l2=%.4f pval@42=%.4f resid_len=%d\n",
            cells[[1]]$att_estimator, cells[[1]]$scaled_l2_imbalance,
            cells[[1]]$pvalue_seed42, length(cells[[1]]$resids)))
