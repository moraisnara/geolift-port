# =============================================================================
# profile_market_selection.R  —  where does the time go in the POWER path,
# and what is the CEILING on the ns lever?
#
# The scaled ns sweep showed lowering ns only buys ~1.4x in GeoLiftMarketSelection.
# Per-call Rprof shows why: only PART of the per-sim work scales with ns (the
# permutation draws + the per-resample test-statistic algebra). The rest — one
# outer augsynth fit, the one-time counterfactual, input re-validation, data
# wrangling — is FIXED per simulation, independent of ns.
#
# Rather than rely on fragile function-name buckets, we measure it directly:
# run the SAME unmodified GeoLiftMarketSelection at several ns and fit
#     time(ns) = fixed + slope * ns
# The intercept `fixed` is the per-sim cost that ns can never remove, so the
# theoretical ceiling on the ns lever is total_time(ns=1000) / fixed.
#
# Writes exploration/results/market_selection_profile.json.
# Run from Geolift/ root:  Rscript exploration/scripts/profile_market_selection.R
# =============================================================================
suppressMessages(library(GeoLift))
suppressMessages(library(jsonlite))

SEED <- 42
resdir <- "exploration/results"; dir.create(resdir, showWarnings = FALSE, recursive = TRUE)

panel <- read.csv("exploration/data/realistic_panel.csv", stringsAsFactors = FALSE)
d <- GeoDataRead(data = panel, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

CFG <- list(data = d, treatment_periods = c(14), N = c(3),
            effect_size = seq(-0.1, 0.1, 0.05), parallel = FALSE, print = FALSE)

time_at <- function(ns_val) {
  set.seed(SEED)
  t0 <- proc.time()[["elapsed"]]
  invisible(suppressMessages(do.call(GeoLift::GeoLiftMarketSelection, c(CFG, list(ns = ns_val)))))
  proc.time()[["elapsed"]] - t0
}

# ns points spanning the range, incl. near-zero to expose the fixed-cost floor.
ns_grid <- c(1000, 500, 250, 50, 10)
times <- sapply(ns_grid, time_at)

# Linear fit time = fixed + slope*ns  (least squares, no extra deps).
fit    <- lm(times ~ ns_grid)
fixed  <- as.numeric(coef(fit)[1])            # intercept = per-sim work ns can't remove
slope  <- as.numeric(coef(fit)[2])            # seconds per unit ns
t_ref  <- times[ns_grid == 1000]
ns_scaled_frac <- (slope * 1000) / t_ref      # fraction of ns=1000 time that scales with ns
ceiling_speedup <- t_ref / fixed              # best possible ns-lever speedup (ns -> 0)

points <- lapply(seq_along(ns_grid), function(i)
  list(ns = ns_grid[i], time_s = round(times[i], 2),
       speedup_vs_1000 = round(t_ref / times[i], 2)))

out <- list(
  context = "GeoLiftMarketSelection on realistic panel, N=3, sequential, time vs ns",
  points = points,
  model = list(
    fixed_per_run_s = round(fixed, 2),
    slope_s_per_ns  = signif(slope, 3),
    ns_scaled_fraction_at_1000 = round(ns_scaled_frac, 3),
    ns_lever_ceiling_speedup   = round(ceiling_speedup, 2),
    interpretation = paste0(
      "Only ~", round(100 * ns_scaled_frac), "% of the ns=1000 runtime scales with ns; ",
      "the other ~", round(100 * (1 - ns_scaled_frac)), "% is fixed per-sim work (outer fit, ",
      "counterfactual, validation, wrangling). So the ns lever alone cannot exceed ~",
      round(ceiling_speedup, 1), "x no matter how low ns goes. Bigger wins require cutting the ",
      "FIXED work: fit-reuse across effect sizes and parallelism.")))
write_json(out, file.path(resdir, "market_selection_profile.json"),
           auto_unbox = TRUE, pretty = TRUE)

cat("\n--- power-path profile: time vs ns (GeoLiftMarketSelection) ---\n")
for (p in points) cat(sprintf("  ns=%-5d %6.2fs  (%.2fx vs ns=1000)\n", p$ns, p$time_s, p$speedup_vs_1000))
cat(sprintf("\n  fixed per-run floor : %.2fs\n", fixed))
cat(sprintf("  ns-scaled fraction  : %.0f%% of ns=1000 time\n", 100 * ns_scaled_frac))
cat(sprintf("  ns-lever ceiling    : %.2fx (even at ns->0)\n", ceiling_speedup))
cat("\nWrote", file.path(resdir, "market_selection_profile.json"), "\n")
