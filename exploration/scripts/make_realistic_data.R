# =============================================================================
# make_realistic_data.R  —  build a larger, realistic geo panel like GeoLift_Test
#
# GeoLift_Test is 40 locations x 105 daily periods, columns (location, Y, date),
# with heterogeneous market sizes and weekly seasonality. We generate a BIGGER
# panel with the same shape and the structure SCM relies on:
#   * heterogeneous market sizes (log-normal, like real DMAs)
#   * shared latent factors (national trend + weekly/monthly seasonality) with
#     market-specific loadings  -> realistic cross-market CORRELATION (donors)
#   * idiosyncratic noise
#   * a `category` (region) column so markets fall into groups
# Integer conversions, deterministic (fixed seed) so re-runs match.
#
# Output: exploration/data/realistic_panel.csv  (+ realistic_panel_info.json)
# Run from Geolift/ root:  Rscript exploration/scripts/make_realistic_data.R
#   optional env: GEO_M (#locations), GEO_T (#days), GEO_CAT (#categories)
# =============================================================================
suppressMessages(library(jsonlite))

SEED <- 42
set.seed(SEED)

M     <- as.integer(Sys.getenv("GEO_M",   "80"))     # locations
Tn    <- as.integer(Sys.getenv("GEO_T",   "220"))    # daily periods
n_cat <- as.integer(Sys.getenv("GEO_CAT", "6"))      # categories / regions

datadir <- "exploration/data"; dir.create(datadir, showWarnings = FALSE, recursive = TRUE)
dates <- as.Date("2022-01-01") + 0:(Tn - 1)

# ---- shared latent factors (induce cross-market correlation) ----------------
trend <- cumsum(rnorm(Tn, mean = 0.03, sd = 0.4))          # slow national trend
week  <- 1.0 * sin(2 * pi * (1:Tn) / 7) + 0.5 * cos(2 * pi * (1:Tn) / 7)  # weekly
month <- 0.7 * sin(2 * pi * (1:Tn) / 30.4)                 # monthly
Fmat  <- cbind(trend, week, month)                         # Tn x 3

# ---- market attributes ------------------------------------------------------
size  <- rlnorm(M, meanlog = log(1), sdlog = 0.7)          # heterogeneous sizes
cat_id <- sample(seq_len(n_cat), M, replace = TRUE)
# category-level loading centers so markets in a region co-move more strongly
cat_load <- matrix(rnorm(n_cat * 3, mean = 1, sd = 0.25), nrow = n_cat)
load  <- cat_load[cat_id, ] + matrix(rnorm(M * 3, 0, 0.15), nrow = M)
base  <- 320 + 40 * rnorm(M)

# ---- assemble panel ---------------------------------------------------------
rows <- vector("list", M)
for (i in seq_len(M)) {
  signal <- base[i] + as.numeric(Fmat %*% load[i, ]) * 12
  y <- size[i] * signal + rnorm(Tn, 0, 18)                 # idiosyncratic noise
  y <- pmax(round(y), 0)
  rows[[i]] <- data.frame(
    location = sprintf("market_%03d", i),
    category = sprintf("region_%d", cat_id[i]),
    date     = dates,
    Y        = y,
    stringsAsFactors = FALSE)
}
panel <- do.call(rbind, rows)

csv <- file.path(datadir, "realistic_panel.csv")
write.csv(panel, csv, row.names = FALSE)

# ---- record properties (so the report can describe the data, not hardcode) --
# average pairwise correlation across markets (on a wide matrix) — sanity that
# donors are usable for SCM.
wide <- reshape(panel[, c("location", "date", "Y")],
                idvar = "date", timevar = "location", direction = "wide")
cmat <- suppressWarnings(cor(wide[, -1]))
avg_cor <- mean(cmat[upper.tri(cmat)], na.rm = TRUE)

info <- list(
  seed = SEED, n_locations = M, n_periods = Tn, n_categories = n_cat,
  date_start = as.character(min(dates)), date_end = as.character(max(dates)),
  n_rows = nrow(panel),
  Y_min = min(panel$Y), Y_median = median(panel$Y), Y_max = max(panel$Y),
  avg_pairwise_correlation = round(avg_cor, 3),
  compared_to_GeoLift_Test = "40 locations x 105 periods",
  file = csv)
write_json(info, file.path(datadir, "realistic_panel_info.json"),
           auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("Wrote %s\n  %d locations x %d periods x %d categories = %d rows\n",
            csv, M, Tn, n_cat, nrow(panel)))
cat(sprintf("  Y range %d..%d (median %d), avg pairwise corr %.2f\n",
            min(panel$Y), max(panel$Y), as.integer(median(panel$Y)), avg_cor))
