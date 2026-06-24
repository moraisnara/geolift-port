# Load the saved GeoLift objects and write a correct benchmark JSON for the Python port.
suppressMessages(library(jsonlite))
outdir <- "exploration/ground_truth"
objs <- readRDS(file.path(outdir, "geolift_objects.rds"))

extract <- function(obj) {
  inf <- obj$inference            # data.frame: ATT, Perc.Lift, pvalue, Lower/Upper.Conf.Int
  r   <- obj$results              # augsynth fit object
  list(
    ATT                 = as.numeric(inf$ATT),
    percent_lift        = as.numeric(inf$Perc.Lift),
    pvalue              = as.numeric(inf$pvalue),
    lower_conf_int      = as.numeric(inf$Lower.Conf.Int),
    upper_conf_int      = as.numeric(inf$Upper.Conf.Int),
    incremental         = as.numeric(obj$incremental),
    l2_imbalance        = as.numeric(r$l2_imbalance),
    scaled_l2_imbalance = as.numeric(r$scaled_l2_imbalance),
    treatment_start     = as.numeric(obj$TreatmentStart),
    treatment_end       = as.numeric(obj$TreatmentEnd),
    progfunc            = r$progfunc,
    weights             = setNames(as.numeric(r$weights[,1]), dimnames(r$weights)[[1]]),
    y_hat               = as.numeric(obj$y_hat),     # counterfactual, 105 periods
    att_series          = as.numeric(obj$ATT)        # per-period ATT, 105 periods
  )
}

bench <- list(
  meta = list(
    example = "GeoLift Walkthrough: locations=c(chicago,portland), t=91..105",
    augsynth_commit = "06415b4 (pre-PR#88)",
    doc_expected = list(percent_lift = 5.4, incremental = 4667,
                        att_per_day = 155.556, pvalue_approx = 0.01)
  ),
  None = extract(objs$none),
  best = extract(objs$best)
)
writeLines(toJSON(bench, auto_unbox = TRUE, pretty = TRUE, digits = 10, na = "null"),
           file.path(outdir, "benchmark.json"))
cat("Wrote benchmark.json\n")
cat(sprintf("None: ATT=%.3f  lift=%.2f%%  incr=%.0f  p=%.3f  L2=%.2f\n",
    bench$None$ATT, bench$None$percent_lift, bench$None$incremental,
    bench$None$pvalue, bench$None$l2_imbalance))
