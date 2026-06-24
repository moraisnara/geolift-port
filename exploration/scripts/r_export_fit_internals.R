# Export the exact augsynth fit internals so the Python port can match them precisely.
suppressMessages(library(jsonlite))
objs <- readRDS("exploration/ground_truth/geolift_objects.rds")
fit  <- objs$none$results          # augsynth object
d    <- objs$none$results$data     # internal data list

internals <- list(
  weight_names = dimnames(fit$weights)[[1]],
  weights      = as.numeric(fit$weights[, 1]),
  t_int        = fit$t_int,
  progfunc     = fit$progfunc,
  fixedeff     = fit$fixedeff,
  scm          = fit$scm,
  data_struct  = names(d)
)
str(d, max.level = 2)
writeLines(toJSON(internals, auto_unbox = TRUE, pretty = TRUE, digits = 12),
           "exploration/ground_truth/fit_internals.json")
cat("\nWrote fit_internals.json with", length(internals$weight_names), "donor weights\n")
cat("treated names check / first donors:\n")
print(head(internals$weight_names, 5))
