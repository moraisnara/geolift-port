# Split the 18s GeoLift() call into: bare SCM fit  vs  permutation inference.
suppressMessages(library(GeoLift))
suppressMessages(library(augsynth))
data(GeoLift_Test)
d <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)
# Mark treatment exactly as GeoLift does.
d$D <- ifelse(d$location %in% c("chicago","portland") & d$time >= 91, 1, 0)

tt <- function(label, expr) { t0<-proc.time()[["elapsed"]]; v<-force(expr)
  cat(sprintf("[time] %-26s %7.3f s\n", label, proc.time()[["elapsed"]]-t0)); v }

fit <- tt("bare augsynth fit", augsynth(Y ~ D, unit = location, time = time,
          data = d, t_int = 91, progfunc = "None", scm = TRUE, fixedeff = TRUE))

tt("predict (counterfactual)", predict(fit, att = TRUE))

tt("permutation p-value", augsynth:::compute_permute_pval(
   data = NULL, ascm = fit, t_int = 91, n_perms = 5000,
   q = 1, side_of_test = "two_sided"))

# Full GeoLift for reference
tt("FULL GeoLift()", GeoLift(Y_id="Y", data=d[,c("location","time","Y")],
   locations=c("chicago","portland"), treatment_start_time=91, treatment_end_time=105))
