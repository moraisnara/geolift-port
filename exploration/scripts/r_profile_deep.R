# Deep profile of a single GeoLift() call using R's sampling profiler.
suppressMessages(library(GeoLift))
data(GeoLift_Test)
d <- GeoDataRead(data = GeoLift_Test, date_id = "date", location_id = "location",
                 Y_id = "Y", X = c(), format = "yyyy-mm-dd", summary = FALSE)

prof <- tempfile(fileext = ".out")
Rprof(prof, interval = 0.01, line.profiling = FALSE)
GeoTest <- GeoLift(Y_id = "Y", data = d,
                   locations = c("chicago", "portland"),
                   treatment_start_time = 91, treatment_end_time = 105)
Rprof(NULL)

s <- summaryRprof(prof)
cat("\n-- TOP BY SELF TIME (where CPU actually sits) --\n")
print(head(s$by.self, 25))
cat("\n-- TOP BY TOTAL TIME (self+children) --\n")
print(head(s$by.total, 30))
cat("\nsampling time:", s$sampling.time, "s\n")
