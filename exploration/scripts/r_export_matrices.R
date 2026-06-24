objs <- readRDS("exploration/ground_truth/geolift_objects.rds")
sd <- objs$none$results$data$synth_data
od <- "exploration/ground_truth"
write.csv(sd$Z0, file.path(od, "Z0_donors_pre.csv"), row.names = FALSE)       # 90 x 38
write.csv(sd$Z1, file.path(od, "Z1_treated_pre.csv"), row.names = FALSE)      # 90 x 1
write.csv(sd$Y0plot, file.path(od, "Y0_donors_full.csv"), row.names = FALSE)  # 105 x 38
write.csv(sd$Y1plot, file.path(od, "Y1_treated_full.csv"), row.names = FALSE) # 105 x 1
cat("exported. dims Z0:", dim(sd$Z0), " Y0:", dim(sd$Y0plot), "\n")
