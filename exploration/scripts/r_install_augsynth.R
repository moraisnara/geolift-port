Sys.setenv(PATH = paste(
  "C:\\RBuildTools\\4.4\\usr\\bin",
  "C:\\RBuildTools\\4.4\\x86_64-w64-mingw32.static.posix\\bin",
  Sys.getenv("PATH"), sep = ";"))
options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

cat("build tools:", pkgbuild::has_build_tools(debug = FALSE), "\n")
# Pin to commit before PR #88 (2026-03-11) which broke summary.augsynth for GeoLift 2.7.5.
remotes::install_github("ebenmichael/augsynth@06415b4da556d5f8c9521abdced7b82b1487947c",
                        upgrade = "never", quiet = FALSE,
                        build_vignettes = FALSE, force = TRUE)
cat("FINAL augsynth available:", requireNamespace("augsynth", quietly = TRUE), "\n")
