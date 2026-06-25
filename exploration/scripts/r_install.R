# =============================================================================
# r_install.R — install the pinned R dependencies for this study, in order:
#   1. augsynth, pinned to commit 06415b4 (the SCM engine GeoLift wraps; pinned
#      to before augsynth PR #88 / 2026-03-11, which broke summary.augsynth for
#      GeoLift 2.7.5).
#   2. GeoLift, from the vendored source in geolift_r_original/.
#
# Run from Geolift/ root:  Rscript exploration/scripts/r_install.R
# =============================================================================

# Rtools is only needed (and only present) on Windows; harmless no-op elsewhere.
if (.Platform$OS.type == "windows") {
  Sys.setenv(PATH = paste(
    "C:\\RBuildTools\\4.4\\usr\\bin",
    "C:\\RBuildTools\\4.4\\x86_64-w64-mingw32.static.posix\\bin",
    Sys.getenv("PATH"), sep = ";"))
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

# ---- 1. augsynth (pinned) ---------------------------------------------------
cat("build tools:", pkgbuild::has_build_tools(debug = FALSE), "\n")
remotes::install_github("ebenmichael/augsynth@06415b4da556d5f8c9521abdced7b82b1487947c",
                        upgrade = "never", quiet = FALSE,
                        build_vignettes = FALSE, force = TRUE)
cat("augsynth available:", requireNamespace("augsynth", quietly = TRUE), "\n")

# ---- 2. GeoLift (from vendored source) --------------------------------------
res <- try(install.packages("geolift_r_original", repos = NULL, type = "source",
                            INSTALL_opts = "--no-staged-install"),
           silent = FALSE)
cat("FINAL GeoLift available:", requireNamespace("GeoLift", quietly = TRUE), "\n")
