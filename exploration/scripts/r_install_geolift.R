Sys.setenv(PATH = paste(
  "C:\\RBuildTools\\4.4\\usr\\bin",
  "C:\\RBuildTools\\4.4\\x86_64-w64-mingw32.static.posix\\bin",
  Sys.getenv("PATH"), sep = ";"))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install directly from source dir, no vignette build, capture errors.
res <- try(install.packages("geolift_r_original", repos = NULL, type = "source",
                            INSTALL_opts = "--no-staged-install"),
           silent = FALSE)
cat("FINAL GeoLift available:", requireNamespace("GeoLift", quietly = TRUE), "\n")
