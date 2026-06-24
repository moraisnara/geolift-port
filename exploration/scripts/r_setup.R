# Ensure Rtools44 toolchain is visible, then install/verify deps.
Sys.setenv(PATH = paste(
  "C:\\RBuildTools\\4.4\\usr\\bin",
  "C:\\RBuildTools\\4.4\\x86_64-w64-mingw32.static.posix\\bin",
  Sys.getenv("PATH"), sep = ";"))
options(repos = c(CRAN = "https://cloud.r-project.org"))

ip <- rownames(installed.packages())
cat("build tools:", pkgbuild::has_build_tools(debug = FALSE), "\n")
cat("augsynth:", "augsynth" %in% ip, "\n")
cat("gsynth:", "gsynth" %in% ip, "\n")
cat("GeoLift:", "GeoLift" %in% ip, "\n")

if (!("GeoLift" %in% ip)) {
  cat("== installing GeoLift from local source ==\n")
  remotes::install_local("geolift_r_original", upgrade = "never", quiet = TRUE)
  cat("GeoLift now:", requireNamespace("GeoLift", quietly = TRUE), "\n")
}
