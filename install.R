# Run once:  Rscript install.R
pkgs <- c("readr", "dplyr", "ggplot2", "MASS", "fitdistrplus", "boot")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
cat("R packages ready.\n")
