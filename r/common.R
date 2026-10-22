# Shared paths and settings for all R scripts.
# Run every R script from the repo root:  Rscript r/10_error_model.R
ROOT <- normalizePath(getwd())
if (!dir.exists(file.path(ROOT, "data")) || !dir.exists(file.path(ROOT, "r"))) {
  stop("Run R scripts from the repository root (the folder that contains data/ and r/).")
}
PROCESSED <- file.path(ROOT, "data", "processed")
FIGURES   <- file.path(ROOT, "figures")
dir.create(FIGURES, showWarnings = FALSE, recursive = TRUE)

SEED <- 42
ALERT_THRESHOLD <- 1e-4       # keep identical to python/common.py
ROC_LABEL_THRESHOLD <- 1e-5   # Student-t Pc at or above this is a "true risk" label for the ROC (r/14)

not_implemented <- function(script, owner) {
  cat(sprintf("%s: not implemented yet (owner: %s)\n", script, owner))
}
