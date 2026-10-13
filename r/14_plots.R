
# 14_plots.R -- OWNER: Vishi
# Final report figures.
# Run from repository root:
# Rscript r/14_plots.R

source("r/common.R")

residual_file <- file.path(PROCESSED, "residuals.csv")
conjunction_file <- file.path(PROCESSED, "conjunctions.csv")
pc_file <- file.path(PROCESSED, "pc_exact.csv")

if (!file.exists(residual_file)) {
  stop("Missing data/processed/residuals.csv")
}

dat <- read.csv(residual_file, stringsAsFactors = FALSE)

required <- c(
  "age_days",
  "d_radial_km",
  "d_along_km",
  "d_cross_km"
)
missing <- setdiff(required, names(dat))
if (length(missing) > 0) {
  stop("Missing columns in residuals.csv: ",
       paste(missing, collapse = ", "))
}

dat <- dat[complete.cases(dat[, required]), required]
dat <- dat[
  is.finite(dat$age_days) & dat$age_days > 0 &
    is.finite(dat$d_radial_km) &
    is.finite(dat$d_along_km) &
    is.finite(dat$d_cross_km),
]

directions <- c(
  radial = "d_radial_km",
  along = "d_along_km",
  cross = "d_cross_km"
)

# 1. Empirical error growth against TLE age.
for (direction in names(directions)) {
  xcol <- directions[[direction]]

  breaks <- unique(as.numeric(
    quantile(dat$age_days, probs = seq(0, 1, length.out = 21))
  ))
  age_bin <- cut(dat$age_days, breaks = breaks, include.lowest = TRUE)

  groups <- split(
    data.frame(age = dat$age_days, error = dat[[xcol]]),
    age_bin
  )

  summaries <- lapply(groups, function(g) {
    if (nrow(g) == 0) return(NULL)
    data.frame(
      age = median(g$age),
      median_abs_error = median(abs(g$error)),
      n = nrow(g)
    )
  })
  summaries <- Filter(Negate(is.null), summaries)
  summaries <- do.call(rbind, summaries)

  png(
    file.path(FIGURES, paste0("final_error_growth_", direction, ".png")),
    width = 1000, height = 700
  )
  plot(
    summaries$age, summaries$median_abs_error,
    log = "xy",
    pch = 19,
    xlab = "TLE age (days)",
    ylab = "Median absolute residual (km)",
    main = paste("Error Growth with TLE Age:", direction)
  )
  dev.off()
}

# 2. Normal Q-Q plots for raw residuals.
# These are descriptive diagnostics, not fitted-model validation by themselves.
for (direction in names(directions)) {
  x <- dat[[directions[[direction]]]]

  png(
    file.path(FIGURES, paste0("final_qq_", direction, ".png")),
    width = 1000, height = 700
  )
  qqnorm(x, main = paste("Residual Normal Q-Q Plot:", direction))
  qqline(x, col = "red", lwd = 2)
  dev.off()
}

# 3. Monte Carlo convergence curves, if a suitable output file exists.
# This section only uses files whose required columns can be identified.
candidate_files <- list.files(
  PROCESSED,
  pattern = "\\.csv$",
  full.names = TRUE
)

convergence_found <- FALSE

for (f in candidate_files) {
  nm <- tolower(basename(f))
  if (!grepl("monte|converg", nm)) next

  d <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || nrow(d) < 2) next

  sample_col <- intersect(
    c("n_samples", "samples", "sample_size", "iteration"),
    names(d)
  )
  estimate_col <- intersect(
    c("pc", "pc_estimate", "estimate", "running_estimate"),
    names(d)
  )

  if (length(sample_col) == 0 || length(estimate_col) == 0) next

  png(
    file.path(FIGURES, "final_monte_carlo_convergence.png"),
    width = 1000, height = 700
  )
  plot(
    d[[sample_col[1]]], d[[estimate_col[1]]],
    type = "l",
    xlab = "Sample size / iteration",
    ylab = "Estimated collision probability",
    main = "Monte Carlo Convergence"
  )
  dev.off()
  convergence_found <- TRUE
  break
}

if (!convergence_found) {
  cat(
    "Skipping convergence curve: no processed CSV with recognized ",
    "sample-size and running-estimate columns was found.\n",
    sep = ""
  )
}

# 4. ROC curve, if a processed CSV contains observed labels and scores.
roc_found <- FALSE

for (f in candidate_files) {
  nm <- tolower(basename(f))
  if (!grepl("roc|classif|label|score", nm)) next

  d <- tryCatch(
    read.csv(f, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(d)) next

  label_col <- intersect(
    c("observed", "label", "actual", "truth", "event"),
    names(d)
  )
  score_col <- intersect(
    c("score", "predicted", "prediction", "pc", "risk_score"),
    names(d)
  )

  if (length(label_col) == 0 || length(score_col) == 0) next

  y <- d[[label_col[1]]]
  score <- d[[score_col[1]]]

  if (is.logical(y)) y <- as.integer(y)
  if (is.factor(y)) y <- as.character(y)
  if (is.character(y)) {
    y <- ifelse(tolower(y) %in% c("1", "true", "yes", "event"), 1,
                ifelse(tolower(y) %in% c("0", "false", "no", "none"), 0, NA))
  }

  valid <- is.finite(y) & is.finite(score) & y %in% c(0, 1)
  y <- y[valid]
  score <- score[valid]

  if (length(unique(y)) < 2) next

  # Compute ROC coordinates by threshold without adding dependencies.
  thresholds <- sort(unique(score), decreasing = TRUE)
  thresholds <- c(Inf, thresholds, -Inf)

  tpr <- vapply(thresholds, function(th) {
    sum(score[y == 1] >= th) / sum(y == 1)
  }, numeric(1))

  fpr <- vapply(thresholds, function(th) {
    sum(score[y == 0] >= th) / sum(y == 0)
  }, numeric(1))

  ord <- order(fpr, tpr)

  png(
    file.path(FIGURES, "final_roc_curve.png"),
    width = 1000, height = 700
  )
  plot(
    fpr[ord], tpr[ord],
    type = "l", lwd = 2,
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "False positive rate",
    ylab = "True positive rate",
    main = "ROC Curve"
  )
  abline(0, 1, lty = 2)
  dev.off()

  roc_found <- TRUE
  break
}

if (!roc_found) {
  cat(
    "Skipping ROC curve: no suitable observed-label and prediction-score ",
    "data were found. Do not create an ROC curve without both.\n",
    sep = ""
  )
}

# 5. Poisson count diagnostic, only if actual event-time data are available.
# A single high-risk event is not sufficient to assess exponential gaps.
poisson_found <- FALSE

if (file.exists(conjunction_file) && file.exists(pc_file)) {
  con <- read.csv(conjunction_file, stringsAsFactors = FALSE)
  pc <- read.csv(pc_file, stringsAsFactors = FALSE)

  if ("conj_id" %in% names(con) && "conj_id" %in% names(pc) &&
      "pc" %in% names(pc) && "tca_utc" %in% names(con)) {

    matched <- merge(
      con[, c("conj_id", "tca_utc")],
      pc[, c("conj_id", "pc")],
      by = "conj_id"
    )

    matched$tca_utc <- as.POSIXct(
      matched$tca_utc,
      tz = "UTC",
      tryFormats = c(
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%OS"
      )
    )

    matched <- matched[
      is.finite(matched$pc) & !is.na(matched$tca_utc) &
        matched$pc >= ALERT_THRESHOLD,
    ]

    if (nrow(matched) >= 2 &&
        length(unique(matched$tca_utc)) >= 2) {
      event_times <- sort(unique(matched$tca_utc))
      gaps_days <- as.numeric(diff(event_times), units = "days")
      gaps_days <- gaps_days[is.finite(gaps_days) & gaps_days > 0]

      if (length(gaps_days) >= 1) {
        rate_per_day <- 1 / mean(gaps_days)

        png(
          file.path(FIGURES, "final_poisson_gap_diagnostic.png"),
          width = 1000, height = 700
        )
        hist(
          gaps_days,
          probability = TRUE,
          xlab = "Gap between high-risk event times (days)",
          main = "High-Risk Event Gap Diagnostic"
        )
        curve(
          rate_per_day * exp(-rate_per_day * x),
          add = TRUE,
          lwd = 2
        )
        dev.off()

        poisson_found <- TRUE
      }
    }
  }
}

if (!poisson_found) {
  cat(
    "Skipping Poisson gap diagnostic: fewer than two distinct high-risk ",
    "event times are available, or required columns are missing.\n",
    sep = ""
  )
}

cat("\nFinal plotting script completed.\n")
cat("Figures saved in:", FIGURES, "\n")
cat("Generated: error-growth and raw-residual Q-Q plots.\n")
cat("Optional plots are generated only when suitable input data exist.\n")
