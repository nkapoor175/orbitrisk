
# 11_gaussian_vs_learned_tests.R -- OWNER: Vishi
# Compare Gaussian and Student-t error models.
# Run from repository root:
# Rscript r/11_gaussian_vs_learned_tests.R

source("r/common.R")

residual_file <- file.path(PROCESSED, "residuals.csv")
params_file <- file.path(PROCESSED, "params.csv")
output_file <- file.path(PROCESSED, "model_comparison_tests.csv")

if (!file.exists(residual_file)) {
  stop("Missing data/processed/residuals.csv. Run the residual pipeline first.")
}

if (!file.exists(params_file)) {
  stop("Missing data/processed/params.csv. Run r/10_error_model.R first.")
}

dat <- read.csv(residual_file, stringsAsFactors = FALSE)
params <- read.csv(params_file, stringsAsFactors = FALSE)

required_data <- c(
  "age_days",
  "d_radial_km",
  "d_along_km",
  "d_cross_km"
)

required_params <- c(
  "direction", "a", "b",
  "student_location", "student_df", "student_scale",
  "gaussian_mean", "gaussian_sd"
)

missing_data <- setdiff(required_data, names(dat))
missing_params <- setdiff(required_params, names(params))

if (length(missing_data) > 0) {
  stop("Missing residual columns: ", paste(missing_data, collapse = ", "))
}

if (length(missing_params) > 0) {
  stop(
    "Missing parameter columns: ",
    paste(missing_params, collapse = ", "),
    ". Rerun r/10_error_model.R using the updated version."
  )
}

directions <- c(
  radial = "d_radial_km",
  along = "d_along_km",
  cross = "d_cross_km"
)

# Calculate empirical distribution-fit diagnostics and paired log scores.
compare_direction <- function(direction, column) {

  p <- params[params$direction == direction, , drop = FALSE]

  if (nrow(p) != 1) {
    stop("Expected exactly one parameter row for: ", direction)
  }

  d <- dat[
    is.finite(dat$age_days) &
      dat$age_days > 0 &
      is.finite(dat[[column]]),
    c("age_days", column),
    drop = FALSE
  ]

  # Recreate the standardization used by r/10_error_model.R.
  sigma_age <- p$a * d$age_days^p$b
  valid <- is.finite(sigma_age) & sigma_age > 0

  z <- d[[column]][valid] / sigma_age[valid]

  if (length(z) < 30) {
    stop("Too few valid standardized observations for: ", direction)
  }

  gm <- p$gaussian_mean
  gs <- p$gaussian_sd
  tl <- p$student_location
  nu <- p$student_df
  ts <- p$student_scale

  if (
    !is.finite(gm) || !is.finite(gs) || gs <= 0 ||
    !is.finite(tl) || !is.finite(nu) || nu <= 0 ||
    !is.finite(ts) || ts <= 0
  ) {
    stop("Invalid fitted distribution parameters for: ", direction)
  }

  # Per-observation log predictive scores.
  log_gaussian <- dnorm(z, mean = gm, sd = gs, log = TRUE)

  log_student <- dt(
    (z - tl) / ts,
    df = nu,
    log = TRUE
  ) - log(ts)

  if (
    any(!is.finite(log_gaussian)) ||
    any(!is.finite(log_student))
  ) {
    stop("Non-finite log scores for direction: ", direction)
  }

  ll_gaussian <- sum(log_gaussian)
  ll_student <- sum(log_student)

  # AIC: Gaussian estimates mean and SD; Student-t estimates
  # location, scale and degrees of freedom.
  aic_gaussian <- 2 * 2 - 2 * ll_gaussian
  aic_student <- 2 * 3 - 2 * ll_student

  # KS statistics compare the empirical standardized residuals
  # with each fitted distribution.
  # Since parameters were estimated from these same observations,
  # the ordinary KS p-values are descriptive, not calibrated tests.
  ks_gaussian <- suppressWarnings(
    ks.test(
      z, "pnorm",
      mean = gm,
      sd = gs
    )
  )

  student_cdf <- function(q) {
    pt((q - tl) / ts, df = nu)
  }

  ks_student <- suppressWarnings(
    ks.test(z, student_cdf)
  )

  # Paired comparison of per-observation log scores.
  # Positive differences favour Student-t.
  score_difference <- log_student - log_gaussian

  wilcox_result <- suppressWarnings(
    wilcox.test(
      score_difference,
      mu = 0,
      alternative = "two.sided",
      exact = FALSE,
      conf.int = FALSE
    )
  )

  # Diagnostic Q-Q plots against both fitted distributions.
  png(
    file.path(FIGURES, paste0("qq_model_comparison_", direction, ".png")),
    width = 1400,
    height = 650
  )

  old_par <- par(mfrow = c(1, 2))
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)

  # Gaussian Q-Q plot.
  probs <- ppoints(length(z))
  observed <- sort(z)

  theoretical_gaussian <- qnorm(
    probs,
    mean = gm,
    sd = gs
  )

  plot(
    theoretical_gaussian,
    observed,
    pch = 16,
    cex = 0.35,
    xlab = "Gaussian theoretical quantiles",
    ylab = "Observed standardized residuals",
    main = paste(direction, "- Gaussian Q-Q")
  )

  abline(
    a = gm,
    b = gs,
    col = "red",
    lwd = 2
  )

  # Student-t Q-Q plot.
  theoretical_student <- tl + ts * qt(probs, df = nu)

  plot(
    theoretical_student,
    observed,
    pch = 16,
    cex = 0.35,
    xlab = "Student-t theoretical quantiles",
    ylab = "Observed standardized residuals",
    main = paste(direction, "- Student-t Q-Q")
  )

  # Reference line through the central quantiles.
  q_obs <- quantile(z, c(0.25, 0.75), names = FALSE)
  q_theory <- tl + ts * qt(c(0.25, 0.75), df = nu)

  if (all(is.finite(q_theory)) && diff(q_theory) != 0) {
    slope <- diff(q_obs) / diff(q_theory)
    intercept <- q_obs[1] - slope * q_theory[1]
    abline(a = intercept, b = slope, col = "red", lwd = 2)
  }

  dev.off()
  on.exit(NULL, add = FALSE)

  data.frame(
    direction = direction,
    n = length(z),
    KS_gaussian_statistic = unname(ks_gaussian$statistic),
    KS_gaussian_p_value = ks_gaussian$p.value,
    KS_student_statistic = unname(ks_student$statistic),
    KS_student_p_value = ks_student$p.value,
    logLik_gaussian = ll_gaussian,
    logLik_student = ll_student,
    AIC_gaussian = aic_gaussian,
    AIC_student = aic_student,
    delta_AIC_student_minus_gaussian = aic_student - aic_gaussian,
    mean_log_score_difference = mean(score_difference),
    wilcoxon_statistic = unname(wilcox_result$statistic),
    wilcoxon_p_value = wilcox_result$p.value
  )
}

results <- do.call(
  rbind,
  lapply(names(directions), function(direction) {
    compare_direction(direction, directions[[direction]])
  })
)

write.csv(results, output_file, row.names = FALSE)

cat("\nGaussian vs Student-t comparison complete.\n")
cat("Results saved to:", output_file, "\n")
cat("Diagnostic plots saved to:", FIGURES, "\n\n")

print(results, row.names = FALSE)

cat(
  "\nInterpretation guide:\n",
  "- Lower AIC indicates the preferred fit among these two models.\n",
  "- Negative delta AIC favours Student-t; positive favours Gaussian.\n",
  "- Positive mean log-score difference favours Student-t.\n",
  "- KS p-values are descriptive because fitted parameters use the same data.\n",
  "- Wilcoxon p-values may be optimistic because residual rows are dependent.\n",
  sep = ""
)

