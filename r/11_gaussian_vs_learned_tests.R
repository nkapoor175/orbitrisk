
# 11_gaussian_vs_learned_tests.R -- OWNER: Vishi
#
# Purpose: Compare Gaussian and Student-t error models and estimate
#          Type I error, Type II error and statistical power.
# Reads: data/processed/residuals.csv,
#        data/processed/error_model_fit.csv
# Writes: data/processed/model_comparison_tests.csv,
#         data/processed/type_error_sim.csv, figures/
# Status: Implemented; simulation results require validation.
# Run from repository root: Rscript r/11_gaussian_vs_learned_tests.R

source("r/common.R")

set.seed(SEED)

residual_file <- file.path(PROCESSED, "residuals.csv")
params_file <- file.path(PROCESSED, "error_model_fit.csv")
output_file <- file.path(PROCESSED, "model_comparison_tests.csv")
simulation_file <- file.path(PROCESSED, "type_error_sim.csv")

if (!file.exists(residual_file)) {
  stop("Missing data/processed/residuals.csv.")
}

if (!file.exists(params_file)) {
  stop("Missing data/processed/error_model_fit.csv. Run r/10_error_model.R first.")
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
  stop("Missing model columns: ", paste(missing_params, collapse = ", "))
}

directions <- c(
  radial = "d_radial_km",
  along = "d_along_km",
  cross = "d_cross_km"
)

# Fit Gaussian and Student-t distributions to a sample.
# Student-t parameters: location, degrees of freedom and scale.
fit_distributions <- function(z) {

  z <- z[is.finite(z)]

  if (length(z) < 10 || !is.finite(sd(z)) || sd(z) <= 0) {
    return(NULL)
  }

  gm <- mean(z)
  gs <- sqrt(mean((z - gm)^2))

  if (!is.finite(gs) || gs <= 0) return(NULL)

  gaussian_loglik <- sum(dnorm(z, gm, gs, log = TRUE))

  student_nll <- function(par) {
    location <- par[1]
    nu <- 2 + exp(par[2])
    scale <- exp(par[3])

    if (!is.finite(nu) || !is.finite(scale) || scale <= 0) {
      return(Inf)
    }

    ll <- sum(
      dt((z - location) / scale, df = nu, log = TRUE) -
        log(scale)
    )

    if (!is.finite(ll)) return(Inf)
    -ll
  }

  opt <- tryCatch(
    optim(
      par = c(median(z), log(4), log(max(sd(z) * sqrt(4 / 6), 1e-6))),
      fn = student_nll,
      method = "BFGS",
      control = list(maxit = 500, reltol = 1e-8)
    ),
    error = function(e) NULL
  )

  if (is.null(opt) || !is.finite(opt$value)) return(NULL)

  list(
    gaussian_mean = gm,
    gaussian_sd = gs,
    gaussian_loglik = gaussian_loglik,
    student_location = opt$par[1],
    student_df = 2 + exp(opt$par[2]),
    student_scale = exp(opt$par[3]),
    student_loglik = -opt$value
  )
}

# Compare observed standardized residuals for one direction.
compare_direction <- function(direction, column) {

  p <- params[params$direction == direction, , drop = FALSE]

  if (nrow(p) != 1) {
    stop("Expected exactly one fitted parameter row for: ", direction)
  }

  d <- dat[
    is.finite(dat$age_days) &
      dat$age_days > 0 &
      is.finite(dat[[column]]),
    c("age_days", column),
    drop = FALSE
  ]

  # Use the same robust age-scale function used in r/10.
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
    stop("Invalid fitted parameters for: ", direction)
  }

  log_gaussian <- dnorm(z, mean = gm, sd = gs, log = TRUE)
  log_student <- dt(
    (z - tl) / ts,
    df = nu,
    log = TRUE
  ) - log(ts)

  ll_gaussian <- sum(log_gaussian)
  ll_student <- sum(log_student)

  aic_gaussian <- 2 * 2 - 2 * ll_gaussian
  aic_student <- 2 * 3 - 2 * ll_student

  ks_gaussian <- suppressWarnings(
    ks.test(z, "pnorm", mean = gm, sd = gs)
  )

  student_cdf <- function(q) {
    pt((q - tl) / ts, df = nu)
  }

  ks_student <- suppressWarnings(
    ks.test(z, student_cdf)
  )

  score_difference <- log_student - log_gaussian

  wilcox_result <- suppressWarnings(
    wilcox.test(
      score_difference,
      mu = 0,
      alternative = "two.sided",
      exact = FALSE
    )
  )

  # Model-comparison Q-Q plots for the report.
  png(
    file.path(FIGURES, paste0("qq_model_comparison_", direction, ".png")),
    width = 1400,
    height = 650
  )

  old_par <- par(mfrow = c(1, 2))

  probs <- ppoints(length(z))
  observed <- sort(z)

  plot(
    qnorm(probs, mean = gm, sd = gs),
    observed,
    pch = 16,
    cex = 0.35,
    xlab = "Gaussian theoretical quantiles",
    ylab = "Observed standardized residuals",
    main = paste(direction, "- Gaussian Q-Q")
  )
  abline(a = gm, b = gs, col = "red", lwd = 2)

  plot(
    tl + ts * qt(probs, df = nu),
    observed,
    pch = 16,
    cex = 0.35,
    xlab = "Student-t theoretical quantiles",
    ylab = "Observed standardized residuals",
    main = paste(direction, "- Student-t Q-Q")
  )

  q_obs <- quantile(z, c(0.25, 0.75), names = FALSE)
  q_theory <- tl + ts * qt(c(0.25, 0.75), df = nu)

  if (all(is.finite(q_theory)) && diff(q_theory) != 0) {
    slope <- diff(q_obs) / diff(q_theory)
    intercept <- q_obs[1] - slope * q_theory[1]
    abline(a = intercept, b = slope, col = "red", lwd = 2)
  }

  par(old_par)
  dev.off()

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
print(results, row.names = FALSE)

# ------------------------------------------------------------------
# Type I / Type II error simulation
#
# H0: Gaussian-generated data.
# H1: Student-t-generated data.
# Each repetition refits both models and applies a paired
# Wilcoxon test to their per-observation log-score differences.
# ------------------------------------------------------------------

sample_sizes <- c(50, 100, 200, 500, 1000)
N_REPS <- 1000
ALPHA <- 0.05

simulate_one <- function(n, truth, fitted) {

  if (truth == "gaussian") {
    z <- rnorm(
      n,
      mean = fitted$gaussian_mean,
      sd = fitted$gaussian_sd
    )
  } else {
    z <- fitted$student_location +
      fitted$student_scale *
      rt(n, df = fitted$student_df)
  }

  fitted_sample <- fit_distributions(z)

  if (is.null(fitted_sample)) return(NA_real_)

  score_g <- dnorm(
    z,
    fitted_sample$gaussian_mean,
    fitted_sample$gaussian_sd,
    log = TRUE
  )

  score_t <- dt(
    (z - fitted_sample$student_location) /
      fitted_sample$student_scale,
    df = fitted_sample$student_df,
    log = TRUE
  ) - log(fitted_sample$student_scale)

  difference <- score_t - score_g

  test <- tryCatch(
    suppressWarnings(
      wilcox.test(
        difference,
        mu = 0,
        alternative = "two.sided",
        exact = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(test) || !is.finite(test$p.value)) {
    return(NA_real_)
  }

  as.numeric(test$p.value < ALPHA)
}

simulation_rows <- list()
row_index <- 1

for (direction in names(directions)) {

  p <- params[params$direction == direction, , drop = FALSE]

  fitted <- list(
    gaussian_mean = p$gaussian_mean,
    gaussian_sd = p$gaussian_sd,
    student_location = p$student_location,
    student_df = p$student_df,
    student_scale = p$student_scale
  )

  cat("\nSimulating Type I/II errors for:", direction, "\n")

  for (n in sample_sizes) {

    for (truth in c("gaussian", "student_t")) {

      cat(
        "  n =", n,
        "| truth =", truth,
        "| repetitions =", N_REPS, "\n"
      )

      rejected <- replicate(
        N_REPS,
        simulate_one(n, truth, fitted)
      )

      valid <- rejected[is.finite(rejected)]

      if (length(valid) == 0) {
        rejection_rate <- NA_real_
      } else {
        rejection_rate <- mean(valid)
      }

      simulation_rows[[row_index]] <- data.frame(
        direction = direction,
        n = n,
        truth = truth,
        reject_rate = rejection_rate,
        valid_repetitions = length(valid)
      )

      row_index <- row_index + 1
    }
  }
}

simulation_results <- do.call(rbind, simulation_rows)

write.csv(
  simulation_results,
  simulation_file,
  row.names = FALSE
)

# Plot rejection rate by sample size under each generating model.
# Gaussian truth estimates Type I error; Student-t truth estimates power.
for (direction in names(directions)) {

  sim <- simulation_results[
    simulation_results$direction == direction,
    ,
    drop = FALSE
  ]

  png(
    file.path(FIGURES, paste0("type_error_power_", direction, ".png")),
    width = 1000,
    height = 750
  )

  plot(
    NA,
    xlim = range(sample_sizes),
    ylim = c(0, 1),
    xlab = "Sample size",
    ylab = "Rejection rate",
    main = paste("Type I Error and Power:", direction),
    xaxt = "n"
  )

  axis(1, at = sample_sizes, labels = sample_sizes)
  abline(h = ALPHA, lty = 2)

  for (truth in c("gaussian", "student_t")) {

    s <- sim[sim$truth == truth, , drop = FALSE]
    s <- s[order(s$n), ]

    lines(
      s$n,
      s$reject_rate,
      type = "b",
      pch = if (truth == "gaussian") 16 else 17,
      lty = if (truth == "gaussian") 2 else 1
    )
  }

  legend(
    "bottomright",
    legend = c("Gaussian truth (Type I rejection)", "Student-t truth (power)"),
    pch = c(16, 17),
    lty = c(2, 1),
    bty = "n"
  )

  dev.off()
}

cat("\nType I/II simulation complete.\n")
cat("Simulation results saved to:", simulation_file, "\n")
cat("Simulation plots saved to:", FIGURES, "\n")
cat(
  "\nInterpretation:\n",
  "- Type I error = rejection rate when Gaussian data are generated.\n",
  "- Type II error = 1 - rejection rate when Student-t data are generated.\n",
  "- Power = rejection rate when Student-t data are generated.\n",
  "- The nominal significance level is ", ALPHA, ".\n",
  "- Residual dependence and in-sample model comparison are limitations.\n",
  sep = ""
)
