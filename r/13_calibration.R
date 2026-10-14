
# 13_calibration.R -- OWNER: Vishi
# Held-out calibration of 95% predictive intervals.
# Run from repository root:
# Rscript r/13_calibration.R

source("r/common.R")
set.seed(SEED)

input_file <- file.path(PROCESSED, "residuals.csv")
if (!file.exists(input_file)) {
  stop("Missing data/processed/residuals.csv")
}

dat <- read.csv(input_file, stringsAsFactors = FALSE)

required <- c(
  "age_days",
  "d_radial_km",
  "d_along_km",
  "d_cross_km"
)
missing <- setdiff(required, names(dat))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

dat <- dat[complete.cases(dat[, required]), required]
dat <- dat[
  is.finite(dat$age_days) & dat$age_days > 0 &
    is.finite(dat$d_radial_km) &
    is.finite(dat$d_along_km) &
    is.finite(dat$d_cross_km),
]

if (nrow(dat) < 100) {
  stop("Too few valid residuals for a train/test calibration check.")
}

directions <- c(
  radial = "d_radial_km",
  along = "d_along_km",
  cross = "d_cross_km"
)

# Random held-out split. Residuals sharing TLEs may be dependent,
# so these intervals are descriptive, not a fully independent test.
test_idx <- sample(seq_len(nrow(dat)), size = floor(0.20 * nrow(dat)))
train <- dat[-test_idx, , drop = FALSE]
test <- dat[test_idx, , drop = FALSE]

# Estimate age-dependent scale using training data only.
fit_scale <- function(age, x) {
  breaks <- unique(as.numeric(
    quantile(age, probs = seq(0, 1, length.out = 21), na.rm = TRUE)
  ))
  if (length(breaks) < 4) {
    stop("Insufficient distinct ages for scale fitting.")
  }

  bins <- cut(age, breaks = breaks, include.lowest = TRUE)
  grouped <- split(data.frame(age = age, x = x), bins)

  estimates <- lapply(grouped, function(g) {
    s <- mad(g$x, center = median(g$x), constant = 1.4826)
    if (!is.finite(s) || s <= 0) s <- sd(g$x)
    if (!is.finite(s) || s <= 0) return(NULL)

    data.frame(age = median(g$age), sigma = s)
  })
  estimates <- Filter(Negate(is.null), estimates)
  estimates <- do.call(rbind, estimates)

  if (is.null(estimates) || nrow(estimates) < 3) {
    stop("Could not estimate age-dependent scale.")
  }

  fit <- lm(log(sigma) ~ log(age), data = estimates)
  a <- unname(exp(coef(fit)[1]))
  b <- unname(coef(fit)[2])

  if (!is.finite(a) || !is.finite(b) || a <= 0) {
    stop("Invalid scale model.")
  }

  list(a = a, b = b)
}

# Fit Student-t location, degrees of freedom and scale by MLE.
student_nll <- function(par, z) {
  location <- par[1]
  nu <-  exp(par[2])
  scale <- exp(par[3])
  -sum(dt((z - location) / scale, df = nu, log = TRUE) - log(scale))
}

fit_distribution <- function(z) {
  gaussian_mean <- mean(z)
  gaussian_sd <- sqrt(mean((z - gaussian_mean)^2))

  opt <- optim(
    par = c(median(z), log(4), log(max(sd(z), 1e-6))),
    fn = student_nll,
    z = z,
    method = "BFGS",
    control = list(maxit = 2000, reltol = 1e-9)
  )

  if (opt$convergence != 0 || !is.finite(opt$value)) {
    warning("Student-t optimization may not have converged.")
  }

  list(
    gaussian_mean = gaussian_mean,
    gaussian_sd = gaussian_sd,
    student_location = opt$par[1],
    student_df = exp(opt$par[2]),
    student_scale = exp(opt$par[3])
  )
}

# Standardized two-sided central 95% predictive interval.
results <- list()

for (direction in names(directions)) {
  column <- directions[[direction]]

  scale_fit <- fit_scale(train$age_days, train[[column]])
  train_sigma <- scale_fit$a * train$age_days^scale_fit$b
  test_sigma <- scale_fit$a * test$age_days^scale_fit$b

  if (any(!is.finite(train_sigma)) || any(train_sigma <= 0) ||
      any(!is.finite(test_sigma)) || any(test_sigma <= 0)) {
    stop("Invalid predicted scale for direction: ", direction)
  }

  z_train <- train[[column]] / train_sigma
  z_test <- test[[column]] / test_sigma

  dist_fit <- fit_distribution(z_train)

  # Gaussian central 95% interval.
  g_lower <- qnorm(
    0.025,
    mean = dist_fit$gaussian_mean,
    sd = dist_fit$gaussian_sd
  )
  g_upper <- qnorm(
    0.975,
    mean = dist_fit$gaussian_mean,
    sd = dist_fit$gaussian_sd
  )

  # Student-t central 95% interval.
  t_lower <- dist_fit$student_location +
    dist_fit$student_scale * qt(0.025, df = dist_fit$student_df)
  t_upper <- dist_fit$student_location +
    dist_fit$student_scale * qt(0.975, df = dist_fit$student_df)

  gaussian_covered <- z_test >= g_lower & z_test <= g_upper
  student_covered <- z_test >= t_lower & z_test <= t_upper

  n_test <- length(z_test)
  g_count <- sum(gaussian_covered)
  t_count <- sum(student_covered)

  results[[direction]] <- data.frame(
    direction = direction,
    n_train = length(z_train),
    n_test = n_test,
    target_coverage = 0.95,
    gaussian_coverage = mean(gaussian_covered),
    gaussian_covered = g_count,
    gaussian_lower = g_lower,
    gaussian_upper = g_upper,
    student_t_coverage = mean(student_covered),
    student_t_covered = t_count,
    student_t_lower = t_lower,
    student_t_upper = t_upper,
    scale_a = scale_fit$a,
    scale_b = scale_fit$b,
    student_location = dist_fit$student_location,
    student_df = dist_fit$student_df,
    student_scale = dist_fit$student_scale
  )

  # Save held-out coverage plot for this direction.
  png(
    file.path(FIGURES, paste0("calibration_", direction, ".png")),
    width = 1000,
    height = 700
  )
  barplot(
    c(
      Target = 0.95,
      Gaussian = mean(gaussian_covered),
      Student_t = mean(student_covered)
    ),
    ylim = c(0, 1),
    ylab = "Held-out coverage",
    main = paste("95% Predictive Interval Coverage:", direction)
  )
  abline(h = 0.95, lty = 2, lwd = 2)
  dev.off()
}

calibration <- do.call(rbind, results)
rownames(calibration) <- NULL

output_file <- file.path(PROCESSED, "calibration_results.csv")
write.csv(calibration, output_file, row.names = FALSE)

cat("\nHeld-out calibration check complete.\n")
cat("Training observations:", nrow(train), "\n")
cat("Held-out observations:", nrow(test), "\n")
cat("Results saved to:", output_file, "\n")
cat("Calibration plots saved in:", FIGURES, "\n\n")
print(calibration, row.names = FALSE)

cat(
  "\nInterpretation: coverage near 0.95 is consistent with the nominal ",
  "95% interval on this held-out split. Coverage alone does not prove ",
  "that a model is well calibrated. Residuals sharing TLEs may be ",
  "dependent, and the random split may place related observations in ",
  "both sets. Consider an object-level or time-based split for stronger ",
  "validation.\n",
  sep = ""
)

