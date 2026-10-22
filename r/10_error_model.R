
# 10_error_model.R -- OWNER: Vishi
#
# Purpose: Fit age-dependent Gaussian and Student-t orbital error models.
# Reads: data/processed/residuals.csv
# Writes: data/processed/error_model_fit.csv,
#         data/processed/params.csv, figures/
# Status: DONE
# Run from repository root: Rscript r/10_error_model.R

source("r/common.R")

set.seed(SEED)

input_file <- file.path(PROCESSED, "residuals.csv")
error_model_file <- file.path(PROCESSED, "error_model_fit.csv")
output_file <- file.path(PROCESSED, "params.csv")

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

if (nrow(dat) < 30) {
  stop("Too few valid observations for model fitting.")
}

directions <- c(
  radial = "d_radial_km",
  along = "d_along_km",
  cross = "d_cross_km"
)

# Create age bins for robust empirical scale estimation.
age_breaks <- unique(as.numeric(
  quantile(
    dat$age_days,
    probs = seq(0, 1, length.out = 21),
    na.rm = TRUE
  )
))

if (length(age_breaks) < 4) {
  stop("Not enough distinct TLE ages to estimate scale growth.")
}

dat$age_bin <- cut(
  dat$age_days,
  breaks = age_breaks,
  include.lowest = TRUE
)

# Student-t negative log-likelihood.
# par[1] = location
# par[2] = log(nu)
# par[3] = log(scale)
student_nll <- function(par, z) {
  location <- par[1]
  nu <- exp(par[2])
  scale <- exp(par[3])

  if (!is.finite(nu) || nu <= 0 ||
      !is.finite(scale) || scale <= 0) {
    return(Inf)
  }

  value <- -sum(
    dt(
      (z - location) / scale,
      df = nu,
      log = TRUE
    ) - log(scale)
  )

  if (!is.finite(value)) {
    return(Inf)
  }

  value
}

fit_direction <- function(column, direction) {

  x <- dat[[column]]
  age <- dat$age_days

  # Estimate empirical residual scale within age bins.
  grouped <- split(
    data.frame(age = age, x = x),
    dat$age_bin
  )

  bins <- lapply(grouped, function(g) {

    s <- mad(
      g$x,
      center = median(g$x),
      constant = 1.4826
    )

    if (!is.finite(s) || s <= 0) {
      s <- sd(g$x)
    }

    if (!is.finite(s) || s <= 0) {
      return(NULL)
    }

    data.frame(
      age = median(g$age),
      sigma = s,
      n = nrow(g)
    )
  })

  bins <- Filter(Negate(is.null), bins)

  if (length(bins) == 0) {
    stop("No usable age bins for direction: ", direction)
  }

  bins <- do.call(rbind, bins)

  if (is.null(bins) || nrow(bins) < 3) {
    stop(
      "Insufficient usable age bins for direction: ",
      direction
    )
  }

  # Fit sigma(age) = a * age^b.
  growth_fit <- lm(
    log(sigma) ~ log(age),
    data = bins
  )

  a <- unname(exp(coef(growth_fit)[1]))
  b <- unname(coef(growth_fit)[2])

  if (!is.finite(a) || a <= 0 || !is.finite(b)) {
    stop("Invalid scale-growth fit for direction: ", direction)
  }

  sigma_age <- a * age^b

  if (any(!is.finite(sigma_age)) || any(sigma_age <= 0)) {
    stop("Invalid age-dependent scales for direction: ", direction)
  }

  # Standardize residuals using the fitted age-dependent scale.
  z <- x / sigma_age

  # Gaussian maximum-likelihood estimates.
  gaussian_mean <- mean(z)
  gaussian_sd <- sqrt(mean((z - gaussian_mean)^2))

  if (!is.finite(gaussian_sd) || gaussian_sd <= 0) {
    stop("Invalid Gaussian scale for direction: ", direction)
  }

  ll_gaussian <- sum(
    dnorm(
      z,
      mean = gaussian_mean,
      sd = gaussian_sd,
      log = TRUE
    )
  )

  # Student-t maximum-likelihood fit.
  # Estimate location, degrees of freedom and scale.
  initial_scale <- sd(z) * sqrt(4 / 6)

  if (!is.finite(initial_scale) || initial_scale <= 0) {
    initial_scale <- 1
  }

  opt <- optim(
    par = c(
      median(z),
      log(6 - 2),
      log(max(initial_scale, 1e-6))
    ),
    fn = student_nll,
    z = z,
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 2000, reltol = 1e-10)
  )

  if (opt$convergence != 0 || !is.finite(opt$value)) {
    warning(
      "Student-t optimization may not have converged for: ",
      direction
    )
  }

  student_location <- opt$par[1]
  nu <- exp(opt$par[2])
  t_scale <- exp(opt$par[3])
  ll_student <- -opt$value

  if (
    !is.finite(student_location) ||
    !is.finite(nu) ||
    !is.finite(t_scale) ||
    t_scale <= 0 ||
    !is.finite(ll_student)
  ) {
    stop("Invalid Student-t fit for direction: ", direction)
  }

  # AIC includes the number of fitted distribution parameters.
  # Gaussian: mean and standard deviation (2 parameters).
  # Student-t: location, scale and degrees of freedom (3 parameters).
  aic_gaussian <- 2 * 2 - 2 * ll_gaussian
  aic_student <- 2 * 3 - 2 * ll_student

  # Normal Q-Q diagnostic plot of standardized residuals.
  # Keep this as a normal Q-Q plot; the Student-t comparison
  # diagnostic belongs in r/11_gaussian_vs_learned_tests.R.
  png(
    file.path(FIGURES, paste0("qq_", direction, ".png")),
    width = 1000,
    height = 800
  )

  qqnorm(
    z,
    main = paste("Normal Q-Q Plot:", direction)
  )
  qqline(z, col = "red", lwd = 2)

  dev.off()

  # Empirical scale versus TLE age.
  png(
    file.path(FIGURES, paste0("scale_age_", direction, ".png")),
    width = 1000,
    height = 800
  )

  plot(
    bins$age,
    bins$sigma,
    log = "xy",
    xlab = "TLE age (days)",
    ylab = "Empirical residual scale (km)",
    main = paste("Residual Scale vs TLE Age:", direction),
    pch = 19
  )

  curve(
    a * x^b,
    add = TRUE,
    col = "red",
    lwd = 2
  )

  dev.off()

  data.frame(
    direction = direction,
    a = a,
    b = b,
    student_location = student_location,
    student_df = nu,
    student_scale = t_scale,
    gaussian_mean = gaussian_mean,
    gaussian_sd = gaussian_sd,
    logLik_gaussian = ll_gaussian,
    logLik_student = ll_student,
    AIC_gaussian = aic_gaussian,
    AIC_student = aic_student,
    n = length(z)
  )
}

# Fit all three orbital directions.
results <- do.call(
  rbind,
  lapply(names(directions), function(direction) {
    fit_direction(directions[[direction]], direction)
  })
)

# Save the wide model-fitting results for the R comparison script.
write.csv(
  results,
  error_model_file,
  row.names = FALSE
)

# Build the six-row parameter table expected by the Python pipeline.
gaussian_params <- data.frame(
  component = results$direction,
  model = "gaussian",
  a = results$a * results$gaussian_sd,
  b = results$b,
  nu = NA_real_,
  aic = results$AIC_gaussian,
  n_obs = results$n
)

student_params <- data.frame(
  component = results$direction,
  model = "student_t",
  a = results$a * results$student_scale,
  b = results$b,
  nu = results$student_df,
  aic = results$AIC_student,
  n_obs = results$n
)

params <- rbind(gaussian_params, student_params)

# Save parameters in the shared Python data-contract format.
write.csv(
  params,
  output_file,
  row.names = FALSE,
  na = ""
)

cat("\nError-model fitting complete.\n")
cat("Valid observations:", nrow(dat), "\n")
cat("Wide model results saved to:", error_model_file, "\n")
cat("Python contract parameters saved to:", output_file, "\n\n")

print(params, row.names = FALSE)

cat("\nDiagnostic plots saved to:", FIGURES, "\n")
