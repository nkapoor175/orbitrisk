
# 10_error_model.R -- OWNER: Vishi
# Fit age-dependent Student-t and Gaussian error models.
# Run from repository root:
# Rscript r/10_error_model.R

source("r/common.R")

set.seed(SEED)

input_file <- file.path(PROCESSED, "residuals.csv")
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

# Fit a power law to the residual scale versus age.
# Use robust empirical scale estimates within age bins.
age_breaks <- unique(as.numeric(
  quantile(dat$age_days, probs = seq(0, 1, length.out = 21))
))

if (length(age_breaks) < 4) {
  stop("Not enough distinct TLE ages to estimate scale growth.")
}

dat$age_bin <- cut(
  dat$age_days,
  breaks = age_breaks,
  include.lowest = TRUE
)

# Negative log-likelihood for a standardized Student-t distribution.
student_nll <- function(par, z) {
  nu <- 2 + exp(par[1])
  scale <- exp(par[2])
  -sum(dt(z / scale, df = nu, log = TRUE) - log(scale))
}

fit_direction <- function(column, direction) {
  x <- dat[[column]]
  age <- dat$age_days

  # Estimate empirical scale per age bin using MAD.
  grouped <- split(data.frame(age = age, x = x), dat$age_bin)

  bins <- lapply(grouped, function(g) {
    s <- mad(g$x, center = median(g$x), constant = 1.4826)
    if (!is.finite(s) || s <= 0) {
      s <- sd(g$x)
    }
    if (!is.finite(s) || s <= 0) return(NULL)

    data.frame(
      age = median(g$age),
      sigma = s,
      n = nrow(g)
    )
  })

  bins <- Filter(Negate(is.null), bins)
  bins <- do.call(rbind, bins)

  if (nrow(bins) < 3) {
    stop("Insufficient usable age bins for direction: ", direction)
  }

  growth_fit <- lm(log(sigma) ~ log(age), data = bins)
  a <- exp(coef(growth_fit)[1])
  b <- unname(coef(growth_fit)[2])

  sigma_age <- a * age^b
  z <- x / sigma_age

  # Gaussian fit: mean and standard deviation by MLE.
  gaussian_mean <- mean(z)
  gaussian_sd <- sqrt(mean((z - gaussian_mean)^2))

  if (!is.finite(gaussian_sd) || gaussian_sd <= 0) {
    stop("Invalid Gaussian scale for direction: ", direction)
  }

  ll_gaussian <- sum(dnorm(
    z, mean = gaussian_mean, sd = gaussian_sd, log = TRUE
  ))

  # Student-t fit: estimate degrees of freedom and scale by MLE.
  opt <- optim(
    par = c(log(6), log(1)),
    fn = student_nll,
    z = z,
    method = "BFGS",
    hessian = TRUE
  )

  nu <- 2 + exp(opt$par[1])
  t_scale <- exp(opt$par[2])
  ll_student <- -opt$value

  # AIC = 2*k - 2*log-likelihood.
  aic_gaussian <- 2 * 2 - 2 * ll_gaussian
  aic_student <- 2 * 2 - 2 * ll_student

  # Diagnostic plots.
  png(
    file.path(FIGURES, paste0("qq_", direction, ".png")),
    width = 1000, height = 800
  )
  qqnorm(z, main = paste("Normal Q-Q Plot:", direction))
  qqline(z, col = "red", lwd = 2)
  dev.off()

  png(
    file.path(FIGURES, paste0("scale_age_", direction, ".png")),
    width = 1000, height = 800
  )
  plot(
    bins$age, bins$sigma,
    log = "xy",
    xlab = "TLE age (days)",
    ylab = "Empirical residual scale (km)",
    main = paste("Residual Scale vs TLE Age:", direction),
    pch = 19
  )
  curve(a * x^b, add = TRUE, col = "red", lwd = 2)
  dev.off()

  data.frame(
    direction = direction,
    a = a,
    b = b,
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

results <- do.call(rbind, lapply(names(directions), function(direction) {
  fit_direction(directions[[direction]], direction)
}))

write.csv(results, output_file, row.names = FALSE)

cat("\nError-model fitting complete.\n")
cat("Valid observations:", nrow(dat), "\n")
cat("Parameters saved to:", output_file, "\n\n")
print(results)
cat("\nDiagnostic plots saved to:", FIGURES, "\n")



