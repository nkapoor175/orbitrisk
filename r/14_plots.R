# 14_plots.R -- OWNER: Vishi
# Purpose: Generate final report figures for error growth, model diagnostics,
#          Monte Carlo convergence, ROC analysis, and fleet-risk diagnostics.
# Reads: data/processed/residuals.csv, params.csv, pc_exact.csv,
#        conjunctions.csv, and available Monte Carlo outputs.
# Writes: figures/
# Status: Implemented; optional plots require suitable input data.
# Run from repository root: Rscript r/14_plots.R

source("r/common.R")

residual_file <- file.path(PROCESSED, "residuals.csv")
params_file <- file.path(PROCESSED, "params.csv")
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
  stop(
    "Missing columns in residuals.csv: ",
    paste(missing, collapse = ", ")
  )
}

dat <- dat[
  complete.cases(dat[, required]),
  required
]

dat <- dat[
  is.finite(dat$age_days) &
    dat$age_days > 0 &
    is.finite(dat$d_radial_km) &
    is.finite(dat$d_along_km) &
    is.finite(dat$d_cross_km),
]

if (nrow(dat) == 0) {
  stop("No valid residual observations remain.")
}

directions <- c(
  radial = "d_radial_km",
  along = "d_along_km",
  cross = "d_cross_km"
)

# 1. Empirical error growth against TLE age.
for (direction in names(directions)) {

  xcol <- directions[[direction]]

  breaks <- unique(as.numeric(
    quantile(
      dat$age_days,
      probs = seq(0, 1, length.out = 21),
      na.rm = TRUE
    )
  ))

  if (length(breaks) < 3) {
    warning("Skipping error-growth plot for ", direction,
            ": insufficient distinct age values.")
    next
  }

  age_bin <- cut(
    dat$age_days,
    breaks = breaks,
    include.lowest = TRUE
  )

  groups <- split(
    data.frame(
      age = dat$age_days,
      error = dat[[xcol]]
    ),
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

  if (length(summaries) == 0) next

  summaries <- do.call(rbind, summaries)

  summaries <- summaries[
    is.finite(summaries$age) &
      summaries$age > 0 &
      is.finite(summaries$median_abs_error) &
      summaries$median_abs_error > 0,
  ]

  if (nrow(summaries) < 2) {
    warning("Skipping error-growth plot for ", direction,
            ": insufficient positive scale values.")
    next
  }

  png(
    file.path(
      FIGURES,
      paste0("final_error_growth_", direction, ".png")
    ),
    width = 1000,
    height = 700
  )

  plot(
    summaries$age,
    summaries$median_abs_error,
    log = "xy",
    pch = 19,
    xlab = "TLE age (days)",
    ylab = "Median absolute residual (km)",
    main = paste("Error Growth with TLE Age:", direction)
  )

  dev.off()
}

# 2. Fitted-model Q-Q comparison plots.
# Use the standardized-residual model-comparison plots from r/11.
# Raw residuals pooled across TLE ages are deliberately not plotted here.

# 3. Monte Carlo convergence curves, when suitable output exists.
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

  if (length(sample_col) == 0 ||
      length(estimate_col) == 0) next

  sample_values <- suppressWarnings(
    as.numeric(d[[sample_col[1]]])
  )

  estimate_values <- suppressWarnings(
    as.numeric(d[[estimate_col[1]]])
  )

  valid <- is.finite(sample_values) &
    is.finite(estimate_values)

  if (sum(valid) < 2) next

  png(
    file.path(FIGURES, "final_monte_carlo_convergence.png"),
    width = 1000,
    height = 700
  )

  plot(
    sample_values[valid],
    estimate_values[valid],
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

# 4. ROC analysis using Parthvi's pc_exact.csv directly.
# Labels: Student-t Pc >= ALERT_THRESHOLD.
# Scores: Gaussian exact Pc.
# These are model-defined alert labels, not observed collision outcomes.

roc_found <- FALSE

if (file.exists(pc_file)) {

  pc <- read.csv(
    pc_file,
    stringsAsFactors = FALSE
  )

  cat("\npc_exact.csv columns:\n")
  print(names(pc))

  # Identify model, method, probability and conjunction ID columns.
  id_candidates <- c(
    "conj_id", "conjunction_id", "pair_id", "event_id", "id"
  )

  model_candidates <- c(
    "error_model", "model", "distribution"
  )

  method_candidates <- c(
    "method", "pc_method", "estimator"
  )

  probability_candidates <- c(
    "pc", "Pc", "collision_probability",
    "probability", "pc_estimate"
  )

  id_col <- intersect(id_candidates, names(pc))
  model_col <- intersect(model_candidates, names(pc))
  method_col <- intersect(method_candidates, names(pc))
  prob_col <- intersect(probability_candidates, names(pc))

  if (length(id_col) > 0 &&
      length(model_col) > 0 &&
      length(prob_col) > 0) {

    id_col <- id_col[1]
    model_col <- model_col[1]
    prob_col <- prob_col[1]

    pc[[prob_col]] <- suppressWarnings(
      as.numeric(pc[[prob_col]])
    )

    model_names <- tolower(trimws(as.character(pc[[model_col]])))

    student_rows <- model_names %in% c(
      "student_t", "student-t", "studentt", "student"
    )

    gaussian_rows <- model_names %in% c(
      "gaussian", "normal"
    )

    student <- pc[
      student_rows,
      c(id_col, prob_col),
      drop = FALSE
    ]

    gaussian <- pc[
      gaussian_rows,
      c(id_col, prob_col),
      drop = FALSE
    ]

    names(student) <- c("join_id", "student_pc")
    names(gaussian) <- c("join_id", "gaussian_pc")

    # If the file has multiple methods, use a consistent method
    # when a recognized method column is available.
    if (length(method_col) > 0) {

      method_col <- method_col[1]

      student_method <- tolower(
        trimws(as.character(pc[[method_col]][student_rows]))
      )

      gaussian_method <- tolower(
        trimws(as.character(pc[[method_col]][gaussian_rows]))
      )

      common_methods <- intersect(
        unique(student_method),
        unique(gaussian_method)
      )

      if (length(common_methods) > 0) {
        chosen_method <- common_methods[1]

        student <- student[
          student_method == chosen_method,
          ,
          drop = FALSE
        ]

        gaussian <- gaussian[
          gaussian_method == chosen_method,
          ,
          drop = FALSE
        ]
      }
    }

    student <- student[
      is.finite(student$student_pc) &
        student$student_pc >= 0 &
        student$student_pc <= 1,
    ]

    gaussian <- gaussian[
      is.finite(gaussian$gaussian_pc) &
        gaussian$gaussian_pc >= 0 &
        gaussian$gaussian_pc <= 1,
    ]

    # Keep one probability per conjunction and model.
    # If duplicate rows remain, do not silently inflate the sample.
    if (anyDuplicated(student$join_id) ||
        anyDuplicated(gaussian$join_id)) {

      warning(
        "ROC skipped: duplicate conjunction IDs remain after ",
        "filtering pc_exact.csv to the selected model/method."
      )

    } else {

      matched <- merge(
        student,
        gaussian,
        by = "join_id",
        all = FALSE
      )

      if (nrow(matched) > 0) {

        matched$label <- as.integer(
          matched$student_pc >= ALERT_THRESHOLD
        )

        matched$score <- matched$gaussian_pc

        valid <- is.finite(matched$label) &
          is.finite(matched$score)

        matched <- matched[valid, ]

        if (length(unique(matched$label)) == 2) {

          y <- matched$label
          score <- matched$score

          thresholds <- c(
            Inf,
            sort(unique(score), decreasing = TRUE),
            -Inf
          )

          tpr <- vapply(thresholds, function(th) {
            sum(score[y == 1] >= th) / sum(y == 1)
          }, numeric(1))

          fpr <- vapply(thresholds, function(th) {
            sum(score[y == 0] >= th) / sum(y == 0)
          }, numeric(1))

          ord <- order(fpr, tpr)

          roc_auc <- sum(
            diff(fpr[ord]) *
              (head(tpr[ord], -1) + tail(tpr[ord], -1)) / 2
          )

          # Report false-alarm and missed-alert rates at the
          # project's shared alert threshold.
          gaussian_alert <- score >= ALERT_THRESHOLD

          false_alarm_rate <- sum(
            gaussian_alert & y == 0
          ) / sum(y == 0)

          missed_alert_rate <- sum(
            !gaussian_alert & y == 1
          ) / sum(y == 1)

          png(
            file.path(FIGURES, "final_roc_curve.png"),
            width = 1000,
            height = 700
          )

          plot(
            fpr[ord],
            tpr[ord],
            type = "l",
            lwd = 2,
            xlim = c(0, 1),
            ylim = c(0, 1),
            xlab = "False positive rate",
            ylab = "True positive rate",
            main = paste0(
              "Gaussian Pc vs Student-t Alerts (AUC = ",
              round(roc_auc, 3),
              ")"
            )
          )

          abline(0, 1, lty = 2)
          dev.off()

          roc_summary <- data.frame(
            n_matched = nrow(matched),
            alert_threshold = ALERT_THRESHOLD,
            false_alarm_rate = false_alarm_rate,
            missed_alert_rate = missed_alert_rate,
            auc = roc_auc
          )

          write.csv(
            roc_summary,
            file.path(PROCESSED, "roc_summary.csv"),
            row.names = FALSE
          )

          cat("\nROC analysis complete.\n")
          print(roc_summary, row.names = FALSE)

          roc_found <- TRUE
        }
      }
    }

  } else {
    cat(
      "Skipping ROC curve: pc_exact.csv lacks a recognized model, ",
      "conjunction ID, or probability column.\n",
      sep = ""
    )
  }

} else {
  cat("Skipping ROC curve: pc_exact.csv was not found.\n")
}

if (!roc_found) {
  cat(
    "ROC analysis skipped: both alert classes and matching model ",
    "probabilities are required.\n",
    sep = ""
  )
}

# 5. Poisson gap diagnostic, only with enough distinct high-risk times.
poisson_found <- FALSE

if (file.exists(conjunction_file) && file.exists(pc_file)) {

  con <- read.csv(
    conjunction_file,
    stringsAsFactors = FALSE
  )

  pc <- read.csv(
    pc_file,
    stringsAsFactors = FALSE
  )

  id_con <- intersect(
    c("conj_id", "conjunction_id", "pair_id", "event_id", "id"),
    names(con)
  )

  id_pc <- intersect(
    c("conj_id", "conjunction_id", "pair_id", "event_id", "id"),
    names(pc)
  )

  time_col <- intersect(
    c("tca_utc", "tca", "event_time", "time_utc"),
    names(con)
  )

  model_col <- intersect(
    c("error_model", "model", "distribution"),
    names(pc)
  )

  prob_col <- intersect(
    c("pc", "Pc", "collision_probability", "probability", "pc_estimate"),
    names(pc)
  )

  if (length(id_con) > 0 &&
      length(id_pc) > 0 &&
      length(time_col) > 0 &&
      length(prob_col) > 0) {

    id_con <- id_con[1]
    id_pc <- id_pc[1]
    time_col <- time_col[1]
    prob_col <- prob_col[1]

    if (length(model_col) > 0) {

      model_col <- model_col[1]

      model_names <- tolower(trimws(as.character(pc[[model_col]])))

      selected <- model_names %in% c(
        "student_t", "student-t", "studentt", "student"
      )

      if (any(selected)) {
        pc <- pc[selected, , drop = FALSE]
      }
    }

    pc[[prob_col]] <- suppressWarnings(
      as.numeric(pc[[prob_col]])
    )

    pc_subset <- pc[, c(id_pc, prob_col), drop = FALSE]
    names(pc_subset) <- c("join_id", "Pc")

    if (!anyDuplicated(pc_subset$join_id)) {

      con_subset <- con[, c(id_con, time_col), drop = FALSE]
      names(con_subset) <- c("join_id", "event_time")

      matched <- merge(
        con_subset,
        pc_subset,
        by = "join_id",
        all = FALSE
      )

      matched$event_time <- as.POSIXct(
        matched$event_time,
        tz = "UTC",
        tryFormats = c(
          "%Y-%m-%dT%H:%M:%S",
          "%Y-%m-%d %H:%M:%S",
          "%Y-%m-%dT%H:%M:%OS"
        )
      )

      matched <- matched[
        is.finite(matched$Pc) &
          matched$Pc >= ALERT_THRESHOLD &
          !is.na(matched$event_time),
      ]

      event_times <- sort(unique(matched$event_time))

      if (length(event_times) >= 3) {

        gaps_days <- as.numeric(
          diff(event_times),
          units = "days"
        )

        gaps_days <- gaps_days[
          is.finite(gaps_days) & gaps_days > 0
        ]

        if (length(gaps_days) >= 2) {

          rate_per_day <- 1 / mean(gaps_days)

          png(
            file.path(FIGURES, "final_poisson_gap_diagnostic.png"),
            width = 1000,
            height = 700
          )

          hist(
            gaps_days,
            probability = TRUE,
            xlab = "Gap between high-risk event times (days)",
            main = "High-Risk Event Gap Diagnostic"
          )

          curve(
            rate_per_day * exp(-rate_per_day * x),
            from = 0,
            to = max(gaps_days),
            add = TRUE,
            lwd = 2
          )

          dev.off()
          poisson_found <- TRUE
        }
      }
    }
  }
}

if (!poisson_found) {
  cat(
    "Skipping Poisson gap diagnostic: insufficient distinct high-risk ",
    "event times or required columns are missing.\n",
    sep = ""
  )
}

cat("\nFinal plotting script completed.\n")
cat("Figures saved in:", FIGURES, "\n")
cat("Generated: error-growth plots.\n")
cat("Model-comparison Q-Q plots are produced by r/11.\n")
cat("Optional plots are generated only when suitable input data exist.\n")