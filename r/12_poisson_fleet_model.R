# 12_poisson_fleet_model.R -- OWNER: Vishi
# Purpose: Fleet-level Poisson counts and exponential inter-event gaps.
# Reads: data/processed/conjunctions.csv, probability results
# Writes: data/processed/fleet_risk/, figures/fleet_inter_event_gaps.png
# Status: Implemented
# Run from repository root: Rscript r/12_poisson_fleet_model.R

source("r/common.R")

set.seed(SEED)

conjunction_file <- file.path(PROCESSED, "conjunctions.csv")
pc_file <- file.path(PROCESSED, "pc_results.csv")

if (!file.exists(pc_file)) {
  alternative_pc_file <- file.path(PROCESSED, "pc_exact.csv")
  if (file.exists(alternative_pc_file)) {
    pc_file <- alternative_pc_file
  }
}

out_dir <- file.path(PROCESSED, "fleet_risk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(conjunction_file)) {
  stop("Missing data/processed/conjunctions.csv")
}
if (!file.exists(pc_file)) {
  stop("Missing pc_results.csv and pc_exact.csv in data/processed/")
}

conj <- read.csv(conjunction_file, stringsAsFactors = FALSE)
pc <- read.csv(pc_file, stringsAsFactors = FALSE)

# Identify the collision-probability column.
pc_candidates <- c(
  "Pc", "pc", "Pc_exact", "pc_exact",
  "collision_probability", "probability"
)
pc_col <- intersect(pc_candidates, names(pc))

if (length(pc_col) == 0) {
  stop("No recognized collision-probability column in ", basename(pc_file))
}
pc_col <- pc_col[1]

# Filter to one method and error model to avoid duplicate conjunctions.
if (all(c("method", "error_model") %in% names(pc))) {
  combinations <- unique(pc[, c("method", "error_model"), drop = FALSE])

  exact_match <- pc$method %in% c("exact", "gaussian_exact", "exact_gaussian") &
    pc$error_model %in% c("gaussian", "Gaussian")

  if (any(exact_match)) {
    pc <- pc[exact_match, , drop = FALSE]
  } else {
    combinations <- combinations[
      order(combinations$method, combinations$error_model),
      ,
      drop = FALSE
    ]

    chosen_method <- combinations$method[1]
    chosen_model <- combinations$error_model[1]

    pc <- pc[
      pc$method == chosen_method &
        pc$error_model == chosen_model,
      ,
      drop = FALSE
    ]

    warning(
      "Exact Gaussian results were not identified. Using: ",
      chosen_method, " / ", chosen_model
    )
  }
} else if ("method" %in% names(pc)) {
  chosen_method <- sort(unique(pc$method))[1]
  pc <- pc[pc$method == chosen_method, , drop = FALSE]
  warning("Filtered probabilities to method: ", chosen_method)
}

# Match conjunctions using a shared event ID.
id_candidates <- c("conj_id", "conjunction_id", "pair_id", "event_id", "id")
id_conj <- intersect(id_candidates, names(conj))
id_pc <- intersect(id_candidates, names(pc))

if (length(id_conj) > 0 && length(id_pc) > 0) {
  names(conj)[names(conj) == id_conj[1]] <- "join_id"
  names(pc)[names(pc) == id_pc[1]] <- "join_id"

  if (anyDuplicated(pc$join_id)) {
    stop(
      "Duplicate conjunction IDs remain after filtering. ",
      "Check method/error_model selection."
    )
  }

  probability_data <- pc[, c("join_id", pc_col), drop = FALSE]
  names(probability_data)[2] <- "Pc"

  dat <- merge(
    conj,
    probability_data,
    by = "join_id",
    all = FALSE,
    sort = FALSE
  )
} else {
  stop(
    "Cannot safely match conjunctions to probabilities: ",
    "a common event ID is required. Row-order matching is disabled."
  )
}

dat$Pc <- suppressWarnings(as.numeric(dat$Pc))
dat <- dat[
  is.finite(dat$Pc) & dat$Pc >= 0 & dat$Pc <= 1,
  ,
  drop = FALSE
]

if (nrow(dat) == 0) {
  stop("No valid collision probabilities after matching.")
}

# Parse conjunction time.
time_candidates <- c(
  "tca_utc", "tca", "event_time", "epoch_utc", "time_utc"
)
time_col <- intersect(time_candidates, names(dat))

if (length(time_col) > 0) {
  dat$event_time <- as.POSIXct(
    dat[[time_col[1]]],
    tz = "UTC",
    tryFormats = c(
      "%Y-%m-%dT%H:%M:%OS",
      "%Y-%m-%d %H:%M:%OS",
      "%Y-%m-%dT%H:%M:%S",
      "%Y-%m-%d %H:%M:%S"
    )
  )
} else {
  dat$event_time <- as.POSIXct(
    rep(NA_real_, nrow(dat)),
    origin = "1970-01-01",
    tz = "UTC"
  )
  warning("No recognized event-time column found.")
}

dat$high_risk <- dat$Pc >= ALERT_THRESHOLD

# Expected collision count under the rare-event approximation.
expected_collision_count <- sum(dat$Pc)
prob_at_least_one_collision <- 1 - exp(-sum(dat$Pc))

cat(
  "Probability of at least one fleet collision:",
  prob_at_least_one_collision,
  "\n"
)
observed_high_risk <- sum(dat$high_risk)


summary_rows <- list(
  data.frame(
    metric = c(
      "matched_valid_conjunctions",
      "high_risk_conjunctions",
      "alert_threshold",
      "sum_collision_probabilities",
      "probability_at_least_one_fleet_collision"
    ),
    value = c(
      nrow(dat),
      observed_high_risk,
      ALERT_THRESHOLD,
      expected_collision_count,
      prob_at_least_one_collision
    )
  )
)


# Time-window analysis using valid timestamps.
valid_times <- dat$event_time[!is.na(dat$event_time)]

if (length(valid_times) >= 2) {
  window_start <- min(valid_times)
  window_end <- max(valid_times)

  window_days <- as.numeric(
    difftime(window_end, window_start, units = "days")
  )

  if (is.finite(window_days) && window_days > 0) {
    observed_rate <- observed_high_risk / window_days

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      metric = c(
        "observed_time_window_days",
        "observed_high_risk_rate_per_day",
        "expected_collision_count_from_probability_sum"
      ),
      value = c(
        window_days,
        observed_rate,
        expected_collision_count
      )
    )

    # Inter-event gaps among high-risk conjunctions.
    high_times <- sort(
      dat$event_time[dat$high_risk & !is.na(dat$event_time)]
    )

    gaps_days <- as.numeric(diff(high_times), units = "days")
    gaps_days <- gaps_days[is.finite(gaps_days) & gaps_days > 0]

    if (length(gaps_days) >= 2) {
      mean_gap <- mean(gaps_days)
      lambda_gap <- 1 / mean_gap

      gap_summary <- data.frame(
        n_gaps = length(gaps_days),
        mean_gap_days = mean_gap,
        median_gap_days = median(gaps_days),
        exponential_rate_per_day = lambda_gap
      )

      write.csv(
        gap_summary,
        file.path(out_dir, "exponential_gap_summary.csv"),
        row.names = FALSE
      )

      png(
        file.path(FIGURES, "fleet_inter_event_gaps.png"),
        width = 1000,
        height = 750
      )

      hist(
        gaps_days,
        breaks = "FD",
        main = "Gaps Between High-Risk Conjunctions",
        xlab = "Inter-event gap (days)",
        ylab = "Frequency"
      )

      dev.off()
    } else {
      warning("Not enough distinct high-risk event times to fit gap model.")
    }
  } else {
    warning("Event timestamps do not span a positive time window.")
  }
} else {
  warning("Too few valid timestamps for time-window analysis.")
}

summary_table <- do.call(rbind, summary_rows)

write.csv(
  summary_table,
  file.path(out_dir, "fleet_risk_summary.csv"),
  row.names = FALSE
)

write.csv(
  dat,
  file.path(out_dir, "matched_conjunction_probabilities.csv"),
  row.names = FALSE
)

cat("\nFleet-level analysis completed.\n")
cat("Input conjunctions:", basename(conjunction_file), "\n")
cat("Input probabilities:", basename(pc_file), "\n")
cat("Matched valid events:", nrow(dat), "\n")
cat("High-risk events (Pc >= alert threshold):", observed_high_risk, "\n")
cat("Expected collision count (sum of Pc):", expected_collision_count, "\n")
cat("Results saved in:", out_dir, "\n\n")

print(summary_table, row.names = FALSE)