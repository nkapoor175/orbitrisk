
# 12_poisson_fleet_model.R -- OWNER: Vishi
# Fleet-level Poisson counts and exponential inter-event gaps.
# Run from repository root:
# Rscript r/12_poisson_fleet_model.R

source("r/common.R")

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

cat("Conjunction columns:\n")
print(names(conj))
cat("\nProbability columns:\n")
print(names(pc))

# Identify a probability column from common project naming conventions.
pc_candidates <- c(
  "Pc", "pc", "Pc_exact", "pc_exact",
  "collision_probability", "probability"
)
pc_col <- pc_candidates[pc_candidates %in% names(pc)][1]

if (is.na(pc_col)) {
  stop(
    "Could not identify collision-probability column in ",
    basename(pc_file),
    ". Inspect the column names printed above."
  )
}

# Identify common conjunction identifiers for joining.
id_candidates <- c(
  "conj_id", "conjunction_id", "pair_id", "event_id", "id"
)
id_col_conj <- id_candidates[id_candidates %in% names(conj)][1]
id_col_pc <- id_candidates[id_candidates %in% names(pc)][1]

if (!is.na(id_col_conj) && !is.na(id_col_pc)) {
  names(conj)[names(conj) == id_col_conj] <- "join_id"
  names(pc)[names(pc) == id_col_pc] <- "join_id"

  dat <- merge(conj, pc[, c("join_id", pc_col)],
               by = "join_id", all = FALSE)
  names(dat)[names(dat) == pc_col] <- "Pc"
} else if (nrow(conj) == nrow(pc)) {
  # Row-order matching is allowed only when row counts agree.
  warning(
    "No shared event ID found. Matching by row order because ",
    "the two files have equal row counts; verify this assumption."
  )
  dat <- conj
  dat$Pc <- pc[[pc_col]]
} else {
  stop(
    "Cannot safely match conjunctions to probabilities: ",
    "no common event ID and row counts differ."
  )
}

dat$Pc <- suppressWarnings(as.numeric(dat$Pc))
dat <- dat[is.finite(dat$Pc) & dat$Pc >= 0 & dat$Pc <= 1, ]

if (nrow(dat) == 0) {
  stop("No valid collision probabilities after matching the input files.")
}

# Find a usable event-time column.
time_candidates <- c(
  "tca_utc", "tca", "event_time", "epoch_utc", "time_utc"
)
time_col <- time_candidates[time_candidates %in% names(dat)][1]

if (!is.na(time_col)) {
  dat$event_time <- as.POSIXct(
    dat[[time_col]],
    tz = "UTC",
    tryFormats = c(
      "%Y-%m-%dT%H:%M:%S",
      "%Y-%m-%d %H:%M:%S",
      "%Y-%m-%dT%H:%M:%OS"
    )
  )
} else {
  dat$event_time <- as.POSIXct(rep(NA, nrow(dat)), origin = "1970-01-01", tz = "UTC")
  warning(
    "No recognized event-time column found. ",
    "Poisson counts will use observed rows, not time-based rates."
  )
}

# High-risk events use the project's shared alert threshold.
dat$high_risk <- dat$Pc >= ALERT_THRESHOLD

summary_rows <- list()

# Probability sum is an expected event count only under the relevant
# rare-event approximation and if the probabilities describe the same
# specified observation window.
expected_count <- sum(dat$Pc)
observed_high_risk <- sum(dat$high_risk)

summary_rows[[1]] <- data.frame(
  metric = c(
    "matched_valid_conjunctions",
    "high_risk_conjunctions",
    "alert_threshold",
    "sum_collision_probabilities"
  ),
  value = c(
    nrow(dat),
    observed_high_risk,
    ALERT_THRESHOLD,
    expected_count
  ),
  stringsAsFactors = FALSE
)

# Estimate a count rate only when actual timestamps span a positive window.
if (all(!is.na(dat$event_time))) {
  dat <- dat[order(dat$event_time), ]
  window_days <- as.numeric(
    difftime(max(dat$event_time), min(dat$event_time), units = "days")
  )

  if (is.finite(window_days) && window_days > 0) {
    rate <- observed_high_risk / window_days
    poisson_mean <- rate * window_days

    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      metric = c(
        "observed_time_window_days",
        "observed_high_risk_rate_per_day",
        "poisson_expected_high_risk_count"
      ),
      value = c(window_days, rate, poisson_mean),
      stringsAsFactors = FALSE
    )

    # Gaps between consecutive high-risk event timestamps.
    high_times <- sort(dat$event_time[dat$high_risk])
    gaps_days <- as.numeric(diff(high_times), units = "days")
    gaps_days <- gaps_days[is.finite(gaps_days) & gaps_days >= 0]

    if (length(gaps_days) >= 2 && mean(gaps_days) > 0) {
      lambda_gap <- 1 / mean(gaps_days)

      gap_summary <- data.frame(
        n_gaps = length(gaps_days),
        mean_gap_days = mean(gaps_days),
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
cat("Results saved in:", out_dir, "\n")
cat("\nFleet summary:\n")
print(summary_table, row.names = FALSE)
