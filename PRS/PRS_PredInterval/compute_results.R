#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# CONFIG
OUTBASE <- Sys.getenv("OUTBASE", unset = "")
if (OUTBASE == "") {
  stop("Environment variable OUTBASE must be set (e.g., source config.sh first)")
}

orig_dir <- file.path(OUTBASE, "original_scale")
bc_dir   <- file.path(OUTBASE, "bc_scale")

# Add coarser + finer tails here
tail_pcts <- c(10, 5, 2, 1, 0.5, 0.1)

# HELPERS
safe_read_results <- function(scale_dir, pheno) {
  f <- file.path(scale_dir, pheno, paste0(pheno, "_full_results.txt"))
  if (!file.exists(f)) return(NULL)
  dt <- fread(f)
  need <- c("id", "y", "pgs", "lower_bound", "upper_bound", "width", "covered")
  miss <- setdiff(need, names(dt))
  if (length(miss) > 0) {
    stop("Missing columns in ", f, ": ", paste(miss, collapse = ", "))
  }
  dt[]
}

compute_tail_metrics <- function(dt, tail_pct, direction) {
  p <- tail_pct / 100

  if (direction == "upper") {
    cutoff <- as.numeric(quantile(dt$y, probs = 1 - p, na.rm = TRUE, type = 7))
    true_tail <- dt$y >= cutoff
    flagged   <- dt$upper_bound >= cutoff
  } else if (direction == "lower") {
    cutoff <- as.numeric(quantile(dt$y, probs = p, na.rm = TRUE, type = 7))
    true_tail <- dt$y <= cutoff
    flagged   <- dt$lower_bound <= cutoff
  } else {
    stop("direction must be 'upper' or 'lower'")
  }

  n_true   <- sum(true_tail, na.rm = TRUE)
  n_flag   <- sum(flagged, na.rm = TRUE)
  n_both   <- sum(true_tail & flagged, na.rm = TRUE)

  success_rate <- if (n_true > 0) n_both / n_true else NA_real_
  ppv          <- if (n_flag > 0) n_both / n_flag else NA_real_

  data.table(
    tail_direction = direction,
    threshold_pct = tail_pct,
    cutoff = cutoff,
    n_true_tail = n_true,
    n_flagged = n_flag,
    n_true_and_flagged = n_both,
    success_rate = success_rate,
    ppv = ppv
  )
}

summarize_scale <- function(dt, pheno, scale_name) {
  data.table(
    trait = pheno,
    scale = scale_name,
    n_test = nrow(dt),
    observed_coverage = mean(dt$covered, na.rm = TRUE),
    mean_width = mean(dt$width, na.rm = TRUE),
    sd_width = sd(dt$width, na.rm = TRUE),
    sd_y = sd(dt$y, na.rm = TRUE),
    mean_width_over_sd_y = mean(dt$width, na.rm = TRUE) / sd(dt$y, na.rm = TRUE),
    mean_pgs = mean(dt$pgs, na.rm = TRUE),
    sd_pgs = sd(dt$pgs, na.rm = TRUE)
  )
}

# FIND TRAITS
orig_traits <- basename(list.dirs(orig_dir, recursive = FALSE, full.names = TRUE))
bc_traits   <- basename(list.dirs(bc_dir, recursive = FALSE, full.names = TRUE))
traits <- sort(unique(c(orig_traits, bc_traits)))

# MAIN
coverage_results <- list()
tail_results <- list()

for (pheno in traits) {
  dt_orig <- safe_read_results(orig_dir, pheno)
  dt_bc   <- safe_read_results(bc_dir, pheno)

  if (!is.null(dt_orig)) {
    coverage_results[[paste0(pheno, "_orig")]] <- summarize_scale(dt_orig, pheno, "original")

    tmp <- rbindlist(lapply(tail_pcts, function(tp) compute_tail_metrics(dt_orig, tp, "upper")))
    tmp[, `:=`(trait = pheno, scale = "original")]
    tail_results[[paste0(pheno, "_orig_upper")]] <- tmp

    tmp <- rbindlist(lapply(tail_pcts, function(tp) compute_tail_metrics(dt_orig, tp, "lower")))
    tmp[, `:=`(trait = pheno, scale = "original")]
    tail_results[[paste0(pheno, "_orig_lower")]] <- tmp
  }

  if (!is.null(dt_bc)) {
    coverage_results[[paste0(pheno, "_bc")]] <- summarize_scale(dt_bc, pheno, "bc")

    tmp <- rbindlist(lapply(tail_pcts, function(tp) compute_tail_metrics(dt_bc, tp, "upper")))
    tmp[, `:=`(trait = pheno, scale = "bc")]
    tail_results[[paste0(pheno, "_bc_upper")]] <- tmp

    tmp <- rbindlist(lapply(tail_pcts, function(tp) compute_tail_metrics(dt_bc, tp, "lower")))
    tmp[, `:=`(trait = pheno, scale = "bc")]
    tail_results[[paste0(pheno, "_bc_lower")]] <- tmp
  }
}

coverage_dt <- if (length(coverage_results)) rbindlist(coverage_results, fill = TRUE) else data.table()
tail_dt     <- if (length(tail_results)) rbindlist(tail_results, fill = TRUE) else data.table()

# Order columns nicely
if (nrow(tail_dt)) {
  setcolorder(
    tail_dt,
    c("trait", "scale", "tail_direction", "threshold_pct", "cutoff",
      "n_true_tail", "n_flagged", "n_true_and_flagged", "success_rate", "ppv")
  )
}

# SAVE MAIN OUTPUTS
fwrite(coverage_dt, file.path(OUTBASE, "results_coverage.txt"), sep = "\t")
fwrite(tail_dt, file.path(OUTBASE, "results_all_thresholds.txt"), sep = "\t")

# PAIRED COMPARISONS: BC vs ORIGINAL
if (nrow(coverage_dt)) {
  cov_wide <- dcast(
    coverage_dt,
    trait ~ scale,
    value.var = c("observed_coverage", "mean_width", "mean_width_over_sd_y")
  )
  fwrite(cov_wide, file.path(OUTBASE, "results_coverage_wide.txt"), sep = "\t")
}

if (nrow(tail_dt)) {
  tail_wide <- dcast(
    tail_dt,
    trait + tail_direction + threshold_pct ~ scale,
    value.var = c("success_rate", "ppv", "cutoff", "n_true_tail", "n_flagged", "n_true_and_flagged")
  )
  fwrite(tail_wide, file.path(OUTBASE, "results_all_thresholds_wide.txt"), sep = "\t")
}

# HUMAN-READABLE SUMMARY
summary_file <- file.path(OUTBASE, "results_summary.txt")
con <- file(summary_file, open = "wt")

writeLines("PredInterval results summary", con)
writeLines("", con)

if (nrow(coverage_dt)) {
  writeLines("Coverage/width summary by trait and scale:", con)
  writeLines(capture.output(print(coverage_dt)), con)
  writeLines("", con)
}

if (nrow(tail_dt)) {
  writeLines("Tail-identification summary by trait, scale, tail direction, and threshold:", con)
  writeLines(capture.output(print(tail_dt)), con)
  writeLines("", con)
}

if (exists("cov_wide") && nrow(cov_wide)) {
  writeLines("BC vs original coverage/width comparison:", con)
  writeLines(capture.output(print(cov_wide)), con)
  writeLines("", con)
}

if (exists("tail_wide") && nrow(tail_wide)) {
  writeLines("BC vs original tail comparison:", con)
  writeLines(capture.output(print(tail_wide)), con)
  writeLines("", con)
}

close(con)

cat("Wrote:\n")
cat(" -", file.path(OUTBASE, "results_coverage.txt"), "\n")
cat(" -", file.path(OUTBASE, "results_all_thresholds.txt"), "\n")
cat(" -", file.path(OUTBASE, "results_coverage_wide.txt"), "\n")
cat(" -", file.path(OUTBASE, "results_all_thresholds_wide.txt"), "\n")
cat(" -", file.path(OUTBASE, "results_summary.txt"), "\n")
