#!/usr/bin/env Rscript
# multiconf_results.R
#
# Produces a 2D grid:
#   confidence_level x tail_threshold x tail_direction x scale
#
# Tail thresholds:
#   upper = top 50%, 25%, 10%, 5%, 1%, 0.1%
#   lower = bottom 50%, 25%, 10%, 5%, 1%, 0.1%
#
# Confidence levels:
#   75%, 80%, 85%, 90%, 95%

suppressPackageStartupMessages(library(data.table))

OUTBASE <- Sys.getenv("OUTBASE", unset = "")
if (OUTBASE == "") {
  stop("Environment variable OUTBASE must be set (e.g., source config.sh first)")
}

conf_levels <- c(0.75, 0.80, 0.85, 0.90, 0.95)
tail_pcts   <- c(50, 25, 10, 5, 1, 0.1)

# LOAD ALL RESULTS
all_results <- list()

for (scale in c("original", "bc")) {
    scale_dir <- file.path(OUTBASE, paste0(scale, "_scale"))
    if (!dir.exists(scale_dir)) next

    traits <- basename(list.dirs(scale_dir, recursive = FALSE))

    for (trait in traits) {
        multiconf_dir <- file.path(scale_dir, trait, "multiconf")
        if (!dir.exists(multiconf_dir)) next

        for (conf in conf_levels) {
            f <- file.path(multiconf_dir,
                sprintf("%s_conf%.2f_full_results.txt", trait, conf))
            if (!file.exists(f)) next

            dt <- fread(f)
            need <- c("id", "y", "pgs", "lower_bound", "upper_bound", "width", "covered")
            if (!all(need %in% names(dt))) next

            # Coverage
            all_results[[length(all_results) + 1]] <- data.table(
                trait = trait,
                scale = scale,
                conf_level = conf,
                tail_direction = NA_character_,
                tail_pct = NA_real_,
                metric = "coverage",
                value = mean(dt$covered, na.rm = TRUE)
            )

            # Tail metrics at each threshold (upper + lower)
            for (tp in tail_pcts) {
                p <- tp / 100

                # upper tail
                cutoff_up <- as.numeric(quantile(dt$y, probs = 1 - p, na.rm = TRUE))
                true_up   <- dt$y >= cutoff_up
                flag_up   <- dt$upper_bound >= cutoff_up

                n_true_up <- sum(true_up, na.rm = TRUE)
                n_both_up <- sum(true_up & flag_up, na.rm = TRUE)
                sr_up     <- if (n_true_up > 0) n_both_up / n_true_up else NA_real_

                all_results[[length(all_results) + 1]] <- data.table(
                    trait = trait,
                    scale = scale,
                    conf_level = conf,
                    tail_direction = "upper",
                    tail_pct = tp,
                    metric = "success_rate",
                    value = sr_up
                )

                # lower tail
                cutoff_low <- as.numeric(quantile(dt$y, probs = p, na.rm = TRUE))
                true_low   <- dt$y <= cutoff_low
                flag_low   <- dt$lower_bound <= cutoff_low

                n_true_low <- sum(true_low, na.rm = TRUE)
                n_both_low <- sum(true_low & flag_low, na.rm = TRUE)
                sr_low     <- if (n_true_low > 0) n_both_low / n_true_low else NA_real_

                all_results[[length(all_results) + 1]] <- data.table(
                    trait = trait,
                    scale = scale,
                    conf_level = conf,
                    tail_direction = "lower",
                    tail_pct = tp,
                    metric = "success_rate",
                    value = sr_low
                )
            }
        }
    }
}

results <- rbindlist(all_results, fill = TRUE)

if (nrow(results) == 0) {
    cat("No results found. Run multi_conf.sh first.\n")
    quit(status = 1)
}

fwrite(results, file.path(OUTBASE, "multiconf_all_raw.txt"), sep = "\t")

# COMPUTE HEATMAP DATA: BC - Original difference
sr_data <- results[metric == "success_rate"]

wide <- dcast(
    sr_data,
    trait + tail_direction + conf_level + tail_pct ~ scale,
    value.var = "value"
)

if (!all(c("bc", "original") %in% names(wide))) {
    cat("Need both scales to compute differences.\n")
    quit(status = 1)
}

wide[, diff := bc - original]

# Average across traits
heatmap <- wide[, .(
    mean_sr_bc = mean(bc, na.rm = TRUE),
    mean_sr_orig = mean(original, na.rm = TRUE),
    mean_diff = mean(diff, na.rm = TRUE),
    sd_diff = sd(diff, na.rm = TRUE),
    n_traits = sum(!is.na(diff)),
    p_value = tryCatch(
        t.test(diff[!is.na(diff)])$p.value,
        error = function(e) NA_real_
    )
), by = .(tail_direction, conf_level, tail_pct)]

fwrite(heatmap, file.path(OUTBASE, "multiconf_heatmap.txt"), sep = "\t")
fwrite(wide, file.path(OUTBASE, "multiconf_wide.txt"), sep = "\t")

# COVERAGE CHECK
cov_data <- results[metric == "coverage"]
cov_wide <- dcast(
    cov_data,
    trait + conf_level ~ scale,
    value.var = "value"
)

if (all(c("bc", "original") %in% names(cov_wide))) {
    cov_summary <- cov_wide[, .(
        mean_cov_bc = mean(bc, na.rm = TRUE),
        mean_cov_orig = mean(original, na.rm = TRUE)
    ), by = conf_level]
    fwrite(cov_summary, file.path(OUTBASE, "multiconf_coverage_check.txt"), sep = "\t")
}

# HUMAN-READABLE SUMMARY
summary_file <- file.path(OUTBASE, "multiconf_summary.txt")
sink(summary_file, split = TRUE)

cat("MULTI-CONFIDENCE LEVEL x TAIL THRESHOLD ANALYSIS\n")
cat("Success rate difference (BC - Original), averaged across traits\n")

cat("--- COVERAGE CHECK (mean across traits) ---\n\n")
cat(sprintf("%10s %12s %12s\n", "Conf", "BC", "Original"))
cat(paste(rep("-", 36), collapse = ""), "\n")
if (exists("cov_summary")) {
    for (i in 1:nrow(cov_summary)) {
        cat(sprintf("%10.0f%% %11.4f %11.4f\n",
            cov_summary$conf_level[i] * 100,
            cov_summary$mean_cov_bc[i],
            cov_summary$mean_cov_orig[i]))
    }
}

# helper printer for each tail direction
print_tail_block <- function(dir_name) {
    sub_heat <- heatmap[tail_direction == dir_name]
    if (nrow(sub_heat) == 0) return(invisible(NULL))

    cat(sprintf("\n\n--- %s TAIL: SUCCESS RATE DIFFERENCE (BC - Orig, percentage points) ---\n\n",
        toupper(dir_name)))
    cat(sprintf("%10s", "Conf\\Tail"))
    for (tp in tail_pcts) cat(sprintf(" %8s", paste0(tp, "%")))
    cat("\n")
    cat(paste(rep("-", 10 + 9 * length(tail_pcts)), collapse = ""), "\n")

    for (conf in conf_levels) {
        cat(sprintf("%10.0f%%", conf * 100))
        for (tp in tail_pcts) {
            row <- sub_heat[conf_level == conf & tail_pct == tp]
            if (nrow(row) == 1 && !is.na(row$mean_diff)) {
                cat(sprintf(" %+7.1f%%", row$mean_diff * 100))
            } else {
                cat(sprintf(" %8s", "N/A"))
            }
        }
        cat("\n")
    }

    cat(sprintf("\n--- %s TAIL: SIGNIFICANCE (paired t-test across traits) ---\n\n",
        toupper(dir_name)))
    cat(sprintf("%10s", "Conf\\Tail"))
    for (tp in tail_pcts) cat(sprintf(" %8s", paste0(tp, "%")))
    cat("\n")
    cat(paste(rep("-", 10 + 9 * length(tail_pcts)), collapse = ""), "\n")

    for (conf in conf_levels) {
        cat(sprintf("%10.0f%%", conf * 100))
        for (tp in tail_pcts) {
            row <- sub_heat[conf_level == conf & tail_pct == tp]
            if (nrow(row) == 1 && !is.na(row$p_value)) {
                if (row$p_value < 0.001) {
                    star <- "***"
                } else if (row$p_value < 0.01) {
                    star <- "**"
                } else if (row$p_value < 0.05) {
                    star <- "*"
                } else {
                    star <- "n.s."
                }
                cat(sprintf(" %8s", star))
            } else {
                cat(sprintf(" %8s", "N/A"))
            }
        }
        cat("\n")
    }

    cat(sprintf("\n--- %s TAIL: BC SUCCESS RATE (%%) ---\n\n",
        toupper(dir_name)))
    cat(sprintf("%10s", "Conf\\Tail"))
    for (tp in tail_pcts) cat(sprintf(" %8s", paste0(tp, "%")))
    cat("\n")
    cat(paste(rep("-", 10 + 9 * length(tail_pcts)), collapse = ""), "\n")

    for (conf in conf_levels) {
        cat(sprintf("%10.0f%%", conf * 100))
        for (tp in tail_pcts) {
            row <- sub_heat[conf_level == conf & tail_pct == tp]
            if (nrow(row) == 1 && !is.na(row$mean_sr_bc)) {
                cat(sprintf(" %7.1f%%", row$mean_sr_bc * 100))
            } else {
                cat(sprintf(" %8s", "N/A"))
            }
        }
        cat("\n")
    }

    cat(sprintf("\n--- %s TAIL: ORIGINAL SUCCESS RATE (%%) ---\n\n",
        toupper(dir_name)))
    cat(sprintf("%10s", "Conf\\Tail"))
    for (tp in tail_pcts) cat(sprintf(" %8s", paste0(tp, "%")))
    cat("\n")
    cat(paste(rep("-", 10 + 9 * length(tail_pcts)), collapse = ""), "\n")

    for (conf in conf_levels) {
        cat(sprintf("%10.0f%%", conf * 100))
        for (tp in tail_pcts) {
            row <- sub_heat[conf_level == conf & tail_pct == tp]
            if (nrow(row) == 1 && !is.na(row$mean_sr_orig)) {
                cat(sprintf(" %7.1f%%", row$mean_sr_orig * 100))
            } else {
                cat(sprintf(" %8s", "N/A"))
            }
        }
        cat("\n")
    }
}

print_tail_block("upper")
print_tail_block("lower")

cat("Analysis complete.\n")

sink()

cat("Wrote:\n")
cat(" -", file.path(OUTBASE, "multiconf_all_raw.txt"), "\n")
cat(" -", file.path(OUTBASE, "multiconf_heatmap.txt"), "\n")
cat(" -", file.path(OUTBASE, "multiconf_wide.txt"), "\n")
cat(" -", file.path(OUTBASE, "multiconf_coverage_check.txt"), "\n")
cat(" -", summary_file, "\n")
