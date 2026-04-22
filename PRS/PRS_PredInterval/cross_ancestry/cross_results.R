#!/usr/bin/env Rscript
# cross_results.R
# Aggregate cross-ancestry PredInterval results
# Compare BC vs Original scale across ancestries

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) > 0) args[1] else Sys.getenv("OUTBASE", unset = "")
if (base_dir == "") {
  stop("Usage: Rscript cross_results.R <base_dir>, or set OUTBASE env var")
}

output_file <- file.path(base_dir, "cross_ancestry_results_summary.txt")
sink(output_file, split = TRUE)

cat("Cross-Ancestry PredInterval Results\n")
cat("Base directory:", base_dir, "\n")
cat("Date:", format(Sys.time()), "\n\n")


results_all <- data.table()
coverage_results <- data.table()
high_risk_thresholds <- c(0.95, 0.99, 0.995, 0.999)

ancestries <- c("asn", "afr", "white_euro")

for (ancestry in ancestries) {
    for (scale in c("bc", "original")) {
        scale_dir <- file.path(base_dir, ancestry, paste0(scale, "_scale"))

        if (!dir.exists(scale_dir)) {
            cat("Directory not found:", scale_dir, "\n")
            next
        }

        traits <- list.dirs(scale_dir, recursive = FALSE, full.names = FALSE)

        for (trait in traits) {
            trait_dir <- file.path(scale_dir, trait)
            full_file <- file.path(trait_dir, paste0(trait, "_full_results.txt"))

            if (!file.exists(full_file)) {
                cat("  Skipping", trait, "(", ancestry, ",", scale, ")\n")
                next
            }

            full_results <- fread(full_file)

            coverage <- mean(full_results$covered, na.rm = TRUE)
            mean_width <- mean(full_results$width, na.rm = TRUE)
            sd_pheno <- sd(full_results$y, na.rm = TRUE)
            n_test <- nrow(full_results)

            coverage_results <- rbind(coverage_results, data.table(
                ancestry = ancestry, trait = trait, scale = scale,
                n_test = n_test, coverage = coverage,
                mean_width = mean_width, sd_pheno = sd_pheno
            ))

            for (thresh in high_risk_thresholds) {
                thresh_val <- quantile(full_results$y, thresh, na.rm = TRUE)

                full_results[, true_high_risk := y >= thresh_val]
                full_results[, identified_high_risk := upper_bound >= thresh_val]

                n_true_hr <- sum(full_results$true_high_risk)

                if (n_true_hr > 0) {
                    success_rate <- mean(full_results[true_high_risk == TRUE]$identified_high_risk,
                                         na.rm = TRUE)
                } else {
                    success_rate <- NA
                }

                n_identified <- sum(full_results$identified_high_risk)
                if (n_identified > 0) {
                    ppv <- sum(full_results$true_high_risk & full_results$identified_high_risk) / n_identified
                } else {
                    ppv <- NA
                }

                thresh_pct <- (1 - thresh) * 100
                thresh_label <- sprintf("%.1f", thresh_pct)

                results_all <- rbind(results_all, data.table(
                    ancestry = ancestry, trait = trait, scale = scale,
                    n_test = n_test, threshold_pct = thresh_label,
                    threshold_value = thresh_val, n_true_high_risk = n_true_hr,
                    n_identified = n_identified, coverage = coverage,
                    mean_width = mean_width, success_rate = success_rate, ppv = ppv
                ))
            }
        }
    }
}

fwrite(results_all, file.path(base_dir, "cross_ancestry_all_thresholds.txt"), sep = "\t")
fwrite(coverage_results, file.path(base_dir, "cross_ancestry_coverage.txt"), sep = "\t")


cat("COVERAGE COMPARISON BY ANCESTRY\n")

for (anc in ancestries) {
    cat(sprintf("--- %s ---\n", toupper(anc)))

    bc_cov <- coverage_results[ancestry == anc & scale == "bc"]
    orig_cov <- coverage_results[ancestry == anc & scale == "original"]

    if (nrow(bc_cov) == 0 || nrow(orig_cov) == 0) {
        cat("  Insufficient data\n\n")
        next
    }

    combined <- merge(bc_cov, orig_cov, by = "trait", suffixes = c("_bc", "_orig"))

    cat(sprintf("  Mean coverage - BC: %.4f, Original: %.4f (target = 0.95)\n",
                mean(combined$coverage_bc, na.rm = TRUE),
                mean(combined$coverage_orig, na.rm = TRUE)))
    cat(sprintf("  Mean width   - BC: %.4f, Original: %.4f\n\n",
                mean(combined$mean_width_bc, na.rm = TRUE),
                mean(combined$mean_width_orig, na.rm = TRUE)))
}


cat("HIGH-RISK IDENTIFICATION BY ANCESTRY\n")

for (thresh_pct in c("5.0", "1.0", "0.5", "0.1")) {
    cat(sprintf("=== TOP %s%% HIGH-RISK ===\n\n", thresh_pct))

    for (anc in ancestries) {
        cat(sprintf("  --- %s ---\n", toupper(anc)))

        subset <- results_all[ancestry == anc & threshold_pct == thresh_pct]
        bc <- subset[scale == "bc"]
        orig <- subset[scale == "original"]

        if (nrow(bc) == 0 || nrow(orig) == 0) {
            cat("    Insufficient data\n")
            next
        }

        combined <- merge(bc, orig, by = "trait", suffixes = c("_bc", "_orig"))

        # Compute success counts: success = success_rate * n_true_high_risk
        combined[, n_success_bc := round(success_rate_bc * n_true_high_risk_bc)]
        combined[, n_success_orig := round(success_rate_orig * n_true_high_risk_orig)]

        cat(sprintf("    Avg N true high-risk:  %.0f (BC) / %.0f (Orig)\n",
                    mean(combined$n_true_high_risk_bc, na.rm = TRUE),
                    mean(combined$n_true_high_risk_orig, na.rm = TRUE)))
        cat(sprintf("    BC success rate:       %.1f%% (SD: %.1f%%)  [avg %d/%d identified]\n",
                    mean(combined$success_rate_bc, na.rm = TRUE) * 100,
                    sd(combined$success_rate_bc, na.rm = TRUE) * 100,
                    round(mean(combined$n_success_bc, na.rm = TRUE)),
                    round(mean(combined$n_true_high_risk_bc, na.rm = TRUE))))
        cat(sprintf("    Original success rate: %.1f%% (SD: %.1f%%)  [avg %d/%d identified]\n",
                    mean(combined$success_rate_orig, na.rm = TRUE) * 100,
                    sd(combined$success_rate_orig, na.rm = TRUE) * 100,
                    round(mean(combined$n_success_orig, na.rm = TRUE)),
                    round(mean(combined$n_true_high_risk_orig, na.rm = TRUE))))

        if (nrow(combined) >= 3) {
            valid <- complete.cases(combined$success_rate_bc, combined$success_rate_orig)
            if (sum(valid) >= 3) {
                t_test <- t.test(combined$success_rate_bc[valid],
                                 combined$success_rate_orig[valid],
                                 paired = TRUE)
                diff <- combined$success_rate_bc[valid] - combined$success_rate_orig[valid]
                cohens_d <- mean(diff) / sd(diff)
                cat(sprintf("    Paired t-test p: %.2e, Cohen's d: %.2f\n",
                            t_test$p.value, cohens_d))
            }
        }
        cat("\n")
    }

    cat("  --- ALL ANCESTRIES COMBINED ---\n")
    subset_all <- results_all[threshold_pct == thresh_pct]
    bc_all <- subset_all[scale == "bc"]
    orig_all <- subset_all[scale == "original"]

    if (nrow(bc_all) > 0 && nrow(orig_all) > 0) {
        combined_all <- merge(bc_all, orig_all,
                              by = c("ancestry", "trait"), suffixes = c("_bc", "_orig"))
        cat(sprintf("    BC success rate:       %.1f%%\n",
                    mean(combined_all$success_rate_bc, na.rm = TRUE) * 100))
        cat(sprintf("    Original success rate: %.1f%%\n",
                    mean(combined_all$success_rate_orig, na.rm = TRUE) * 100))

        if (nrow(combined_all) >= 3) {
            valid <- complete.cases(combined_all$success_rate_bc, combined_all$success_rate_orig)
            if (sum(valid) >= 3) {
                t_test <- t.test(combined_all$success_rate_bc[valid],
                                 combined_all$success_rate_orig[valid],
                                 paired = TRUE)
                cat(sprintf("    Combined paired t-test p: %.2e\n", t_test$p.value))
            }
        }
    }
    cat("\n")
}


cat("INTERVAL WIDTH (RAW) & PGS R² BY ANCESTRY\n")
cat("(Width comparisons only meaningful within same scale)\n\n")

for (scale in c("original", "bc")) {
    cat(sprintf("--- %s SCALE ---\n\n", toupper(scale)))

    scale_data <- coverage_results[scale == (scale)]

    if (nrow(scale_data) == 0) {
        cat("  No data\n\n")
        next
    }

    cat(sprintf("  %-12s %10s %12s %12s %10s\n",
                "Ancestry", "N traits", "Mean width", "Mean width/SD", "Mean R²"))
    cat(paste0("  ", paste(rep("-", 58), collapse = "")), "\n")

    for (anc in ancestries) {
        anc_data <- scale_data[ancestry == anc]
        if (nrow(anc_data) == 0) next

        # Get cor(y, pgs) from summary files to compute R²
        r2_vals <- c()
        trait_dirs <- list.dirs(file.path(base_dir, anc, paste0(scale, "_scale")),
                                recursive = FALSE, full.names = TRUE)
        for (td in trait_dirs) {
            trait_name <- basename(td)
            summ_file <- file.path(td, paste0(trait_name, "_summary.txt"))
            if (file.exists(summ_file)) {
                summ <- fread(summ_file)
                if ("cor_y_pgs" %in% colnames(summ)) {
                    r2_vals <- c(r2_vals, summ$cor_y_pgs^2)
                }
            }
        }

        width_sd <- anc_data$mean_width / anc_data$sd_pheno

        cat(sprintf("  %-12s %10d %12.2f %12.2f %10.4f\n",
                    toupper(anc),
                    nrow(anc_data),
                    mean(anc_data$mean_width, na.rm = TRUE),
                    mean(width_sd, na.rm = TRUE),
                    if (length(r2_vals) > 0) mean(r2_vals, na.rm = TRUE) else NA))
    }
    cat("\n")
}

cat("Note: Width/SD is only comparable across ancestries WITHIN the same scale.\n")
cat("BC and original scales have different units, so cross-scale width/SD comparison is not meaningful.\n")


cat("PER-TRAIT DETAIL TABLES BY THRESHOLD\n")

for (thresh_pct in c("5.0", "1.0", "0.5", "0.1")) {
    cat(sprintf("\n--- TOP %s%% ---\n\n", thresh_pct))

    detail <- results_all[threshold_pct == thresh_pct]
    bc_detail <- detail[scale == "bc"]
    orig_detail <- detail[scale == "original"]
    combined_detail <- merge(bc_detail, orig_detail,
                              by = c("ancestry", "trait"), suffixes = c("_bc", "_orig"))

    if (nrow(combined_detail) == 0) {
        cat("  No data\n")
        next
    }

    cat(sprintf("%-12s %-22s %8s %8s %8s %14s %14s\n",
                "Ancestry", "Trait", "BC(%)", "Orig(%)", "Diff(%)", "BC(n/total)", "Orig(n/total)"))
    cat(paste(rep("-", 88), collapse = ""), "\n")

    for (i in 1:nrow(combined_detail)) {
        diff <- (combined_detail$success_rate_bc[i] - combined_detail$success_rate_orig[i]) * 100
        n_success_bc <- round(combined_detail$success_rate_bc[i] * combined_detail$n_true_high_risk_bc[i])
        n_success_orig <- round(combined_detail$success_rate_orig[i] * combined_detail$n_true_high_risk_orig[i])
        cat(sprintf("%-12s %-22s %7.1f%% %7.1f%% %+7.1f%% %6d/%-6d %6d/%-6d\n",
                    combined_detail$ancestry[i],
                    combined_detail$trait[i],
                    combined_detail$success_rate_bc[i] * 100,
                    combined_detail$success_rate_orig[i] * 100,
                    diff,
                    n_success_bc, combined_detail$n_true_high_risk_bc[i],
                    n_success_orig, combined_detail$n_true_high_risk_orig[i]))
    }
}

cat("Analysis complete!\n")

sink()
cat("\nSummary saved to:", output_file, "\n")
