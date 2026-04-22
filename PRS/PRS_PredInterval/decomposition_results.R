#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# CONFIG
OUTBASE <- Sys.getenv("OUTBASE", unset = "")
if (OUTBASE == "") {
  stop("Environment variable OUTBASE must be set (e.g., source config.sh first)")
}

dirs <- list(
    A = file.path(OUTBASE, "original_scale"),      # original PGS + original resid
    B = file.path(OUTBASE, "hybrid_bc2orig"),      # BC PGS (back-transformed) + original resid
    C = file.path(OUTBASE, "bc_scale"),            # BC PGS + BC resid
    D = file.path(OUTBASE, "hybrid_orig2bc")       # orig PGS (forward-transformed) + BC resid
)

cell_labels <- c(
    A = "Orig PGS + Orig Resid",
    B = "BC PGS + Orig Resid",
    C = "BC PGS + BC Resid",
    D = "Orig PGS + BC Resid"
)

tail_pcts <- c(10, 5, 2, 1, 0.5, 0.1)

# HELPERS
safe_read_results <- function(scale_dir, pheno) {
    f <- file.path(scale_dir, pheno, paste0(pheno, "_full_results.txt"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f)
    need <- c("id", "y", "pgs", "lower_bound", "upper_bound", "width", "covered")
    miss <- setdiff(need, names(dt))
    if (length(miss) > 0) {
        warning("Missing columns in ", f, ": ", paste(miss, collapse = ", "))
        return(NULL)
    }
    dt[]
}

compute_tail_metrics <- function(dt, tail_pct, direction) {
    p <- tail_pct / 100
    if (direction == "upper") {
        cutoff    <- as.numeric(quantile(dt$y, probs = 1 - p, na.rm = TRUE, type = 7))
        true_tail <- dt$y >= cutoff
        flagged   <- dt$upper_bound >= cutoff
    } else {
        cutoff    <- as.numeric(quantile(dt$y, probs = p, na.rm = TRUE, type = 7))
        true_tail <- dt$y <= cutoff
        flagged   <- dt$lower_bound <= cutoff
    }
    n_true <- sum(true_tail, na.rm = TRUE)
    n_flag <- sum(flagged, na.rm = TRUE)
    n_both <- sum(true_tail & flagged, na.rm = TRUE)
    data.table(
        tail_direction = direction, threshold_pct = tail_pct, cutoff = cutoff,
        n_true_tail = n_true, n_flagged = n_flag, n_true_and_flagged = n_both,
        success_rate = if (n_true > 0) n_both / n_true else NA_real_,
        ppv = if (n_flag > 0) n_both / n_flag else NA_real_
    )
}

summarize_scale <- function(dt, pheno, cell_name) {
    data.table(
        trait = pheno, cell = cell_name,
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

# DISCOVER TRAITS
all_traits <- character(0)
for (d in dirs) {
    if (dir.exists(d)) {
        all_traits <- c(all_traits, basename(list.dirs(d, recursive = FALSE)))
    }
}
traits <- sort(unique(all_traits))

# MAIN LOOP
coverage_results <- list()
tail_results     <- list()

for (pheno in traits) {
    cells_data <- list()
    for (cell in names(dirs)) {
        cells_data[[cell]] <- safe_read_results(dirs[[cell]], pheno)
    }

    for (cell in names(cells_data)) {
        dt <- cells_data[[cell]]
        if (is.null(dt)) next

        coverage_results[[paste0(pheno, "_", cell)]] <-
            summarize_scale(dt, pheno, cell)

        for (direction in c("upper", "lower")) {
            tmp <- rbindlist(lapply(tail_pcts, function(tp)
                compute_tail_metrics(dt, tp, direction)))
            tmp[, `:=`(trait = pheno, cell = cell)]
            tail_results[[paste0(pheno, "_", cell, "_", direction)]] <- tmp
        }
    }
}

coverage_dt <- if (length(coverage_results)) rbindlist(coverage_results, fill = TRUE) else data.table()
tail_dt     <- if (length(tail_results))     rbindlist(tail_results, fill = TRUE)     else data.table()

if (nrow(tail_dt)) {
    setcolorder(tail_dt, c(
        "trait", "cell", "tail_direction", "threshold_pct", "cutoff",
        "n_true_tail", "n_flagged", "n_true_and_flagged", "success_rate", "ppv"
    ))
}

# SAVE RAW OUTPUTS
fwrite(coverage_dt, file.path(OUTBASE, "decomp_coverage.txt"),  sep = "\t")
fwrite(tail_dt,     file.path(OUTBASE, "decomp_tail_all.txt"),  sep = "\t")

# WIDE FORMAT: 4-cell comparison
if (nrow(coverage_dt)) {
    cov_wide <- dcast(coverage_dt, trait ~ cell,
        value.var = c("observed_coverage", "mean_width", "mean_width_over_sd_y"))
    fwrite(cov_wide, file.path(OUTBASE, "decomp_coverage_wide.txt"), sep = "\t")
}

if (nrow(tail_dt)) {
    tail_wide <- dcast(tail_dt,
        trait + tail_direction + threshold_pct ~ cell,
        value.var = c("success_rate", "ppv"))

    sr <- function(x) paste0("success_rate_", x)
    have_all <- all(c(sr("A"), sr("B"), sr("C"), sr("D")) %in% names(tail_wide))

    if (have_all) {
        tail_wide[, `:=`(
            # Explicit five pair comparisons
            pair_CA = get(sr("C")) - get(sr("A")),  # total
            pair_DA = get(sr("D")) - get(sr("A")),  # residual effect holding orig PGS fixed
            pair_CB = get(sr("C")) - get(sr("B")),  # residual effect holding BC PGS fixed
            pair_BA = get(sr("B")) - get(sr("A")),  # PGS effect holding orig resid fixed
            pair_CD = get(sr("C")) - get(sr("D")),  # PGS effect holding BC resid fixed

            # Existing decomposition summaries
            row_effect_orig_pgs  = get(sr("D")) - get(sr("A")),  # D - A
            row_effect_bc_pgs    = get(sr("C")) - get(sr("B")),  # C - B
            col_effect_orig_res  = get(sr("B")) - get(sr("A")),  # B - A
            col_effect_bc_res    = get(sr("C")) - get(sr("D")),  # C - D
            gain_total           = get(sr("C")) - get(sr("A")),

            avg_row_effect = (get(sr("D")) - get(sr("A")) + get(sr("C")) - get(sr("B"))) / 2,
            avg_col_effect = (get(sr("B")) - get(sr("A")) + get(sr("C")) - get(sr("D"))) / 2,

            pct_from_residual = fifelse(
                get(sr("C")) - get(sr("A")) != 0,
                ((get(sr("D")) - get(sr("A")) + get(sr("C")) - get(sr("B"))) / 2) /
                    (get(sr("C")) - get(sr("A"))) * 100,
                NA_real_),
            pct_from_pgs = fifelse(
                get(sr("C")) - get(sr("A")) != 0,
                ((get(sr("B")) - get(sr("A")) + get(sr("C")) - get(sr("D"))) / 2) /
                    (get(sr("C")) - get(sr("A"))) * 100,
                NA_real_)
        )]
    }

    fwrite(tail_wide, file.path(OUTBASE, "decomp_tail_wide.txt"), sep = "\t")
}

# HUMAN-READABLE SUMMARY
summary_file <- file.path(OUTBASE, "decomp_summary.txt")
sink(summary_file, split = TRUE)

cat("FULL 2x2 DECOMPOSITION\n")
cat("                  Original PGS              BC PGS\n")
cat("Orig residuals    A (baseline)              B (BC PGS, orig resid)\n")
cat("BC residuals      D (orig PGS, BC resid)    C (full SIQReg)\n\n")
cat("Row effect  = switching residual scale (D-A, C-B)\n")
cat("Col effect  = switching PGS (B-A, C-D)\n")
cat("If residual calibration dominates: C ≈ D >> A ≈ B\n")
cat("If PGS quality dominates:          C ≈ B >> A ≈ D\n\n")

cat("--- COVERAGE ---\n")
if (nrow(coverage_dt)) {
    cov_print <- dcast(coverage_dt, trait ~ cell, value.var = "observed_coverage")
    print(cov_print, digits = 4)
}
cat("\n")

if (exists("tail_wide") && nrow(tail_wide)) {
    for (dir in c("upper", "lower")) {
        cat(sprintf("\n--- %s TAIL SUCCESS RATES ---\n\n", toupper(dir)))

        sub <- tail_wide[tail_direction == dir]
        if (nrow(sub) == 0) next

        for (tp in tail_pcts) {
            sub_tp <- sub[threshold_pct == tp]
            if (nrow(sub_tp) == 0) next

            cat(sprintf("=== %s tail, top %.1f%% ===\n", dir, tp))
            cat(sprintf("%-20s %7s %7s %7s %7s | %8s %8s | %6s %6s\n",
                "Trait", "A", "B", "D", "C", "RowEff", "ColEff", "%Res", "%PGS"))
            cat(paste(rep("-", 100), collapse = ""), "\n")

            for (i in 1:nrow(sub_tp)) {
                r <- sub_tp[i]
                fmt_sr <- function(x) if (is.na(x)) "   N/A" else sprintf("%6.1f%%", x * 100)
                fmt_ef <- function(x) if (is.na(x)) "    N/A" else sprintf("%+7.1f%%", x * 100)
                fmt_pc <- function(x) if (is.na(x)) "  N/A" else sprintf("%5.0f%%", x)

                cat(sprintf("%-20s %7s %7s %7s %7s | %8s %8s | %6s %6s\n",
                    r$trait,
                    fmt_sr(r[[sr("A")]]), fmt_sr(r[[sr("B")]]),
                    fmt_sr(r[[sr("D")]]), fmt_sr(r[[sr("C")]]),
                    fmt_ef(r$avg_row_effect), fmt_ef(r$avg_col_effect),
                    fmt_pc(r$pct_from_residual), fmt_pc(r$pct_from_pgs)))
            }

            cat(paste(rep("-", 100), collapse = ""), "\n")
            cat(sprintf("%-20s %6.1f%% %6.1f%% %6.1f%% %6.1f%% | %+7.1f%% %+7.1f%% | %5.0f%% %5.0f%%\n",
                "MEAN",
                mean(sub_tp[[sr("A")]], na.rm = TRUE) * 100,
                mean(sub_tp[[sr("B")]], na.rm = TRUE) * 100,
                mean(sub_tp[[sr("D")]], na.rm = TRUE) * 100,
                mean(sub_tp[[sr("C")]], na.rm = TRUE) * 100,
                mean(sub_tp$avg_row_effect, na.rm = TRUE) * 100,
                mean(sub_tp$avg_col_effect, na.rm = TRUE) * 100,
                mean(sub_tp$pct_from_residual, na.rm = TRUE),
                mean(sub_tp$pct_from_pgs, na.rm = TRUE)))
            cat("\n\n")
        }
    }
}

cat("\n--- FIVE PAIR COMPARISONS ---\n\n")

if (exists("tail_wide") && "pair_CA" %in% names(tail_wide)) {
    pair_info <- list(
        pair_CA = "C - A : total SIQReg effect",
        pair_DA = "D - A : residual effect, holding original PGS fixed",
        pair_CB = "C - B : residual effect, holding BC PGS fixed",
        pair_BA = "B - A : PGS effect, holding original residual fixed",
        pair_CD = "C - D : PGS effect, holding BC residual fixed"
    )

    for (dir in c("upper", "lower")) {
        cat(sprintf("%s tail\n", toupper(dir)))
        for (tp in tail_pcts) {
            sub <- tail_wide[tail_direction == dir & threshold_pct == tp]
            if (nrow(sub) == 0) next

            cat(sprintf("  threshold %.1f%% (n = %d traits)\n", tp, nrow(sub)))
            for (nm in names(pair_info)) {
                x <- sub[[nm]]
                x <- x[is.finite(x)]
                if (length(x) < 3) next
                tt <- t.test(x)
                cat(sprintf("    %-45s mean = %+6.2f%%, p = %.3e\n",
                    pair_info[[nm]], mean(x) * 100, tt$p.value))
            }
            cat("\n")
        }
        cat("\n")
    }
}

cat("\n--- STATISTICAL TESTS ---\n\n")

if (exists("tail_wide") && "avg_row_effect" %in% names(tail_wide)) {
    for (dir in c("upper", "lower")) {
        for (tp in tail_pcts) {
            sub <- tail_wide[tail_direction == dir & threshold_pct == tp]
            valid <- complete.cases(sub$avg_row_effect, sub$avg_col_effect, sub$gain_total)
            if (sum(valid) < 3) next

            cat(sprintf("%s tail, top %.1f%% (n = %d traits):\n", dir, tp, sum(valid)))

            t_row <- t.test(sub$avg_row_effect[valid])
            t_col <- t.test(sub$avg_col_effect[valid])
            t_tot <- t.test(sub$gain_total[valid])

            cat(sprintf("  Row effect (residual): mean = %+.2f%%, p = %.3e\n",
                mean(sub$avg_row_effect[valid]) * 100, t_row$p.value))
            cat(sprintf("  Col effect (PGS):      mean = %+.2f%%, p = %.3e\n",
                mean(sub$avg_col_effect[valid]) * 100, t_col$p.value))
            cat(sprintf("  Total (C-A):           mean = %+.2f%%, p = %.3e\n",
                mean(sub$gain_total[valid]) * 100, t_tot$p.value))

            diff <- sub$avg_row_effect[valid] - sub$avg_col_effect[valid]
            if (sd(diff) > 0) {
                t_diff <- t.test(diff)
                cat(sprintf("  Row > Col?             mean diff = %+.2f%%, p = %.3e\n",
                    mean(diff) * 100, t_diff$p.value))
            }
            cat("\n")
        }
    }
}

cat("Analysis complete.\n")

sink()

# PRINT FILE LIST
cat("Wrote:\n")
cat(" -", file.path(OUTBASE, "decomp_coverage.txt"), "\n")
cat(" -", file.path(OUTBASE, "decomp_tail_all.txt"), "\n")
cat(" -", file.path(OUTBASE, "decomp_coverage_wide.txt"), "\n")
cat(" -", file.path(OUTBASE, "decomp_tail_wide.txt"), "\n")
cat(" -", summary_file, "\n")
