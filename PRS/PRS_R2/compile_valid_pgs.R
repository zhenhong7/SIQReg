#!/usr/bin/env Rscript
## compile_valid_pgs.R (REVISED, minimal change + add spearman)

BOOTSTRAP_PRS_PATH <- Sys.getenv("BOOTSTRAP_PRS_PATH",
  unset = "/path/to/SI_GWAS/bootstrap_prs.R")
source(BOOTSTRAP_PRS_PATH)
suppressPackageStartupMessages(library(dplyr))

R2_TEST <- "PEARSON"
B_BOOT  <- 1000L

PHENO_DIGITS <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(PHENO_DIGITS) || PHENO_DIGITS == "") {
  stop("Usage: Rscript compile_valid_pgs.R <PHENO_DIGITS>")
}

PHENO_NO_DIGITS <- ifelse(grepl("EA4", PHENO_DIGITS), "EA4",
                   ifelse(grepl("FEV1674206", PHENO_DIGITS), "FEV1",
                   ifelse(grepl("IGF-1674178", PHENO_DIGITS), "IGF-1",
                          gsub("[0-9]+$", "", PHENO_DIGITS))))
PHENO_LOWER <- tolower(PHENO_NO_DIGITS)

PGS_COLS <- c("FID","IID","ALLELE_CT","ALLELE_DOSAGE_SUM","SCORE1_AVG","SCORE1_SUM")

PHENO_DIR <- Sys.getenv("PHENO_DIR", unset = "/path/to/extracted_phenotypes/")
PHENO_FILE <- paste0(PHENO_DIR, PHENO_NO_DIGITS, "/", PHENO_DIGITS, ".pheno")
pheno_table <- read.table(PHENO_FILE, header=TRUE)

if (!(R2_TEST %in% c("SPEARMAN", "PEARSON"))) {
  stop("R2_TEST must be either 'SPEARMAN' or 'PEARSON'")
}

META_TAB <- Sys.getenv("META_TAB",
  unset = "/path/to/SI_GWAS/meta_both_summary_all.txt")
meta <- read.table(META_TAB, header=TRUE, sep="\t", stringsAsFactors=FALSE)
ii <- match(PHENO_NO_DIGITS, meta$phenotype)
if (is.na(ii)) stop("Cannot find phenotype in meta table: ", PHENO_NO_DIGITS)
lambda <- as.numeric(meta$lambda_WM[ii])

inv_boxcox_prs <- function(z, lam) {
  z <- as.numeric(z)
  if (!is.finite(lam)) return(rep(NA_real_, length(z)))
  if (abs(lam) < 1e-12) return(exp(z))
  u <- 1 + lam * z
  out <- rep(NA_real_, length(z))
  ok <- is.finite(u) & (u > 0)
  out[ok] <- u[ok]^(1/lam)
  out
}

cat("pheno:", PHENO_LOWER, " lambda_WM:", lambda, "\n")

BC_R2_BASE   <- Sys.getenv("BC_R2_BASE",
  unset = "/path/to/SI_GWAS/results/pgs_outputs/")
ORIG_R2_BASE <- Sys.getenv("ORIG_R2_BASE",
  unset = "/path/to/SI_GWAS/original_scale/results/pgs_outputs/")

BC_PGS_BASE   <- Sys.getenv("BC_PGS_BASE",
  unset = "/path/to/SI_GWAS/results/pgs_outputs/pgs_out/")
ORIG_PGS_BASE <- Sys.getenv("ORIG_PGS_BASE",
  unset = "/path/to/SI_GWAS/original_scale/results/pgs_outputs/pgs_out/")

OUT_BASE <- Sys.getenv("OUT_BASE",
  unset = "/path/to/SI_GWAS/results/pgs_outputs/compare_bc_vs_orig/")

## MINIMAL CHANGE: decide by metric (spearman vs others), not by R2_TEST
pick_processed_file <- function(base_dir, pheno_lower, metric, lambda_val) {
  is_bc <- identical(base_dir, BC_R2_BASE)

  if (metric == "spearman") {
    if (is_bc) {
      cand <- c("r2_out/", "r2_out_log/")
    } else {
      cand <- c("r2_out_spearman/", "r2_out_spearman_log/")
    }
  } else {
    cand <- c("r2_out/", "r2_out_log/")
  }

  fn <- paste0("processed_", pheno_lower, "_whitebrit_PGS_", metric, ".txt")

  f1 <- file.path(base_dir, cand[1], fn)
  if (file.exists(f1)) return(f1)

  ## keep your "try log dir if lambda==0" behavior
  if (is.finite(lambda_val) && abs(lambda_val) < 1e-12) {
    f2 <- file.path(base_dir, cand[2], fn)
    if (file.exists(f2)) return(f2)
  }

  f2 <- file.path(base_dir, cand[2], fn)
  if (file.exists(f2)) return(f2)

  stop("Cannot find processed file: ", fn, " under ", base_dir)
}

COVAR_FILE <- Sys.getenv("COVAR_FILE",
  unset = paste0(PHENO_DIR, "covariates_sa40PC/covariates_sa40PC674178.pheno"))
covs_table <- read.table(COVAR_FILE, header=TRUE)[,1:14]

## MINIMAL CHANGE: add "spearman"
metrics <- c("std", "spearman", "rankLin", "rankEC")

for (metric in metrics) {

  ## best thresholds (bc + original) for this metric
  f_bc   <- pick_processed_file(BC_R2_BASE,   PHENO_LOWER, metric, lambda)
  f_orig <- pick_processed_file(ORIG_R2_BASE, PHENO_LOWER, metric, lambda)

  best0 <- read.table(f_bc,   header=TRUE)
  best1 <- read.table(f_orig, header=TRUE)

  ## same as your original: first row / first col is threshold
  thresh0 <- format(best0[1,1], scientific=FALSE)
  thresh1 <- format(best1[1,1], scientific=FALSE)

  cat("[", metric, "] thresh_bc=", thresh0, " thresh_orig=", thresh1, "\n", sep="")

  ## output folder per metric
  out_dir <- file.path(OUT_BASE, metric)
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

  for (group in c("all","female","male")) {
    for (valid_pop in c("white_euro","afr","asn")) {

      out_file <- file.path(out_dir,
                            paste0("bc_vs_orig_", PHENO_LOWER, "_", valid_pop, "_", group, ".txt"))
      if (file.exists(out_file)) next()

      ## validation PGS files (same idea: bc vs original, same valid_pop/group)
      ## NOTE: metric="spearman" -> valid_spearman/
      PGS0_file <- file.path(BC_PGS_BASE,
                             paste0("valid_", metric),
                             paste0(PHENO_LOWER, "_", valid_pop, "_pgs_valid_", metric, ".", thresh0, ".sscore"))
      PGS1_file <- file.path(ORIG_PGS_BASE,
                             paste0("valid_", metric),
                             paste0(PHENO_LOWER, "_", valid_pop, "_pgs_valid_", metric, ".", thresh1, ".sscore"))

      if (!file.exists(PGS0_file)) stop("Missing bc sscore: ", PGS0_file)
      if (!file.exists(PGS1_file)) stop("Missing orig sscore: ", PGS1_file)


      ## universe IDs from bc file
      PGS_IDs <- read.table(PGS0_file, header=FALSE, col.names=PGS_COLS)[,1:2]

      covs <- covs_table[covs_table$FID %in% PGS_IDs$FID, ]
      if (group == "female") {
        covs <- covs[covs[["X31.0.0"]] == 1, ]
      } else if (group == "male") {
        covs <- covs[covs[["X31.0.0"]] == 0, ]
      }

      pheno <- pheno_table[pheno_table$FID %in% PGS_IDs$FID, ]
      pheno_cov <- merge(pheno, covs, by=c("FID","IID"))
      colnames(pheno_cov)[3] <- "pheno_code"

      y <- pheno_cov$pheno_code
      X <- pheno_cov[, !colnames(pheno_cov) %in% c("FID","IID","pheno_code")]

      prs0 <- read.table(PGS0_file, header=FALSE, col.names=PGS_COLS)[,c(1,6)]
      prs1 <- read.table(PGS1_file, header=FALSE, col.names=PGS_COLS)[,c(1,6)]

      ## bc PRS: inverse Box-Cox using trait lambda
      prs0 <- inv_boxcox_prs(merge(pheno_cov, prs0, by="FID")$SCORE1_SUM, lambda)

      ## original PRS: raw score
      prs1 <- merge(pheno_cov, prs1, by="FID")$SCORE1_SUM

      ## metric-specific predictive performance
      if (metric == "std") {
        if (R2_TEST == "SPEARMAN") {
          prs.R2 <- r2_diff_spearman(prs0, prs1, y=y, X=X, B=B_BOOT)
        } else {
          prs.R2 <- r2_diff_boot(prs0, prs1, y=y, X=X, B=B_BOOT)
        }
      } else if (metric == "spearman") {
        prs.R2 <- r2_diff_spearman(prs0, prs1, y=y, X=X, B=B_BOOT)
      } else if (metric == "rankLin") {
        prs.R2 <- r2_diff_rankLin(prs0, prs1, y=y, X=X, B=B_BOOT)
      } else if (metric == "rankEC") {
        prs.R2 <- r2_diff_rankEC(prs0, prs1, y=y, X=X, B=B_BOOT)
      } else {
        stop("Unknown metric: ", metric)
      }

      prs.result <- data.frame(
        metric = metric,
        thresh0 = thresh0,
        thresh1 = thresh1,
        r20 = prs.R2$r20,
        r21 = prs.R2$r21,
        ci_95 = prs.R2$ci_95,
        pval = prs.R2$pval
      )

      write.table(prs.result, file=out_file, row.names=FALSE, quote=FALSE)
      cat("Wrote: ", out_file, "\n", sep="")
    }
  }
}
