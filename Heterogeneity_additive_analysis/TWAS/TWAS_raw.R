#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript TWAS_raw.R <chr> <gene_index> <pheno> <tissue>")
chr        <- args[1]
gene_index <- as.integer(args[2])
pheno_name <- args[3]
tissue     <- args[4]

suppressPackageStartupMessages({
  library(data.table)
  library(conquer)
})

pheno_dir <- Sys.getenv("PHENO_DIR",  "phenotypes")
covar_fp  <- Sys.getenv("COVAR_FILE", "phenotypes/covariates_sa40PC.pheno")
expr_dir  <- Sys.getenv("EXPR_DIR",   "TWAS/imputed_expr")
expr_file <- sprintf("%s/%s/%s__chr%s__predict.txt", expr_dir, tissue, tissue, chr)
out_dir   <- file.path(Sys.getenv("TWAS_RESULTS", "results/TWAS"), pheno_name, tissue)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BIRTH_YEAR_BASELINE <- as.integer(Sys.getenv("BIRTH_YEAR_BASELINE", "2025"))

keep_file <- Sys.getenv("KEEP_FILE", "ancestry_ids/whitebrit.fid_iid.txt")
keep_iid  <- as.character(fread(keep_file, header = FALSE)[[2]])

# phenotype
phfile <- list.files(pheno_dir, pattern = paste0("^", pheno_name, ".*\\.pheno$"), full.names = TRUE)[1]
if (is.na(phfile)) stop("Phenotype file not found for: ", pheno_name)
ph <- fread(phfile)
trait_col <- setdiff(names(ph), c("FID","IID"))[1]
setnames(ph, c("FID","IID",trait_col), c("FID","IID","Y"))
ph[, IID := as.character(IID)]

# covariates
cv <- fread(covar_fp)
setnames(cv, names(cv)[3:44], c("sex","age", paste0("PC",1:40)))
cv[, IID := as.character(IID)]

base_data <- merge(ph, cv, by = "IID")[!is.na(Y)]
base_data <- base_data[IID %in% keep_iid]

# gene list
genes <- grep("^ENSG", names(fread(expr_file, nrows = 0)), value = TRUE)
if (gene_index < 1 || gene_index > length(genes)) stop("gene_index out of range")
gene_id  <- genes[gene_index]
expr_col <- "expr"

# expression
expr_dt <- fread(expr_file, select = c("IID", gene_id))
setnames(expr_dt, gene_id, expr_col)
expr_dt[, IID := as.character(sub("_.*", "", IID))]
dat <- merge(base_data, expr_dt, by = "IID", all = FALSE)
if (nrow(dat) == 0) quit(status = 1)

dat[, age := BIRTH_YEAR_BASELINE - age]
mu_age <- mean(dat$age, na.rm = TRUE)
dat[, age := age - mu_age]
dat[, age2 := age^2]

fm_full <- as.formula(paste0("Y ~ ", expr_col, " + sex + age + age2 + ", paste0("PC", 1:40, collapse = " + ")))
fm_lr   <- update(fm_full, Y_new ~ .)
X_full  <- model.matrix(fm_full, dat)

lambda_hat <- NA_real_
dat[, Y_new := Y]

# LR + SE
lr <- lm(fm_lr, dat)
sum_lr <- summary(lr); coeffs <- sum_lr$coefficients
ci_lr  <- tryCatch(confint(lr), error = function(e) NULL)

if (!(expr_col %in% rownames(coeffs))) {
  lr_beta_gene <- NA_real_
  lr_se_gene   <- NA_real_
  lr_p_gene    <- NA_real_
  lr_ci_g_low  <- NA_real_
  lr_ci_g_high <- NA_real_
} else {
  lr_beta_gene <- coeffs[expr_col, "Estimate"]
  lr_se_gene   <- coeffs[expr_col, "Std. Error"]
  lr_p_gene    <- coeffs[expr_col, "Pr(>|t|)"]
  lr_ci_g_low  <- if (!is.null(ci_lr) && (expr_col %in% rownames(ci_lr))) ci_lr[expr_col,1] else NA_real_
  lr_ci_g_high <- if (!is.null(ci_lr) && (expr_col %in% rownames(ci_lr))) ci_lr[expr_col,2] else NA_real_
}

lr_beta_age <- coeffs["age","Estimate"]; lr_se_age <- coeffs["age","Std. Error"]; lr_p_age <- coeffs["age","Pr(>|t|)"]
lr_ci_a_low <- if (!is.null(ci_lr) && ("age" %in% rownames(ci_lr))) ci_lr["age",1] else NA_real_
lr_ci_a_hi  <- if (!is.null(ci_lr) && ("age" %in% rownames(ci_lr))) ci_lr["age",2] else NA_real_

lr_beta_sex <- coeffs["sex","Estimate"]; lr_se_sex <- coeffs["sex","Std. Error"]; lr_p_sex <- coeffs["sex","Pr(>|t|)"]
lr_ci_s_low <- if (!is.null(ci_lr) && ("sex" %in% rownames(ci_lr))) ci_lr["sex",1] else NA_real_
lr_ci_s_hi  <- if (!is.null(ci_lr) && ("sex" %in% rownames(ci_lr))) ci_lr["sex",2] else NA_real_

# QR
quant_levels <- c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99); nQuant <- length(quant_levels)
X_qr <- model.matrix(fm_lr, dat)[,-1,drop=FALSE]; Yv <- dat$Y_new

beta_g  <- setNames(rep(NA_real_, nQuant), quant_levels)
se_g    <- beta_g; p_g <- beta_g
ci_g_lo <- beta_g; ci_g_hi <- beta_g

beta_a  <- beta_g; se_a <- beta_g; p_a <- beta_g; ci_a_lo <- beta_g; ci_a_hi <- beta_g
beta_s  <- beta_g; se_s <- beta_g; p_s <- beta_g; ci_s_lo <- beta_g; ci_s_hi <- beta_g

is_const <- apply(X_qr, 2, function(z) {
  z <- z[is.finite(z)]
  length(z) == 0L || length(unique(z)) <= 1L
})
skip_qr <- any(is_const)

orig <- if (skip_qr) rep(NA_real_, nQuant) else sapply(quant_levels, function(tau) {
  fit <- conquer(X = X_qr, Y = Yv, tau = tau, ci = "asymptotic")
  if (length(fit$coeff) == 0L) return(NA_real_)
  fit$coeff <- setNames(fit$coeff, c("(Intercept)", colnames(X_qr)))
  fit$coeff[expr_col]
})

if (!skip_qr) {
  for (i in seq_along(quant_levels)) {
    tau <- quant_levels[i]
    fit <- conquer(X = X_qr, Y = Yv, tau = tau, ci = "asymptotic")
    if (length(fit$coeff) == 0L || is.null(fit$asyCI) || length(fit$asyCI) == 0L) next
    cn  <- colnames(X_qr)
    fit$coeff <- setNames(fit$coeff, c("(Intercept)", cn))
    rownames(fit$asyCI) <- c("(Intercept)", cn)

    # gene
    beta_g[i]  <- fit$coeff[expr_col]
    ci_g_lo[i] <- fit$asyCI[expr_col,1]; ci_g_hi[i] <- fit$asyCI[expr_col,2]
    se_g[i]    <- (fit$asyCI[expr_col,2] - fit$coeff[expr_col]) / qnorm(0.975)
    p_g[i]     <- 2 * pnorm(-abs(fit$coeff[expr_col] / se_g[i]))

    # age
    beta_a[i]  <- fit$coeff["age"]
    ci_a_lo[i] <- fit$asyCI["age",1];  ci_a_hi[i] <- fit$asyCI["age",2]
    se_a[i]    <- (fit$asyCI["age",2] - fit$coeff["age"]) / qnorm(0.975)
    p_a[i]     <- 2 * pnorm(-abs(fit$coeff["age"] / se_a[i]))

    # sex
    beta_s[i]  <- fit$coeff["sex"]
    ci_s_lo[i] <- fit$asyCI["sex",1];  ci_s_hi[i] <- fit$asyCI["sex",2]
    se_s[i]    <- (fit$asyCI["sex",2] - fit$coeff["sex"]) / qnorm(0.975)
    p_s[i]     <- 2 * pnorm(-abs(fit$coeff["sex"] / se_s[i]))
  }
}

# bootstrap covariance
n_boot <- 100
boot_coefs <- matrix(NA_real_, nrow = n_boot, ncol = nQuant)
if (!skip_qr) {
  for (b in seq_len(n_boot)) {
    db <- dat[sample(nrow(dat), replace = TRUE)]
    Xb <- model.matrix(fm_lr, db)[,-1,drop=FALSE]; Yb <- db$Y_new
    boot_coefs[b,] <- sapply(quant_levels, function(tau){
      fit <- conquer(X = Xb, Y = Yb, tau = tau, ci = "asymptotic")
      if (length(fit$coeff) == 0L) return(NA_real_)
      fit$coeff <- setNames(fit$coeff, c("(Intercept)", colnames(Xb)))
      fit$coeff[expr_col]
    })
  }
}

cov_mat <- tryCatch(
  cov(boot_coefs, use = "pairwise.complete.obs"),
  error = function(e) matrix(NA_real_, nQuant, nQuant)
)

joint_p <- equality_p <- NA_real_
ok_cov  <- (!skip_qr) && is.matrix(cov_mat) && all(dim(cov_mat) == c(nQuant, nQuant)) && all(is.finite(cov_mat))
ok_orig <- all(is.finite(orig))

if (ok_cov && ok_orig) {
  if (det(cov_mat) > 0) {
    W <- t(orig) %*% solve(cov_mat) %*% orig
    joint_p <- pchisq(W, df = nQuant, lower.tail = FALSE)
  }
  R <- diag(nQuant)[-nQuant,] - diag(nQuant)[-1,]
  RsRt <- R %*% cov_mat %*% t(R)
  if (all(is.finite(RsRt)) && det(RsRt) > 0) {
    W2 <- t(R %*% orig) %*% solve(RsRt) %*% (R %*% orig)
    equality_p <- pchisq(W2, df = nQuant-1, lower.tail = FALSE)
  }
}

add_cols <- function(prefix, v) setNames(data.frame(t(v)), paste0(prefix, "_", names(v)))

gene_row <- cbind(
  data.frame(chr=chr, gene=gene_id, lambda_hat=lambda_hat,
             beta_lr=lr_beta_gene, se_lr=lr_se_gene, p_lr=lr_p_gene,
             ci_lr_lower=lr_ci_g_low, ci_lr_upper=lr_ci_g_high,
             joint_p_qr=joint_p, equality_p_qr=equality_p),
  add_cols("beta_qr", setNames(beta_g,  quant_levels)),
  add_cols("se_qr",   setNames(se_g,    quant_levels)),
  add_cols("ci_qr_lower", setNames(ci_g_lo, quant_levels)),
  add_cols("ci_qr_upper", setNames(ci_g_hi, quant_levels)),
  add_cols("p_qr",    setNames(p_g,    quant_levels))
)

age_sex_row <- cbind(
  data.frame(chr=chr, gene=gene_id,
             beta_lr_age=lr_beta_age, se_lr_age=lr_se_age, p_lr_age=lr_p_age,
             ci_lr_age_lower=lr_ci_a_low, ci_lr_age_upper=lr_ci_a_hi,
             beta_lr_sex=lr_beta_sex, se_lr_sex=lr_se_sex, p_lr_sex=lr_p_sex,
             ci_lr_sex_lower=lr_ci_s_low, ci_lr_sex_upper=lr_ci_s_hi),
  add_cols("beta_qr_age", setNames(beta_a, quant_levels)),
  add_cols("se_qr_age",   setNames(se_a,   quant_levels)),
  add_cols("ci_qr_age_lower", setNames(ci_a_lo, quant_levels)),
  add_cols("ci_qr_age_upper", setNames(ci_a_hi, quant_levels)),
  add_cols("p_qr_age",    setNames(p_a,   quant_levels)),
  add_cols("beta_qr_sex", setNames(beta_s, quant_levels)),
  add_cols("se_qr_sex",   setNames(se_s,   quant_levels)),
  add_cols("ci_qr_sex_lower", setNames(ci_s_lo, quant_levels)),
  add_cols("ci_qr_sex_upper", setNames(ci_s_hi, quant_levels)),
  add_cols("p_qr_sex",    setNames(p_s,   quant_levels))
)

fwrite(gene_row, file = file.path(out_dir, sprintf("association_results_raw_chr%s.txt", chr)),
       sep = "\t", append = file.exists(file.path(out_dir, sprintf("association_results_raw_chr%s.txt", chr))),
       col.names = !file.exists(file.path(out_dir, sprintf("association_results_raw_chr%s.txt", chr))))
fwrite(age_sex_row, file = file.path(out_dir, sprintf("association_age_sex_raw_chr%s.txt", chr)),
       sep = "\t", append = file.exists(file.path(out_dir, sprintf("association_age_sex_raw_chr%s.txt", chr))),
       col.names = !file.exists(file.path(out_dir, sprintf("association_age_sex_raw_chr%s.txt", chr))))

# NOTE: intentionally NOT saving all_cov_matrices_raw_chr*.RData to avoid parallel load/save corruption
