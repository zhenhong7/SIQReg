#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(conquer)
})

## inputs
pheno_name <- "LDL"
tissue     <- "Whole_Blood"
chr        <- "1"

## replace with the real ENSG IDs used in your predict file
target_genes <- c(
  "ENSG00000169174.10",  # PCSK9
  "ENSG00000134222.16"   # PSRC1
)

## paths: configurable via environment variables
pheno_dir   <- Sys.getenv("PHENO_DIR",       "phenotypes")
covar_fp    <- Sys.getenv("COVAR_FILE",      "phenotypes/covariates_sa40PC.pheno")
expr_dir    <- Sys.getenv("EXPR_DIR",        "TWAS/imputed_expr")
expr_file   <- sprintf("%s/%s/%s__chr%s__predict.txt", expr_dir, tissue, tissue, chr)
keep_file   <- Sys.getenv("KEEP_FILE",       "ancestry_ids/whitebrit.fid_iid.txt")
statin_fp   <- Sys.getenv("STATIN_FILE",     "phenotypes/Statins/Statins.pheno")
meta_both   <- Sys.getenv("META_BOTH_FILE",  "Lambda_learning/meta_both_summary_all.txt")

BIRTH_YEAR_BASELINE <- as.integer(Sys.getenv("BIRTH_YEAR_BASELINE", "2025"))

out_dir <- file.path(Sys.getenv("TWAS_RESULTS", "results/TWAS"),
                     "statin_stratified", pheno_name, tissue)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## helpers
boxcox_transform <- function(y, lambda) {
  if (abs(lambda) < 1e-8) log(y) else (y^lambda - 1) / lambda
}

add_cols <- function(prefix, v) {
  setNames(data.frame(t(v)), paste0(prefix, "_", names(v)))
}

run_one <- function(dat, gene_id, expr_col = "expr", scale_name = c("raw","bc"), lambda_hat = NA_real_, statin_group = NA_integer_) {
  scale_name <- match.arg(scale_name)

  dat <- copy(dat)

  if (scale_name == "raw") {
    dat[, Y_new := Y]
    lambda_out <- NA_real_
  } else {
    dat <- dat[is.finite(Y) & Y > 0]
    dat[, Y_new := boxcox_transform(Y, lambda_hat)]
    lambda_out <- lambda_hat
  }

  if (nrow(dat) == 0) return(NULL)

  fm_full <- as.formula(paste0("Y ~ ", expr_col, " + sex + age + age2 + ", paste0("PC", 1:40, collapse = " + ")))
  fm_lr   <- update(fm_full, Y_new ~ .)

  ## LR
  lr <- lm(fm_lr, dat)
  sum_lr <- summary(lr)
  coeffs <- sum_lr$coefficients
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
    lr_ci_g_low  <- if (!is.null(ci_lr) && (expr_col %in% rownames(ci_lr))) ci_lr[expr_col, 1] else NA_real_
    lr_ci_g_high <- if (!is.null(ci_lr) && (expr_col %in% rownames(ci_lr))) ci_lr[expr_col, 2] else NA_real_
  }

  lr_beta_age <- coeffs["age", "Estimate"]
  lr_se_age   <- coeffs["age", "Std. Error"]
  lr_p_age    <- coeffs["age", "Pr(>|t|)"]
  lr_ci_a_low <- if (!is.null(ci_lr) && ("age" %in% rownames(ci_lr))) ci_lr["age", 1] else NA_real_
  lr_ci_a_hi  <- if (!is.null(ci_lr) && ("age" %in% rownames(ci_lr))) ci_lr["age", 2] else NA_real_

  lr_beta_sex <- coeffs["sex", "Estimate"]
  lr_se_sex   <- coeffs["sex", "Std. Error"]
  lr_p_sex    <- coeffs["sex", "Pr(>|t|)"]
  lr_ci_s_low <- if (!is.null(ci_lr) && ("sex" %in% rownames(ci_lr))) ci_lr["sex", 1] else NA_real_
  lr_ci_s_hi  <- if (!is.null(ci_lr) && ("sex" %in% rownames(ci_lr))) ci_lr["sex", 2] else NA_real_

  ## QR: exact same setup
  quant_levels <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
  nQuant <- length(quant_levels)
  X_qr <- model.matrix(fm_lr, dat)[, -1, drop = FALSE]
  Yv   <- dat$Y_new

  beta_g  <- setNames(rep(NA_real_, nQuant), quant_levels)
  se_g    <- beta_g
  p_g     <- beta_g
  ci_g_lo <- beta_g
  ci_g_hi <- beta_g

  beta_a  <- beta_g
  se_a    <- beta_g
  p_a     <- beta_g
  ci_a_lo <- beta_g
  ci_a_hi <- beta_g

  beta_s  <- beta_g
  se_s    <- beta_g
  p_s     <- beta_g
  ci_s_lo <- beta_g
  ci_s_hi <- beta_g

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
      cn <- colnames(X_qr)
      fit$coeff <- setNames(fit$coeff, c("(Intercept)", cn))
      rownames(fit$asyCI) <- c("(Intercept)", cn)

      beta_g[i]  <- fit$coeff[expr_col]
      ci_g_lo[i] <- fit$asyCI[expr_col, 1]
      ci_g_hi[i] <- fit$asyCI[expr_col, 2]
      se_g[i]    <- (fit$asyCI[expr_col, 2] - fit$coeff[expr_col]) / qnorm(0.975)
      p_g[i]     <- 2 * pnorm(-abs(fit$coeff[expr_col] / se_g[i]))

      beta_a[i]  <- fit$coeff["age"]
      ci_a_lo[i] <- fit$asyCI["age", 1]
      ci_a_hi[i] <- fit$asyCI["age", 2]
      se_a[i]    <- (fit$asyCI["age", 2] - fit$coeff["age"]) / qnorm(0.975)
      p_a[i]     <- 2 * pnorm(-abs(fit$coeff["age"] / se_a[i]))

      beta_s[i]  <- fit$coeff["sex"]
      ci_s_lo[i] <- fit$asyCI["sex", 1]
      ci_s_hi[i] <- fit$asyCI["sex", 2]
      se_s[i]    <- (fit$asyCI["sex", 2] - fit$coeff["sex"]) / qnorm(0.975)
      p_s[i]     <- 2 * pnorm(-abs(fit$coeff["sex"] / se_s[i]))
    }
  }

  ## same bootstrap covariance
  n_boot <- 100
  boot_coefs <- matrix(NA_real_, nrow = n_boot, ncol = nQuant)

  if (!skip_qr) {
    for (b in seq_len(n_boot)) {
      db <- dat[sample(nrow(dat), replace = TRUE)]
      Xb <- model.matrix(fm_lr, db)[, -1, drop = FALSE]
      Yb <- db$Y_new
      boot_coefs[b, ] <- sapply(quant_levels, function(tau) {
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
    R <- diag(nQuant)[-nQuant, ] - diag(nQuant)[-1, ]
    RsRt <- R %*% cov_mat %*% t(R)
    if (all(is.finite(RsRt)) && det(RsRt) > 0) {
      W2 <- t(R %*% orig) %*% solve(RsRt) %*% (R %*% orig)
      equality_p <- pchisq(W2, df = nQuant - 1, lower.tail = FALSE)
    }
  }

  gene_row <- cbind(
    data.frame(
      chr = chr,
      gene = gene_id,
      statins = statin_group,
      scale = scale_name,
      n = nrow(dat),
      lambda_hat = lambda_out,
      beta_lr = lr_beta_gene,
      se_lr = lr_se_gene,
      p_lr = lr_p_gene,
      ci_lr_lower = lr_ci_g_low,
      ci_lr_upper = lr_ci_g_high,
      joint_p_qr = joint_p,
      equality_p_qr = equality_p
    ),
    add_cols("beta_qr", setNames(beta_g, quant_levels)),
    add_cols("se_qr", setNames(se_g, quant_levels)),
    add_cols("ci_qr_lower", setNames(ci_g_lo, quant_levels)),
    add_cols("ci_qr_upper", setNames(ci_g_hi, quant_levels)),
    add_cols("p_qr", setNames(p_g, quant_levels))
  )

  age_sex_row <- cbind(
    data.frame(
      chr = chr,
      gene = gene_id,
      statins = statin_group,
      scale = scale_name,
      n = nrow(dat),
      beta_lr_age = lr_beta_age,
      se_lr_age = lr_se_age,
      p_lr_age = lr_p_age,
      ci_lr_age_lower = lr_ci_a_low,
      ci_lr_age_upper = lr_ci_a_hi,
      beta_lr_sex = lr_beta_sex,
      se_lr_sex = lr_se_sex,
      p_lr_sex = lr_p_sex,
      ci_lr_sex_lower = lr_ci_s_low,
      ci_lr_sex_upper = lr_ci_s_hi
    ),
    add_cols("beta_qr_age", setNames(beta_a, quant_levels)),
    add_cols("se_qr_age", setNames(se_a, quant_levels)),
    add_cols("ci_qr_age_lower", setNames(ci_a_lo, quant_levels)),
    add_cols("ci_qr_age_upper", setNames(ci_a_hi, quant_levels)),
    add_cols("p_qr_age", setNames(p_a, quant_levels)),
    add_cols("beta_qr_sex", setNames(beta_s, quant_levels)),
    add_cols("se_qr_sex", setNames(se_s, quant_levels)),
    add_cols("ci_qr_sex_lower", setNames(ci_s_lo, quant_levels)),
    add_cols("ci_qr_sex_upper", setNames(ci_s_hi, quant_levels)),
    add_cols("p_qr_sex", setNames(p_s, quant_levels))
  )

  list(gene_row = gene_row, age_sex_row = age_sex_row)
}

## load shared data
keep_iid <- as.character(fread(keep_file, header = FALSE)[[2]])

phfile <- list.files(pheno_dir, pattern = paste0("^", pheno_name, ".*\\.pheno$"), full.names = TRUE)[1]
if (is.na(phfile)) stop("Phenotype file not found for: ", pheno_name)
ph <- fread(phfile)
trait_col <- setdiff(names(ph), c("FID", "IID"))[1]
setnames(ph, c("FID", "IID", trait_col), c("FID", "IID", "Y"))
ph[, IID := as.character(IID)]

cv <- fread(covar_fp)
setnames(cv, names(cv)[3:44], c("sex", "age", paste0("PC", 1:40)))
cv[, IID := as.character(IID)]

st <- fread(statin_fp)
st[, IID := as.character(IID)]

meta_tab <- fread(meta_both)
lambda_hat <- meta_tab[phenotype == pheno_name, lambda_WM][1]
if (!is.finite(lambda_hat)) stop("lambda_WM not found / not finite for phenotype: ", pheno_name)

base_data <- merge(ph, cv, by = "IID")[!is.na(Y)]
base_data <- base_data[IID %in% keep_iid]
base_data <- merge(base_data, st[, .(IID, statins)], by = "IID", all = FALSE)

base_data[, age := BIRTH_YEAR_BASELINE - age]
mu_age <- mean(base_data$age, na.rm = TRUE)
base_data[, age := age - mu_age]
base_data[, age2 := age^2]

expr_header <- names(fread(expr_file, nrows = 0))
need_cols <- c("IID", target_genes)
miss <- setdiff(need_cols, expr_header)
if (length(miss) > 0) stop("Missing genes in expr file: ", paste(miss, collapse = ", "))

expr_dt <- fread(expr_file, select = need_cols)
expr_dt[, IID := as.character(sub("_.*", "", IID))]

dat0 <- merge(base_data, expr_dt, by = "IID", all = FALSE)
if (nrow(dat0) == 0) stop("No merged samples")

## run: 2 genes x 2 statin groups x 2 scales
all_gene_rows <- list()
all_age_sex_rows <- list()
k1 <- 1
k2 <- 1

for (gene_id in target_genes) {
  dat_gene <- copy(dat0)
  setnames(dat_gene, gene_id, "expr")

  for (g in c(0, 1)) {
    dg <- dat_gene[statins == g]
    if (nrow(dg) == 0) next

    out_raw <- run_one(dg, gene_id = gene_id, scale_name = "raw", lambda_hat = lambda_hat, statin_group = g)
    out_bc  <- run_one(dg, gene_id = gene_id, scale_name = "bc",  lambda_hat = lambda_hat, statin_group = g)

    if (!is.null(out_raw)) {
      all_gene_rows[[k1]] <- out_raw$gene_row; k1 <- k1 + 1
      all_age_sex_rows[[k2]] <- out_raw$age_sex_row; k2 <- k2 + 1
    }
    if (!is.null(out_bc)) {
      all_gene_rows[[k1]] <- out_bc$gene_row; k1 <- k1 + 1
      all_age_sex_rows[[k2]] <- out_bc$age_sex_row; k2 <- k2 + 1
    }
  }
}

gene_tab <- rbindlist(all_gene_rows, fill = TRUE)
age_sex_tab <- rbindlist(all_age_sex_rows, fill = TRUE)

fwrite(
  gene_tab,
  file = file.path(out_dir, sprintf("association_results_statin_stratified_chr%s.txt", chr)),
  sep = "\t"
)

fwrite(
  age_sex_tab,
  file = file.path(out_dir, sprintf("association_age_sex_statin_stratified_chr%s.txt", chr)),
  sep = "\t"
)

print(gene_tab[, .(gene, statins, scale, n, beta_lr, p_lr, joint_p_qr, equality_p_qr)])
