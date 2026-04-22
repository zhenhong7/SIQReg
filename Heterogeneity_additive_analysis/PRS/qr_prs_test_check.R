#!/usr/bin/env Rscript
# QR heterogeneity tests for PRS (raw vs Box-Cox)
# Usage: Rscript qr_prs_test_check.R <PHENO_ID> <PRS_ID> <MODE>

suppressPackageStartupMessages({
  library(conquer)
  library(data.table)
  library(MASS)      # ginv
  library(Matrix)    # bdiag
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript qr_prs_test_check.R <PHENO_ID> <PRS_ID> <MODE{mean|median}>")
pheno_id <- args[1]
prs_id   <- args[2]    # basename without .pheno
mode_lambda <- match.arg(args[3], c("mean","median"))

pheno_dir  <- Sys.getenv("PHENO_DIR",       "phenotypes")
covar_file <- Sys.getenv("COVAR_FILE",      "phenotypes/covariates_sa40PC.pheno")
prs_dir    <- Sys.getenv("PRS_DIR",         "PRS/ALL_PRS")

keep_file <- Sys.getenv("KEEP_FILE", "ancestry_ids/whitebrit.fid_iid.txt")
if (!file.exists(keep_file)) stop("keep_file not found: ", keep_file)
keep_iid <- as.character(fread(keep_file, header = FALSE)[[2]])

# lambda source is the meta-both summary (mean vs median)
meta_both_file <- Sys.getenv("META_BOTH_FILE", "Lambda_learning/meta_both_summary_all.txt")

# output base depends on MODE; everything else unchanged
out_base <- file.path(Sys.getenv("OUT_BASE", "results/qr_prs_tests"),
                      if (mode_lambda == "mean") "mean_meta" else "median_meta")
out_dir <- file.path(out_base, pheno_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_txt <- file.path(out_dir, paste0(pheno_id, "_qr_prs_tests.txt"))

taus <- c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
K    <- length(taus)
p    <- 1     # only PRS is the target
d    <- p * K
B    <- 100
set.seed(2025)
write_qr_per_tau <- TRUE

# Age handling toggle: set TRUE if 'age' column is actually YOB
age_is_yob <- FALSE
BIRTH_YEAR_BASELINE <- as.integer(Sys.getenv("BIRTH_YEAR_BASELINE", "2025"))

tau_label <- function(t) paste0("t", sprintf("%03d", as.integer(round(1000*t))))
labs <- sapply(taus, tau_label)
alpha_ne <- 0.10; band_lo <- 1/1.2; band_hi <- 1.2   # non-equivalence band

pheno_path <- list.files(pheno_dir, pattern = paste0("^", pheno_id, ".*\\.pheno$"), full.names = TRUE)
if (length(pheno_path) != 1) stop("Found ", length(pheno_path), " phenotype files for ", pheno_id)
pheno <- fread(pheno_path)
setnames(pheno, c("FID","IID","Pheno"))
pheno[, IID := as.character(IID)]

covar <- fread(covar_file)
setnames(covar, names(covar)[3:44], c("sex","age", paste0("PC",1:40)))
covar[, IID := as.character(IID)]

base_data <- merge(pheno, covar, by = "IID")[!is.na(Pheno)]
rm(pheno, covar); invisible(gc())
setDT(base_data)
base_data[, IID := as.integer(IID)]

prs_file <- file.path(prs_dir, paste0(prs_id, ".pheno"))
if (!file.exists(prs_file)) stop("PRS file not found: ", prs_file)
prs <- fread(prs_file)
setnames(prs, names(prs)[3], "PRS")
prs[, IID := as.integer(IID)]

base_data <- merge(prs, base_data, by = "IID")
rm(prs); invisible(gc())

base_data <- base_data[IID %in% keep_iid]

if (age_is_yob) {
  base_data[, age := BIRTH_YEAR_BASELINE - age]
}
mu_age <- mean(base_data$age, na.rm = TRUE)
base_data[, age := age - mu_age]
base_data[, age2 := age^2]

dat <- base_data

fm_full <- as.formula(paste0(
  "Pheno ~ PRS + sex + age + age2 + ", paste0("PC", 1:40, collapse = " + ")
))

mf_raw <- model.frame(fm_full, data = dat, na.action = na.omit)
X <- model.matrix(fm_full, data = mf_raw)[, -1, drop = FALSE]  # no intercept in X
Y <- model.response(mf_raw)
if (nrow(X) != length(Y)) stop("X/Y mismatch after model.frame: ", nrow(X), " vs ", length(Y))

cat("PRS target:", prs_id, "   MODE:", mode_lambda, "\n")
cat("N(raw after NA drop) =", nrow(X), "\n")

is_const <- vapply(seq_len(ncol(X)), function(j) var(X[, j]) == 0, logical(1))
if (any(is_const)) {
  cat("Dropping constant columns in X:", paste(colnames(X)[is_const], collapse = ", "), "\n")
  X <- X[, !is_const, drop = FALSE]
}

qrX <- qr(X)
if (qrX$rank < ncol(X)) {
  keep_idx <- sort(qrX$pivot[seq_len(qrX$rank)])
  dropped  <- setdiff(seq_len(ncol(X)), keep_idx)
  dropped_names <- colnames(X)[dropped]
  cat("Dropping aliased columns in X:", paste(dropped_names, collapse = ", "), "\n")

  if ("PRS" %in% dropped_names) {
    cat("PRS is aliased/constant; writing NA rows and exiting.\n")

    na_cols <- as.list(setNames(rep(NA_real_, 0), character(0)))
    for (lab in labs) {
      na_cols[[paste0("beta_qr_PRS_", lab)]]         <- NA_real_
      na_cols[[paste0("se_qr_PRS_",   lab)]]         <- NA_real_
      na_cols[[paste0("p_qr_PRS_",    lab)]]         <- NA_real_
      na_cols[[paste0("ci_qr_PRS_lower_", lab)]]     <- NA_real_
      na_cols[[paste0("ci_qr_PRS_upper_", lab)]]     <- NA_real_
      na_cols[[paste0("ratio_qr_over_ols_", lab)]]   <- NA_real_
      na_cols[[paste0("ratio_ci_lower_", lab)]]      <- NA_real_
      na_cols[[paste0("ratio_ci_upper_", lab)]]      <- NA_real_
      na_cols[[paste0("non_equiv_flag_", lab)]]      <- NA_real_
    }

    row_raw <- data.table(
      phenotype = pheno_id, prs_id = prs_id, scale = "raw",
      lambda_used = NA_real_,
      W_joint = NA_real_, p_joint = NA_real_,
      W_equality = NA_real_, p_equality = NA_real_,
      lm_beta_PRS = NA_real_, lm_se_PRS = NA_real_, lm_p_PRS = NA_real_,
      lm_ci_lower_PRS = NA_real_, lm_ci_upper_PRS = NA_real_,
      non_equiv_count = NA_real_
    )
    row_raw <- cbind(row_raw, as.data.table(na_cols), fill = TRUE)

    row_tr <- copy(row_raw); row_tr[, `:=`(scale = "transformed", lambda_used = NA_real_)]

    res_na <- rbind(row_raw, row_tr, fill = TRUE)
    add_header <- !file.exists(out_txt)
    fwrite(res_na, out_txt, sep = "\t", append = !add_header, col.names = add_header)
    quit(save = "no", status = 0)
  }

  X <- X[, keep_idx, drop = FALSE]
}

lambda_used <- NA_real_
meta_tab <- fread(meta_both_file)
if (!all(c("phenotype","lambda_DL","lambda_WM") %in% names(meta_tab))) {
  stop("meta_both_summary_all.txt must contain: phenotype, lambda_DL, lambda_WM")
}
lambda_used <- if (mode_lambda == "mean") {
  meta_tab[phenotype == pheno_id, lambda_DL][1]
} else {
  meta_tab[phenotype == pheno_id, lambda_WM][1]
}
cat("[", mode_lambda, "] lambda_used = ", lambda_used, "\n", sep = "")

qr_fit_one <- function(X, Y, tau, ci = "none") {
  fit <- conquer(X = X, Y = Y, tau = tau, ci = ci, tol = 1e-6)
  nm <- c("(Intercept)", colnames(X))
  setNames(fit$coeff, nm)
}

theta_from_data <- function(X, Y, taus, S_names) {
  betas <- lapply(taus, function(tau) {
    cf <- qr_fit_one(X, Y, tau, ci = "none")
    unname(cf[S_names])
  })
  as.numeric(do.call(c, betas))  # length p*K
}

boot_cov_theta <- function(X, Y, taus, S_names, B, min_valid = 25) {
  n <- nrow(X); d <- length(S_names) * length(taus)
  Theta <- matrix(NA_real_, nrow = 0, ncol = d)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    tb <- try(theta_from_data(X[ii,,drop=FALSE], Y[ii], taus, S_names), silent = TRUE)
    if (!inherits(tb, "try-error") && all(is.finite(tb))) {
      Theta <- rbind(Theta, tb)
    }
  }
  if (nrow(Theta) < min_valid) return(NULL)
  S <- cov(Theta)
  if (qr(S)$rank < ncol(S)) return(NULL)
  S
}

wald_quad <- function(theta, Sigma) {
  invS <- tryCatch(solve(Sigma), error = function(e) MASS::ginv(Sigma))
  as.numeric(t(theta) %*% invS %*% theta)
}

build_R <- function(K, p) {
  D <- diag(K-1, K-1)
  D <- cbind(D, matrix(0, nrow = K-1, ncol = 1)) - cbind(matrix(0, K-1, 1), D)
  as.matrix(Matrix::bdiag(replicate(p, D, simplify = FALSE)))
}

qr_per_tau_stats_prs <- function(X, Y, taus) {
  out <- list()
  for (tau in taus) {
    lab <- tau_label(tau)
    res <- try({
      fit <- conquer(X = X, Y = Y, tau = tau, ci = "asymptotic", tol = 1e-6)
      names(fit$coeff) <- c("(Intercept)", colnames(X))
      if (!is.null(fit$asyCI) && length(fit$asyCI) > 0L) {
        rownames(fit$asyCI) <- c("(Intercept)", colnames(X))
      }
      v <- "PRS"
      if (v %in% names(fit$coeff)) {
        est <- fit$coeff[v]
        if (!is.null(fit$asyCI) && v %in% rownames(fit$asyCI)) {
          se  <- (fit$asyCI[v, 2] - est) / qnorm(0.975)
          z   <- est / se
          out[[paste0("beta_qr_PRS_", lab)]]          <- est
          out[[paste0("se_qr_PRS_",   lab)]]          <- se
          out[[paste0("p_qr_PRS_",    lab)]]          <- 2 * pnorm(-abs(z))
          out[[paste0("ci_qr_PRS_lower_", lab)]]      <- fit$asyCI[v, 1]
          out[[paste0("ci_qr_PRS_upper_", lab)]]      <- fit$asyCI[v, 2]
        } else {
          out[[paste0("beta_qr_PRS_", lab)]]          <- est
          out[[paste0("se_qr_PRS_",   lab)]]          <- NA_real_
          out[[paste0("p_qr_PRS_",    lab)]]          <- NA_real_
          out[[paste0("ci_qr_PRS_lower_", lab)]]      <- NA_real_
          out[[paste0("ci_qr_PRS_upper_", lab)]]      <- NA_real_
        }
      } else {
        out[[paste0("beta_qr_PRS_", lab)]]            <- NA_real_
        out[[paste0("se_qr_PRS_",   lab)]]            <- NA_real_
        out[[paste0("p_qr_PRS_",    lab)]]            <- NA_real_
        out[[paste0("ci_qr_PRS_lower_", lab)]]        <- NA_real_
        out[[paste0("ci_qr_PRS_upper_", lab)]]        <- NA_real_
      }
    }, silent = TRUE)
    if (inherits(res, "try-error")) {
      out[[paste0("beta_qr_PRS_", lab)]]              <- NA_real_
      out[[paste0("se_qr_PRS_",   lab)]]              <- NA_real_
      out[[paste0("p_qr_PRS_",    lab)]]              <- NA_real_
      out[[paste0("ci_qr_PRS_lower_", lab)]]          <- NA_real_
      out[[paste0("ci_qr_PRS_upper_", lab)]]          <- NA_real_
    }
  }
  as.data.table(out)
}

pick_lm_prs <- function(smm, ci_mat = NULL) {
  if ("PRS" %in% rownames(smm)) {
    est <- smm["PRS","Estimate"]
    se  <- smm["PRS","Std. Error"]
    p   <- smm["PRS","Pr(>|t|)"]
    lo  <- if (!is.null(ci_mat) && "PRS" %in% rownames(ci_mat)) ci_mat["PRS",1] else NA_real_
    hi  <- if (!is.null(ci_mat) && "PRS" %in% rownames(ci_mat)) ci_mat["PRS",2] else NA_real_
    c(lm_beta_PRS = est, lm_se_PRS = se, lm_p_PRS = p, lm_ci_lower_PRS = lo, lm_ci_upper_PRS = hi)
  } else {
    c(lm_beta_PRS = NA_real_, lm_se_PRS = NA_real_, lm_p_PRS = NA_real_,
      lm_ci_lower_PRS = NA_real_, lm_ci_upper_PRS = NA_real_)
  }
}

S_names <- "PRS"
theta_raw <- theta_from_data(X, Y, taus, S_names)
Sigma_raw <- boot_cov_theta(X, Y, taus, S_names, B = B)

W_joint_raw <- NA_real_; p_joint_raw <- NA_real_
W_eq_raw    <- NA_real_; p_eq_raw    <- NA_real_

if (!is.null(Sigma_raw) && all(is.finite(Sigma_raw))) {
  W_joint_raw <- wald_quad(theta_raw, Sigma_raw)          # df = K
  p_joint_raw <- pchisq(W_joint_raw, df = d, lower.tail = FALSE)
  R <- build_R(K, p)                                      # df = K-1
  RtSR <- R %*% Sigma_raw %*% t(R)
  if (qr(RtSR)$rank == ncol(RtSR)) {
    W_eq_raw <- wald_quad(as.numeric(R %*% theta_raw), RtSR)
    p_eq_raw <- pchisq(W_eq_raw, df = p*(K-1), lower.tail = FALSE)
  } else {
    W_eq_raw <- NA_real_; p_eq_raw <- NA_real_
  }
}

lm_raw <- lm(fm_full, data = mf_raw)   # SAME rows used for X,Y
summ_raw <- summary(lm_raw)$coefficients
ci_raw   <- tryCatch(confint(lm_raw), error = function(e) NULL)
lm_vals_raw <- pick_lm_prs(summ_raw, ci_raw)

qr_cols_raw <- if (write_qr_per_tau) qr_per_tau_stats_prs(X, Y, taus) else data.table()

ratio_cols_raw <- list(); ne_flags_raw <- numeric(length(taus))
if (write_qr_per_tau) {
  Rboot <- matrix(NA_real_, nrow = B, ncol = K)
  for (b in seq_len(B)) {
    db <- mf_raw[sample.int(nrow(mf_raw), nrow(mf_raw), replace = TRUE), , drop = FALSE]
    Xb <- model.matrix(fm_full, db)[,-1,drop=FALSE]
    Yb <- model.response(db)
    lb <- lm(fm_full, data = db)
    beta_ols_b <- summary(lb)$coefficients["PRS","Estimate"]
    for (i in seq_along(taus)) {
      fitb <- conquer(X = Xb, Y = Yb, tau = taus[i], ci = "none", tol = 1e-6)
      names(fitb$coeff) <- c("(Intercept)", colnames(Xb))
      Rboot[b,i] <- unname(fitb$coeff["PRS"] / beta_ols_b)
    }
  }
  for (i in seq_along(taus)) {
    lab <- labs[i]
    ratio_hat <- as.numeric(qr_cols_raw[[paste0("beta_qr_PRS_", lab)]]) / as.numeric(lm_vals_raw["lm_beta_PRS"])
    ci_lo <- as.numeric(quantile(Rboot[,i], probs = alpha_ne/2, na.rm = TRUE))
    ci_hi <- as.numeric(quantile(Rboot[,i], probs = 1 - alpha_ne/2, na.rm = TRUE))
    flag <- as.numeric((ci_hi < band_lo) | (ci_lo > band_hi))
    ratio_cols_raw[[paste0("ratio_qr_over_ols_", lab)]] <- ratio_hat
    ratio_cols_raw[[paste0("ratio_ci_lower_",    lab)]] <- ci_lo
    ratio_cols_raw[[paste0("ratio_ci_upper_",    lab)]] <- ci_hi
    ratio_cols_raw[[paste0("non_equiv_flag_",    lab)]] <- flag
    ne_flags_raw[i] <- flag
  }
}
non_equiv_count_raw <- sum(ne_flags_raw, na.rm = TRUE)

theta_tr <- Sigma_tr <- NULL
W_joint_tr <- p_joint_tr <- W_eq_tr <- p_eq_tr <- NA_real_
lm_tr_vals <- setNames(rep(NA_real_, length(lm_vals_raw)), names(lm_vals_raw))
qr_cols_tr <- data.table()
ratio_cols_tr <- list(); ne_flags_tr <- numeric(length(taus))

if (is.finite(lambda_used)) {
  Y_base <- Y
  Ytr <- if (abs(lambda_used) < 1e-8) log(Y_base) else (Y_base^lambda_used - 1)/lambda_used
  keep <- is.finite(Ytr)
  X_tr <- X[keep, , drop = FALSE]; Ytr <- Ytr[keep]
  cat("N(transformed after finite-drop) =", length(Ytr), "\n")

  theta_tr <- theta_from_data(X_tr, Ytr, taus, S_names)
  Sigma_tr <- boot_cov_theta(X_tr, Ytr, taus, S_names, B = B)
  if (!is.null(Sigma_tr) && all(is.finite(Sigma_tr))) {
    W_joint_tr <- wald_quad(theta_tr, Sigma_tr)
    p_joint_tr <- pchisq(W_joint_tr, df = d, lower.tail = FALSE)
    R <- build_R(K, p)
    RtSR <- R %*% Sigma_tr %*% t(R)
    if (qr(RtSR)$rank == ncol(RtSR)) {
      W_eq_tr <- wald_quad(as.numeric(R %*% theta_tr), RtSR)
      p_eq_tr <- pchisq(W_eq_tr, df = p*(K-1), lower.tail = FALSE)
    } else {
      W_eq_tr <- NA_real_; p_eq_tr <- NA_real_
    }
  }

  mf_dt_tr <- as.data.table(mf_raw)[keep]
  mf_dt_tr[, Y_tr := Ytr]
  lm_tr  <- lm(update(fm_full, Y_tr ~ .), data = mf_dt_tr)
  summ_tr <- summary(lm_tr)$coefficients
  ci_tr   <- tryCatch(confint(lm_tr), error = function(e) NULL)
  lm_tr_vals <- pick_lm_prs(summ_tr, ci_tr)

  if (write_qr_per_tau) qr_cols_tr <- qr_per_tau_stats_prs(X_tr, Ytr, taus)

  if (write_qr_per_tau) {
    Rboot_tr <- matrix(NA_real_, nrow = B, ncol = K)
    for (b in seq_len(B)) {
      db <- mf_dt_tr[sample.int(nrow(mf_dt_tr), nrow(mf_dt_tr), replace = TRUE)]
      Xb <- model.matrix(update(fm_full, Y_tr ~ .), db)[,-1,drop=FALSE]
      Yb <- db$Y_tr
      lb <- lm(update(fm_full, Y_tr ~ .), data = db)
      beta_ols_b <- summary(lb)$coefficients["PRS","Estimate"]
      for (i in seq_along(taus)) {
        fitb <- conquer(X = Xb, Y = Yb, tau = taus[i], ci = "none", tol = 1e-6)
        names(fitb$coeff) <- c("(Intercept)", colnames(Xb))
        Rboot_tr[b,i] <- unname(fitb$coeff["PRS"] / beta_ols_b)
      }
    }
    for (i in seq_along(taus)) {
      lab <- labs[i]
      ratio_hat <- as.numeric(qr_cols_tr[[paste0("beta_qr_PRS_", lab)]]) / as.numeric(lm_tr_vals["lm_beta_PRS"])
      ci_lo <- as.numeric(quantile(Rboot_tr[,i], probs = alpha_ne/2, na.rm = TRUE))
      ci_hi <- as.numeric(quantile(Rboot_tr[,i], probs = 1 - alpha_ne/2, na.rm = TRUE))
      flag <- as.numeric((ci_hi < band_lo) | (ci_lo > band_hi))
      ratio_cols_tr[[paste0("ratio_qr_over_ols_", lab)]] <- ratio_hat
      ratio_cols_tr[[paste0("ratio_ci_lower_",    lab)]] <- ci_lo
      ratio_cols_tr[[paste0("ratio_ci_upper_",    lab)]] <- ci_hi
      ratio_cols_tr[[paste0("non_equiv_flag_",    lab)]] <- flag
      ne_flags_tr[i] <- flag
    }
  }
}
non_equiv_count_tr <- sum(ne_flags_tr, na.rm = TRUE)

row_raw <- data.table(
  phenotype = pheno_id, prs_id = prs_id, scale = "raw",
  lambda_used = NA_real_,
  W_joint = W_joint_raw, p_joint = p_joint_raw,
  W_equality = W_eq_raw, p_equality = p_eq_raw,
  lm_beta_PRS = lm_vals_raw["lm_beta_PRS"],
  lm_se_PRS   = lm_vals_raw["lm_se_PRS"],
  lm_p_PRS    = lm_vals_raw["lm_p_PRS"],
  lm_ci_lower_PRS = lm_vals_raw["lm_ci_lower_PRS"],
  lm_ci_upper_PRS = lm_vals_raw["lm_ci_upper_PRS"],
  non_equiv_count = non_equiv_count_raw
)
row_raw <- cbind(row_raw, qr_cols_raw, as.data.table(ratio_cols_raw))

row_tr <- data.table(
  phenotype = pheno_id, prs_id = prs_id, scale = "transformed",
  lambda_used = lambda_used,
  W_joint = W_joint_tr, p_joint = p_joint_tr,
  W_equality = W_eq_tr, p_equality = p_eq_tr,
  lm_beta_PRS = lm_tr_vals["lm_beta_PRS"],
  lm_se_PRS   = lm_tr_vals["lm_se_PRS"],
  lm_p_PRS    = lm_tr_vals["lm_p_PRS"],
  lm_ci_lower_PRS = lm_tr_vals["lm_ci_lower_PRS"],
  lm_ci_upper_PRS = lm_tr_vals["lm_ci_upper_PRS"],
  non_equiv_count = non_equiv_count_tr
)
row_tr <- cbind(row_tr, qr_cols_tr, as.data.table(ratio_cols_tr))

res <- rbind(row_raw, row_tr, fill = TRUE)

add_header <- !file.exists(out_txt)
fwrite(res, out_txt, sep = "\t", append = !add_header, col.names = add_header)
cat("Appended rows to:\n  ", out_txt, "\n")
