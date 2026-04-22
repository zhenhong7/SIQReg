#!/usr/bin/env Rscript
# QR heterogeneity tests for a single covariate (raw vs Box-Cox)
# Usage: Rscript qr_cov_tests_check.R <PHENO_ID> <COVAR_NAME> <MODE>

suppressPackageStartupMessages({
  library(conquer)
  library(data.table)
  library(MASS)      # ginv
  library(Matrix)    # bdiag
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript qr_cov_tests_check.R <PHENO_ID> <COVAR_NAME> <MODE{mean|median}>")
}
pheno_id     <- args[1]
covar_target <- args[2]            # one of: sex, age, age2, PC1..PC40
mode_lambda  <- match.arg(args[3], c("mean","median"))

pheno_dir    <- Sys.getenv("PHENO_DIR",       "phenotypes")
covar_file   <- Sys.getenv("COVAR_FILE",      "phenotypes/covariates_sa40PC.pheno")
prs_dir      <- Sys.getenv("PRS_DIR",         "PRS/ALL_PRS")
top1_map_csv <- Sys.getenv("TOP1_PRS_CSV",    "PRS/R2/top1_prs_by_pheno.csv")

meta_both_file <- Sys.getenv("META_BOTH_FILE", "Lambda_learning/meta_both_summary_all.txt")

keep_file <- Sys.getenv("KEEP_FILE", "ancestry_ids/whitebrit.fid_iid.txt")
if (!file.exists(keep_file)) stop("keep_file not found: ", keep_file)
keep_iid <- as.integer(fread(keep_file, header = FALSE)[[2]])

out_base <- file.path(Sys.getenv("OUT_BASE", "results/qr_cov_tests"),
                      if (mode_lambda == "mean") "mean_meta" else "median_meta")
out_dir <- file.path(out_base, pheno_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_txt <- file.path(out_dir, paste0(pheno_id, "_qr_cov_tests.txt"))

taus <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
K    <- length(taus)
S_covars <- covar_target
p    <- length(S_covars)          # = 1
d    <- p * K                     # = K
B    <- 100
set.seed(2025)
write_qr_per_tau <- TRUE

# Age handling toggle: TRUE if 'age' in covariates is actually Year-of-Birth
age_is_yob <- FALSE
BIRTH_YEAR_BASELINE <- as.integer(Sys.getenv("BIRTH_YEAR_BASELINE", "2025"))

tau_label <- function(t) paste0("t", sprintf("%03d", as.integer(round(1000*t))))

pheno_path <- list.files(pheno_dir, pattern = paste0("^", pheno_id, ".*\\.pheno$"),
                         full.names = TRUE)
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
base_data <- base_data[IID %in% keep_iid]

top1 <- fread(top1_map_csv)
if (!all(c("phenotype","top_prs") %in% names(top1))) {
  stop("top1_prs_by_pheno.csv must have columns: phenotype, top_prs")
}
prs_id <- top1[phenotype == pheno_id, top_prs][1]
if (is.na(prs_id)) stop("No top PRS mapping for phenotype: ", pheno_id)

prs_file <- file.path(prs_dir, paste0(prs_id, ".pheno"))
if (!file.exists(prs_file)) stop("PRS file not found: ", prs_file)
prs <- fread(prs_file)
setnames(prs, names(prs)[3], "PRS")
prs[, IID := as.integer(IID)]

base_data <- merge(prs, base_data, by = "IID")
rm(prs); invisible(gc())

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
if (nrow(X) != length(Y)) stop("X/Y mismatch after model.frame: ",
                               nrow(X), " vs ", length(Y))

if (!covar_target %in% colnames(X)) stop("Target covariate not in design matrix: ", covar_target)

cat("Target covariate:", covar_target, "\n")
cat("N(raw after NA drop) =", nrow(X), "\n")

is_const <- vapply(seq_len(ncol(X)), function(j) var(X[, j]) == 0, logical(1))
if (any(is_const)) {
  cat("Dropping constant columns in X:", paste(colnames(X)[is_const], collapse = ", "), "\n")
  X <- X[, !is_const, drop = FALSE]
}

qrX <- qr(X)
if (qrX$rank < ncol(X)) {
  keep_idx <- sort(qrX$pivot[seq_len(qrX$rank)])
  dropped_names <- colnames(X)[setdiff(seq_len(ncol(X)), keep_idx)]
  cat("Dropping aliased columns in X:", paste(dropped_names, collapse = ", "), "\n")
  if (covar_target %in% dropped_names) {
    cat("Target", covar_target, "is aliased/constant; writing NA rows and exiting.\n")
    row_na_raw <- data.table(
      phenotype = pheno_id, prs_id = prs_id, target = covar_target, scale = "raw",
      lambda_used = NA_real_,
      W_joint = NA_real_, p_joint = NA_real_,
      W_equality = NA_real_, p_equality = NA_real_,
      lm_beta_PRS = NA_real_, lm_p_PRS = NA_real_,
      lm_beta_target = NA_real_, lm_p_target = NA_real_
    )
    row_na_tr <- copy(row_na_raw); row_na_tr[, scale := "transformed"]
    res_na <- rbind(row_na_raw, row_na_tr, fill = TRUE)
    add_header <- !file.exists(out_txt)
    fwrite(res_na, out_txt, sep = "\t", append = !add_header, col.names = add_header)
    quit(save = "no", status = 0)
  }
  X <- X[, keep_idx, drop = FALSE]
}

lambda_used <- NA_real_
meta_tab <- fread(meta_both_file)
if (!all(c("phenotype","lambda_DL","lambda_WM") %in% names(meta_tab))) {
  stop("meta_both_summary_all.txt must have columns: phenotype, lambda_DL, lambda_WM")
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

qr_per_tau_stats <- function(X, Y, taus, target_var) {
  out <- list()
  for (tau in taus) {
    lab <- tau_label(tau)
    res <- try({
      fit <- conquer(X = X, Y = Y, tau = tau, ci = "asymptotic", tol = 1e-6)
      names(fit$coeff) <- c("(Intercept)", colnames(X))
      if (!is.null(fit$asyCI) && length(fit$asyCI) > 0L) {
        rownames(fit$asyCI) <- c("(Intercept)", colnames(X))
      }

      for (v in c("PRS", target_var)) {
        nm <- if (v == "PRS") "PRS" else "target"

        if (!v %in% names(fit$coeff)) {
          out[[paste0("beta_qr_", nm, "_", lab)]]          <- NA_real_
          out[[paste0("se_qr_",   nm, "_", lab)]]          <- NA_real_
          out[[paste0("p_qr_",    nm, "_", lab)]]          <- NA_real_
          out[[paste0("ci_qr_",   nm, "_lower_", lab)]]    <- NA_real_
          out[[paste0("ci_qr_",   nm, "_upper_", lab)]]    <- NA_real_
          next
        }

        est <- fit$coeff[v]

        if (!is.null(fit$asyCI) && v %in% rownames(fit$asyCI)) {
          lo <- fit$asyCI[v, 1]
          hi <- fit$asyCI[v, 2]
          se <- (hi - est) / qnorm(0.975)
          z  <- est / se

          out[[paste0("beta_qr_", nm, "_", lab)]]          <- est
          out[[paste0("se_qr_",   nm, "_", lab)]]          <- se
          out[[paste0("p_qr_",    nm, "_", lab)]]          <- 2 * pnorm(-abs(z))
          out[[paste0("ci_qr_",   nm, "_lower_", lab)]]    <- lo
          out[[paste0("ci_qr_",   nm, "_upper_", lab)]]    <- hi
        } else {
          out[[paste0("beta_qr_", nm, "_", lab)]]          <- est
          out[[paste0("se_qr_",   nm, "_", lab)]]          <- NA_real_
          out[[paste0("p_qr_",    nm, "_", lab)]]          <- NA_real_
          out[[paste0("ci_qr_",   nm, "_lower_", lab)]]    <- NA_real_
          out[[paste0("ci_qr_",   nm, "_upper_", lab)]]    <- NA_real_
        }
      }
    }, silent = TRUE)

    if (inherits(res, "try-error")) {
      for (nm in c("PRS", "target")) {
        out[[paste0("beta_qr_", nm, "_", lab)]]          <- NA_real_
        out[[paste0("se_qr_",   nm, "_", lab)]]          <- NA_real_
        out[[paste0("p_qr_",    nm, "_", lab)]]          <- NA_real_
        out[[paste0("ci_qr_",   nm, "_lower_", lab)]]    <- NA_real_
        out[[paste0("ci_qr_",   nm, "_upper_", lab)]]    <- NA_real_
      }
    }
  }
  as.data.table(out)
}

pick_lm <- function(smm, nm) {
  if (nm %in% rownames(smm)) smm[nm, c("Estimate","Pr(>|t|)")] else c(NA, NA)
}

theta_raw <- theta_from_data(X, Y, taus, S_covars)
Sigma_raw <- boot_cov_theta(X, Y, taus, S_covars, B = B)

W_joint_raw <- NA_real_; p_joint_raw <- NA_real_
W_eq_raw    <- NA_real_; p_eq_raw    <- NA_real_

if (!is.null(Sigma_raw) && all(is.finite(Sigma_raw))) {
  W_joint_raw <- wald_quad(theta_raw, Sigma_raw)
  p_joint_raw <- pchisq(W_joint_raw, df = d, lower.tail = FALSE)

  R <- build_R(K, p)
  RtSR <- R %*% Sigma_raw %*% t(R)
  if (qr(RtSR)$rank == ncol(RtSR)) {
    W_eq_raw <- wald_quad(as.numeric(R %*% theta_raw), RtSR)
    p_eq_raw <- pchisq(W_eq_raw, df = p*(K-1), lower.tail = FALSE)
  }
}

lm_raw <- lm(fm_full, data = mf_raw)
summ_raw <- summary(lm_raw)$coefficients
lm_vals_raw <- c(
  pick_lm(summ_raw, "PRS"),
  pick_lm(summ_raw, covar_target)
)
names(lm_vals_raw) <- c("lm_beta_PRS","lm_p_PRS","lm_beta_target","lm_p_target")

qr_cols_raw <- if (write_qr_per_tau) qr_per_tau_stats(X, Y, taus, covar_target) else data.table()

theta_tr <- Sigma_tr <- NULL
W_joint_tr <- p_joint_tr <- W_eq_tr <- p_eq_tr <- NA_real_
lm_tr_vals <- setNames(rep(NA_real_, length(lm_vals_raw)), names(lm_vals_raw))
qr_cols_tr <- data.table()

if (is.finite(lambda_used)) {
  Y_base <- Y
  Ytr <- if (abs(lambda_used) < 1e-8) log(Y_base) else (Y_base^lambda_used - 1)/lambda_used
  keep <- is.finite(Ytr)
  X_tr <- X[keep, , drop = FALSE]; Ytr <- Ytr[keep]
  cat("N(transformed after finite-drop) =", length(Ytr), "\n")

  theta_tr <- theta_from_data(X_tr, Ytr, taus, S_covars)
  Sigma_tr <- boot_cov_theta(X_tr, Ytr, taus, S_covars, B = B)
  if (!is.null(Sigma_tr) && all(is.finite(Sigma_tr))) {
    W_joint_tr <- wald_quad(theta_tr, Sigma_tr)
    p_joint_tr <- pchisq(W_joint_tr, df = d, lower.tail = FALSE)
    R <- build_R(K, p)
    RtSR <- R %*% Sigma_tr %*% t(R)
    if (qr(RtSR)$rank == ncol(RtSR)) {
      W_eq_tr <- wald_quad(as.numeric(R %*% theta_tr), RtSR)
      p_eq_tr <- pchisq(W_eq_tr, df = p*(K-1), lower.tail = FALSE)
    }
  }

  mf_dt_tr <- as.data.table(mf_raw)[keep]
  mf_dt_tr[, Y_tr := Ytr]
  lm_tr  <- lm(update(fm_full, Y_tr ~ .), data = mf_dt_tr)
  summ_tr <- summary(lm_tr)$coefficients
  lm_tr_vals <- c(
    pick_lm(smm = summ_tr, "PRS"),
    pick_lm(smm = summ_tr, covar_target)
  )
  names(lm_tr_vals) <- c("lm_beta_PRS","lm_p_PRS","lm_beta_target","lm_p_target")

  if (write_qr_per_tau) qr_cols_tr <- qr_per_tau_stats(X_tr, Ytr, taus, covar_target)
}

row_raw <- data.table(
  phenotype = pheno_id, prs_id = prs_id, target = covar_target, scale = "raw",
  lambda_used = NA_real_,
  W_joint = W_joint_raw, p_joint = p_joint_raw,
  W_equality = W_eq_raw, p_equality = p_eq_raw
)
for (nm in names(lm_vals_raw)) row_raw[[nm]] <- lm_vals_raw[[nm]]
row_raw <- cbind(row_raw, qr_cols_raw)

row_tr <- data.table(
  phenotype = pheno_id, prs_id = prs_id, target = covar_target, scale = paste0("transformed_", mode_lambda),
  lambda_used = lambda_used,
  W_joint = W_joint_tr, p_joint = p_joint_tr,
  W_equality = W_eq_tr, p_equality = p_eq_tr
)
for (nm in names(lm_tr_vals)) row_tr[[nm]] <- lm_tr_vals[[nm]]
row_tr <- cbind(row_tr, qr_cols_tr)

res <- rbind(row_raw, row_tr, fill = TRUE)

add_header <- !file.exists(out_txt)
fwrite(res, out_txt, sep = "\t", append = !add_header, col.names = add_header)
cat("Appended rows to:\n  ", out_txt, "\n", sep = "")
