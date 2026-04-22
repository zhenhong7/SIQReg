#!/usr/bin/env Rscript
# Covariate lambda-learning via QR-ratio for PHENO ~ top1PRS + covariates
# Usage: Rscript lambda_covariate.R <PHENO_ID> <chunk_i> <TARGET_VAR>

suppressPackageStartupMessages({
  library(conquer)
  library(data.table)
  library(dplyr)
  library(MASS)
})

BIRTH_YEAR_BASELINE <- as.integer(Sys.getenv("BIRTH_YEAR_BASELINE", "2025"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3)
  stop("Usage: Rscript lambda_covariate.R <PHENO_ID> <chunk_i> <TARGET_VAR>")
pheno_id   <- args[1]
chunk_i    <- as.integer(args[2])
n_chunks <- as.integer(Sys.getenv("N_CHUNKS", "25"))
if (!is.finite(n_chunks) || n_chunks < 1) stop("Bad N_CHUNKS: ", Sys.getenv("N_CHUNKS"))
if (chunk_i < 1 || chunk_i > n_chunks) stop("chunk_i must be in 1..", n_chunks)

target_var <- args[3]  # one of: sex, age, age2, PC1..PC40

top1_csv <- Sys.getenv("TOP1_CSV", "PRS/top1_prs_by_pheno.csv")

anc <- trimws(Sys.getenv("ANCESTRY", "whitebrit"))  # whitebrit | white_euro | afr | asn
keep_dir <- Sys.getenv("KEEP_DIR", "ancestry_ids")
keep_file <- file.path(keep_dir, paste0(anc, ".fid_iid.txt"))
if (!file.exists(keep_file)) stop("keep_file not found: ", keep_file)
keep_iid <- as.character(fread(keep_file, header = FALSE)[[2]])

.err <- new.env(parent = emptyenv())
.err$n   <- 0L
.err$max <- 50L
log_fail <- function(stage, lambda, tau, msg) {
  if (.err$n < .err$max) {
    message(sprintf("[WARN][%s] pheno=%s target=%s chunk=%d lambda=%.4f tau=%s :: %s",
                    stage, pheno_id, target_var, chunk_i, lambda,
                    format(tau, digits = 3), msg))
    .err$n <- .err$n + 1L
  }
}

# 1) Load phenotype and covariates
pheno_dir <- Sys.getenv("PHENO_DIR", "phenotypes")
pheno_path <- list.files(pheno_dir,
                         pattern = paste0("^", pheno_id, ".*\\.pheno$"),
                         full.names = TRUE)
if (length(pheno_path) != 1)
  stop("Found ", length(pheno_path), " pheno files for ", pheno_id)

pheno <- fread(pheno_path)
setnames(pheno, c("FID","IID","Pheno"))
pheno[, IID := as.character(IID)]

covar_file <- Sys.getenv("COVAR_FILE", "covariates_sa40PC.pheno")
covar <- fread(covar_file)
setnames(covar, names(covar)[3:44], c("sex","age", paste0("PC",1:40)))
covar[, IID := as.character(IID)]

base_data <- merge(pheno, covar, by = "IID")[!is.na(Pheno)]
rm(pheno, covar); invisible(gc())
setDT(base_data)

# 2) Load top-1 PRS for this phenotype
top1 <- fread(top1_csv)   # columns: phenotype, top_prs, r2
if (!"phenotype" %in% names(top1) || !"top_prs" %in% names(top1))
  stop("top1_prs_by_pheno.csv must have columns: phenotype, top_prs")

prs_id <- top1[phenotype == pheno_id, top_prs][1]
if (is.na(prs_id)) stop("No top PRS mapping found for phenotype: ", pheno_id)

prs_dir <- Sys.getenv("PRS_DIR", "PRS/ALL_PRS")
prs_file <- file.path(prs_dir, paste0(prs_id, ".pheno"))
if (!file.exists(prs_file)) stop("PRS file not found: ", prs_file)

prs <- fread(prs_file)
setnames(prs, names(prs)[3], "PRS")
prs[, IID := as.character(IID)]

base_data <- merge(prs, base_data, by = "IID")
rm(prs); invisible(gc())

base_data <- base_data[IID %in% keep_iid]

# 3) Prepare data  (add age2)
base_data[, age := BIRTH_YEAR_BASELINE - age]
mu_age <- mean(base_data$age, na.rm = TRUE)
base_data[, age := age - mu_age]
base_data[, age2 := age^2]

dat <- na.omit(base_data)  # IID, PRS, Pheno, sex, age, age2, PC1..PC40

fm_full <- as.formula(
  paste0("Pheno ~ PRS + sex + age + age2 + ", paste0("PC", 1:40, collapse = " + "))
)

covar_names <- c("sex","age","age2", paste0("PC",1:40))
if (!target_var %in% covar_names)
  stop("TARGET_VAR must be one of: ", paste(covar_names, collapse=", "))

taus <- seq(0.1, 0.9, by = 0.1)

# 4) QR-ratio machinery (only conquer() wrapped)
boxcox_transform <- function(y, lambda) {
  if (abs(lambda) < 1e-8) log(y) else (y^lambda - 1)/lambda
}

slopes_at_lambda <- function(lambda, y, G, C, taus) {
  Ystar <- boxcox_transform(y, lambda)
  if (any(is.nan(Ystar)) || any(is.infinite(Ystar))) {
    log_fail("Ystar_nonfinite", lambda, NA_real_, "Ystar has NaN/Inf")
    return(rep(NA_real_, length(taus)))
  }

  # robust center/scale (stabilizes quantile fits during lambda search)
  m <- suppressWarnings(median(Ystar, na.rm = TRUE))
  s <- suppressWarnings(mad(Ystar, constant = 1.4826, na.rm = TRUE))
  if (!is.finite(s) || s < 1e-8) s <- sd(Ystar)
  if (!is.finite(s) || s < 1e-8) s <- 1.0
  Ystar <- (Ystar - m) / s

  sapply(taus, function(tau) {
    fit <- tryCatch(
      conquer(X = cbind(G, C), Y = Ystar, tau = tau, ci = "none", tol = 1e-6),
      error = function(e) {
        log_fail("conquer", lambda, tau, conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) return(NA_real_)
    b <- tryCatch(
      fit[["coeff"]][2],
      error = function(e) {
        log_fail("coef_extract", lambda, tau, conditionMessage(e))
        NA_real_
      }
    )
    as.numeric(b)
  })
}

D_ratio <- function(lambda, y, G, C, taus, B) {
  n <- length(y); K <- length(taus)
  beta0 <- slopes_at_lambda(lambda, y, G, C, taus)
  if (any(is.na(beta0))) return(1e8)
  beta_bar <- mean(beta0)
  if (beta_bar == 0)  return(1e8)
  r0 <- 1 - beta0 / beta_bar
  boot_r <- replicate(B, {
    idx  <- sample.int(n, n, replace = TRUE)
    bvec <- slopes_at_lambda(lambda, y[idx], G[idx], C[idx,], taus)
    if (any(is.na(bvec))) return(rep(NA_real_, K))
    bbar <- mean(bvec)
    if (bbar == 0) return(rep(NA_real_, K))
    1 - bvec/bbar
  })
  Rmat <- t(boot_r)
  if (!is.matrix(Rmat) || all(!is.finite(Rmat))) return(1e8)
  Sigma_r <- cov(Rmat, use = "complete.obs")
  invR <- tryCatch(solve(Sigma_r), error = function(e) NULL)
  if (is.null(invR)) invR <- MASS::ginv(Sigma_r)
  as.numeric(t(r0) %*% invR %*% r0)
}

estimate_with_ratio <- function(y, G, C, taus, lambda_lower, lambda_upper, B,
                                seed = NULL) {
  objfun <- function(lam) {
    if (!is.null(seed)) set.seed(seed)   # deterministic bootstrap across lam
    D_ratio(lam, y, G, C, taus, B = B)
  }
  sol  <- optimize(objfun, lower = lambda_lower, upper = lambda_upper, tol = 1e-4)
  lam <- sol$minimum; val <- sol$objective
  list(lambda_hat = lam, obj_value = val)
}

# 5) Subsample lambda-learning with timing
N_full <- nrow(dat); n_sub <- ceiling(N_full / n_chunks)
B_for_est <- 500
set.seed(2025)
perm <- sample.int(N_full)
chunk_id <- rep(1:n_chunks, each = n_sub, length.out = N_full)
idx_sub  <- perm[ chunk_id == chunk_i ]
message(sprintf("[%s] Chunk %d: %d samples", pheno_id, chunk_i, length(idx_sub)))

dat_sub  <- dat[idx_sub]
X_sub    <- model.matrix(fm_full, dat_sub)[,-1, drop = FALSE]  # PRS, sex, age, age2, PC1..PC40
y_sub    <- dat_sub$Pheno
j <- match(target_var, colnames(X_sub))
if (is.na(j)) stop("TARGET_VAR not found in design matrix: ", target_var)
G_sub <- X_sub[, j]                 # target covariate
C_sub <- X_sub[, -j, drop = FALSE]  # PRS + the other covariates

# Deterministic seed per (pheno, chunk, target)
seed_base <- (sum(as.integer(charToRaw(pheno_id))) %% 1e6) +
  1009L * chunk_i +
  17L   * match(target_var, covar_names)
seed_base <- as.integer(abs(seed_base) %% .Machine$integer.max)

t0 <- Sys.time()
est_sub <- estimate_with_ratio(y_sub, G_sub, C_sub, taus,
                               lambda_lower = -5, lambda_upper = 5, B = B_for_est,
                               seed = seed_base)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

out_root <- Sys.getenv("OUTROOT", "results_covariates")
out_dir  <- file.path(out_root, anc, pheno_id)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, paste0(pheno_id, "_", target_var, "_lambda.csv"))

fwrite(data.table(
  phenotype  = pheno_id,
  prs_id     = prs_id,
  target     = target_var,
  chunk      = chunk_i,
  lambda_hat = est_sub$lambda_hat,
  obj_value  = est_sub$obj_value,
  seconds    = elapsed,
  n_sub      = length(y_sub)
), out_file, append = file.exists(out_file))
