#!/usr/bin/env Rscript
# PRS coefs with subsample-based lambda-learning and timing
# Usage: Rscript lam_learner.R <PHENO_ID> <PRS_ID> <chunk_i>

suppressPackageStartupMessages({
  library(conquer)
  library(data.table)
  library(dplyr)
  library(MASS)
})

BIRTH_YEAR_BASELINE <- as.integer(Sys.getenv("BIRTH_YEAR_BASELINE", "2025"))
n_chunks <- as.integer(Sys.getenv("N_CHUNKS", "25"))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript lam_learner.R <PHENO_ID> <PRS_ID> <chunk_i>")
pheno_id <- args[[1]]
prs_id   <- args[[2]]
chunk_i  <- as.integer(args[3])

keep_file <- Sys.getenv("KEEP_FILE", "ancestry_ids/whitebrit.fid_iid.txt")
if (!file.exists(keep_file)) stop("keep_file not found: ", keep_file)
keep_iid <- as.integer(fread(keep_file, header = FALSE)[[2]])

# 1) Load phenotype and covariates
pheno_dir <- Sys.getenv("PHENO_DIR", "phenotypes")
pheno_path <- list.files(pheno_dir,
                         pattern = paste0("^", pheno_id, ".*\\.pheno$"),
                         full.names = TRUE)
if (length(pheno_path) != 1) stop("Found ", length(pheno_path), " phenos for ", pheno_id)
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
base_data[, IID := as.integer(IID)]

base_data <- base_data[IID %in% keep_iid]

# 2) Load PRS
prs_dir <- Sys.getenv("PRS_DIR", "PRS/ALL_PRS")
prs_file <- file.path(prs_dir, paste0(prs_id, ".pheno"))
prs <- fread(prs_file)
setnames(prs, names(prs)[3], "PRS")
prs[, IID := as.integer(IID)]

prs <- prs[IID %in% keep_iid]

base_data <- merge(prs, base_data, by = "IID")
rm(prs)

# 3) Prepare data
base_data[, age := BIRTH_YEAR_BASELINE - age]
dat <- na.omit(base_data)  # data.table with columns IID, PRS, Pheno, sex, age, PC1..PC40
fm_full <- as.formula(
  paste0("Pheno ~ PRS + sex + age + ", paste0("PC", 1:40, collapse = " + "))
)

# Construct X_qr, y, G, C for full data (used later)
X_qr_full <- model.matrix(fm_full, dat)[,-1, drop = FALSE]
y_full <- dat$Pheno
G_full <- X_qr_full[, 1]
C_full <- X_qr_full[, -1, drop = FALSE]
# Full tau grid
taus <- c(seq(0.1,0.9,by=0.1))

# 4) Define functions
boxcox_transform <- function(y, lambda) {
  if (abs(lambda) < 1e-8) log(y) else (y^lambda - 1)/lambda
}

slopes_at_lambda <- function(lambda, y, G, C, taus) {
  Ystar <- boxcox_transform(y, lambda)
  if (any(is.nan(Ystar)) || any(is.infinite(Ystar))) return(rep(NA_real_, length(taus)))
  # robust center/scale (affine-invariant for r, but numerically stabilizes)
  m <- suppressWarnings(median(Ystar, na.rm = TRUE))
  s <- suppressWarnings(mad(Ystar, constant = 1.4826, na.rm = TRUE))
  if (!is.finite(s) || s < 1e-8) s <- sd(Ystar)
  if (!is.finite(s) || s < 1e-8) s <- 1.0
  Ystar <- (Ystar - m) / s
  sapply(taus, function(tau) {
    fit <- conquer(X = cbind(G, C), Y = Ystar, tau = tau, ci = "none", tol = 1e-6)
    fit[["coeff"]][2]
  })
}

huber_k <- 1.345
eps     <- 1e-8
D_ratio <- function(lambda, y, G, C, taus, B) {
  n <- length(y); K <- length(taus)
  beta0 <- slopes_at_lambda(lambda, y, G, C, taus)
  if (any(is.na(beta0))) return(1e8)
  # Huber M-center instead of mean
  beta_bar <- MASS::huber(beta0, k = huber_k)$mu
  if (!is.finite(beta_bar) || abs(beta_bar) < eps) return(1e8)
  r0 <- 1 - beta0 / beta_bar
  boot_r <- replicate(B, {
    idx  <- sample.int(n, n, replace = TRUE)
    bvec <- slopes_at_lambda(lambda, y[idx], G[idx], C[idx,], taus)
    if (any(is.na(bvec))) return(rep(NA_real_, K))
    bbar <- MASS::huber(bvec, k = huber_k)$mu
    if (!is.finite(bbar) || abs(bbar) < eps) return(rep(NA_real_, K))
    1 - bvec/bbar
  })
  # ensure matrix & finite
  Rmat <- t(boot_r)
  if (!is.matrix(Rmat) || all(!is.finite(Rmat))) return(1e8)
  Sigma_r <- cov(Rmat, use = "complete.obs")
  invR <- tryCatch(solve(Sigma_r), error = function(e) NULL)
  if (is.null(invR)) invR <- MASS::ginv(Sigma_r)
  as.numeric(t(r0) %*% invR %*% r0)
}

estimate_with_ratio <- function(y, G, C, taus,
                                lambda_lower, lambda_upper, B) {
  objfun <- function(lam) D_ratio(lam, y, G, C, taus, B = B)
  sol  <- optimize(objfun, lower = lambda_lower, upper = lambda_upper, tol = 1e-4)
  lam <- sol$minimum
  val <- sol$objective

  halfwidth <- 0.5
  polish_iters = 3
  for (k in seq_len(polish_iters)) {
    a <- max(lambda_lower, lam - halfwidth)
    b <- min(lambda_upper, lam + halfwidth)
    sol2 <- optimize(objfun, lower = a, upper = b, tol = 1e-5)
    if (sol2$objective < val) {
      lam <- sol2$minimum
      val <- sol2$objective
    }
    halfwidth <- halfwidth / 2  # tighten window
  }
  list(lambda_hat = lam, obj_value = val)
}

# 5) Subsample lambda-learning with timing
N_full <- nrow(dat); n_sub <- ceiling(N_full / n_chunks)
B_for_est <- 500
set.seed(2025)
perm <- sample.int(N_full)
chunk_id <- rep(1:n_chunks, each = n_sub, length.out = N_full)
idx_sub  <- perm[ chunk_id == chunk_i ]
message(sprintf("Chunk %d: %d samples", chunk_i, length(idx_sub)))

dat_sub <- dat[idx_sub]
X_qr_sub  <- model.matrix(fm_full, dat_sub)[,-1, drop = FALSE]
y_sub <- dat_sub$Pheno
G_sub <- X_qr_sub[, 1]
C_sub <- X_qr_sub[,-1, drop = FALSE]

t0 <- Sys.time()
est_sub <- estimate_with_ratio(y_sub, G_sub, C_sub, taus,
                               lambda_lower = -5, lambda_upper = 5, B = B_for_est)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat(sprintf("%s,%d,%.6f,%.3f\n", prs_id, chunk_i, est_sub$lambda_hat, elapsed))
