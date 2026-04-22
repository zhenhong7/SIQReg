#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(conquer))
suppressPackageStartupMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript sim_prs_hetero_gamma_power_job.R <gamma> <rep_id>\n",
       call. = FALSE)
}
gamma  <- as.numeric(args[1])   # 0..2 by 0.2
rep_id <- as.integer(args[2])   # 1..1000

lambda0 <- 0.5
n       <- 5000L
B       <- 300L
pW      <- 9L
taus    <- seq(0.1, 0.9, by = 0.1)
K       <- length(taus)
lam_lo  <- -5
lam_hi  <-  5

b0 <- 0
bG <- 0.4
aT <- 0.5
aW <- rep(0, pW)

## heteroskedastic baseline (same scale as your old eps_sd)
sigma0 <- 0.03

equality_p_qr <- function(Y, G,
                          quant_levels = seq(0.1, 0.9, by = 0.1),
                          n_boot = 300) {
  Y <- as.numeric(Y); G <- as.numeric(G)
  good <- is.finite(Y) & is.finite(G)
  Y <- Y[good]; G <- G[good]
  n <- length(Y)
  if (n < 50L) return(NA_real_)

  K <- length(quant_levels)
  X <- cbind(G = G)

  orig <- sapply(quant_levels, function(tau) {
    fit <- conquer(X = X, Y = Y, tau = tau, ci = "none")
    if (length(fit$coeff) == 0L) return(NA_real_)
    fit$coeff[2]
  })
  if (any(!is.finite(orig))) return(NA_real_)

  boot_coefs <- matrix(NA_real_, nrow = n_boot, ncol = K)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    Xb  <- X[idx, , drop = FALSE]
    Yb  <- Y[idx]
    boot_coefs[b, ] <- sapply(quant_levels, function(tau) {
      fit <- conquer(X = Xb, Y = Yb, tau = tau, ci = "none")
      if (length(fit$coeff) == 0L) return(NA_real_)
      fit$coeff[2]
    })
  }

  cov_mat <- cov(boot_coefs, use = "complete.obs")
  if (any(!is.finite(cov_mat))) return(NA_real_)

  R    <- diag(K)[-K, ] - diag(K)[-1, ]
  RsRt <- R %*% cov_mat %*% t(R)

  inv_RsRt <- tryCatch(solve(RsRt), error = function(e) NULL)
  if (is.null(inv_RsRt)) inv_RsRt <- MASS::ginv(RsRt)

  diff_beta <- R %*% orig
  stat <- as.numeric(t(diff_beta) %*% inv_RsRt %*% diff_beta)
  df   <- K - 1L

  pchisq(stat, df = df, lower.tail = FALSE)
}

rint_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

objfun_factory <- function(y, X){
  function(lam){
    Yst <- if (abs(lam) < 1e-8) log(y) else (y^lam - 1)/lam
    if (any(!is.finite(Yst))) return(1e8)

    m <- suppressWarnings(median(Yst))
    s <- suppressWarnings(mad(Yst, constant = 1.4826))
    if (!is.finite(s) || s < 1e-8) s <- sd(Yst)
    if (!is.finite(s) || s < 1e-8) s <- 1
    Yst <- (Yst - m)/s

    beta0 <- numeric(K)
    for (k in seq_len(K)) {
      fitk <- conquer(X = X, Y = Yst, tau = taus[k], ci = "none", tol = 1e-6)
      beta0[k] <- fitk$coeff[2]  # slope on G=tgt
    }
    if (any(!is.finite(beta0))) return(1e8)

    beta_bar <- mean(beta0)
    if (!is.finite(beta_bar) || abs(beta_bar) < 1e-12) return(1e8)
    r0 <- 1 - beta0 / beta_bar

    set.seed(2025)
    bootR <- matrix(NA_real_, nrow = B, ncol = K)

    for (b in seq_len(B)) {
      idx <- sample.int(length(Yst), length(Yst), TRUE)
      Yb  <- Yst[idx]
      Xb  <- X[idx, , drop = FALSE]

      bb <- numeric(K); ok <- TRUE
      for (k in seq_len(K)) {
        fitkb <- try(conquer(X = Xb, Y = Yb, tau = taus[k],
                             ci = "none", tol = 1e-6),
                     silent = TRUE)
        if (inherits(fitkb, "try-error")) { ok <- FALSE; break }
        bb[k] <- fitkb$coeff[2]
      }
      if (!ok || any(!is.finite(bb))) next

      bbm <- mean(bb)
      if (!is.finite(bbm) || abs(bbm) < 1e-12) next
      bootR[b, ] <- 1 - bb / bbm
    }

    Rmat <- bootR[stats::complete.cases(bootR), , drop = FALSE]
    if (nrow(Rmat) < (K + 5)) return(1e8)

    Sigma_r <- stats::cov(Rmat, use = "complete.obs")
    invR <- tryCatch(solve(Sigma_r), error = function(e) NULL)
    if (is.null(invR)) invR <- MASS::ginv(Sigma_r)

    as.numeric(t(r0) %*% invR %*% r0)
  }
}

base_seed <- 123456L
set.seed(base_seed)
tgt <- runif(n)
W   <- replicate(pW, runif(n))

seed_PRS <- (base_seed + 999999L + rep_id) %% .Machine$integer.max
set.seed(seed_PRS)
PRS <- runif(n)

seed_z <- (base_seed + 2000000L + rep_id) %% .Machine$integer.max
set.seed(seed_z)
z <- rnorm(n, 0, 1)

prs_c <- PRS - mean(PRS)
sigma <- sigma0 * exp(gamma * prs_c)
eps   <- z * sigma

y_star <- b0 + bG * PRS + aT * tgt + as.vector(W %*% aW) + eps

## keep your range-control (same as your reference script)
if (min(y_star) <= 0 || max(y_star) >= 0.5) {
  rng <- range(y_star)
  y_star <- 0.05 + (y_star - rng[1]) * (0.40 / (rng[2] - rng[1] + 1e-12))
}

y <- (1 + lambda0 * y_star)^(1 / lambda0)

log_y  <- log(y)
rint_y <- rint_transform(y)

G <- tgt
C <- cbind(PRS, W)
X <- cbind(G, C)

objfun <- objfun_factory(y, X)
sol    <- optimize(objfun, lower = lam_lo, upper = lam_hi, tol = 1e-4)
lambda_hat <- sol$minimum
D_value    <- sol$objective

eq_p_y_star <- equality_p_qr(y_star, PRS)
eq_p_y      <- equality_p_qr(y, PRS)
eq_p_log_y  <- equality_p_qr(log_y, PRS)
eq_p_rint_y <- equality_p_qr(rint_y, PRS)

y_star_est <- if (abs(lambda_hat) < 1e-8) log(y) else (y^lambda_hat - 1) / lambda_hat
eq_p_y_star_est <- equality_p_qr(y_star_est, PRS)

vals <- c(
  n, lambda0, gamma, rep_id, lambda_hat, D_value,
  eq_p_y, eq_p_log_y, eq_p_rint_y, eq_p_y_star, eq_p_y_star_est
)
cat(paste(format(vals, digits = 16, scientific = TRUE), collapse = ","), "\n")
