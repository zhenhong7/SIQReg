#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(conquer))
suppressPackageStartupMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript sim_lambda_one_rep.R <lambda0> <n> <rep_id> [B=300]")
}
lambda0 <- as.numeric(args[1])
n       <- as.integer(args[2])
rep_id  <- as.integer(args[3])
B       <- if (length(args) >= 4) as.integer(args[4]) else 300L

## Fixed knobs
pW     <- 9L
taus   <- seq(0.1, 0.9, by = 0.1)
K      <- length(taus)
lam_lo <- -5
lam_hi <-  5

t_df    <- 1.5
t_scale <- 0.06

max_tries     <- 200L
domain_margin <- 1e-6

## Fixed design per (lambda0, n)
seed_X <- as.integer(abs(1e6 * (lambda0 + 10)) + n) %% .Machine$integer.max
set.seed(seed_X)

PRS <- runif(n)
tgt <- runif(n)
W   <- replicate(pW, runif(n))

b0 <- -0.1
bG <- 0.1
aT <- 0.50
aW <- rep(0.003, pW)

G <- tgt
C <- cbind(PRS, W)
X <- cbind(G, C)

X_full   <- cbind(1, X)
XtX      <- crossprod(X_full)
XtX_inv  <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))

boxcox_stable <- function(y, lam){
  logy <- log(y)
  if (any(!is.finite(logy))) return(rep(NA_real_, length(y)))
  if (abs(lam) < 1e-8) return(logy)
  t <- lam * logy
  t <- pmin(pmax(t, -700), 700)
  ylam <- exp(t)
  (ylam - 1) / lam
}

objfun_factory <- function(y){
  function(lam){
    Yst <- boxcox_stable(y, lam)
    if (any(!is.finite(Yst))) return(1e8)

    m <- suppressWarnings(median(Yst))
    s <- suppressWarnings(mad(Yst, constant = 1.4826))
    if (!is.finite(s) || s < 1e-8) s <- sd(Yst)
    if (!is.finite(s) || s < 1e-8) s <- 1
    Yst <- (Yst - m)/s

    beta0 <- numeric(K)
    for (k in seq_len(K)) {
      fitk <- conquer(X = X, Y = Yst, tau = taus[k], ci = "none", tol = 1e-6)
      beta0[k] <- fitk$coeff[2]
    }
    if (any(!is.finite(beta0))) return(1e8)

    beta_bar <- mean(beta0)
    if (!is.finite(beta_bar) || abs(beta_bar) < 1e-12) return(1e8)
    r0 <- 1 - beta0 / beta_bar

    set.seed(2025)
    bootR <- matrix(NA_real_, nrow = B, ncol = K)

    for (b in 1:B) {
      idx <- sample.int(n, n, TRUE)
      Yb  <- Yst[idx]
      Xb  <- X[idx, , drop = FALSE]

      bb <- numeric(K); ok <- TRUE
      for (k in seq_len(K)) {
        fitkb <- try(conquer(X = Xb, Y = Yb, tau = taus[k], ci = "none", tol = 1e-6),
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
    invR <- tryCatch(solve(Sigma_r), error = function(e) MASS::ginv(Sigma_r))

    as.numeric(t(r0) %*% invR %*% r0)
  }
}

nll_factory <- function(y){
  ylog <- log(y); nn <- length(y)
  function(lam){
    Yst <- boxcox_stable(y, lam)
    if (any(!is.finite(Yst))) return(1e8)

    beta_hat <- XtX_inv %*% crossprod(X_full, Yst)
    mu_hat   <- as.vector(X_full %*% beta_hat)
    res      <- Yst - mu_hat
    sse      <- sum(res^2)
    if (!is.finite(sse) || sse <= 0) return(1e8)

    (nn/2) * log(sse/nn) - (lam - 1) * sum(ylog)
  }
}

## One rep: only error changes
set.seed(seed_X + 100000L + rep_id)

mu <- b0 + bG * PRS + aT * tgt + as.vector(W %*% aW)
eps <- rt(n, df = t_df) * t_scale
y_star <- mu + eps

## resample ONLY violators so 1 + lambda0*y_star > 0
if (abs(lambda0) >= 1e-8) {
  bound <- -1 / lambda0
  tries <- 0L
  repeat {
    tries <- tries + 1L
    bad <- if (lambda0 > 0) (y_star <= (bound + domain_margin)) else (y_star >= (bound - domain_margin))
    if (!any(bad)) break
    if (tries >= max_tries) { y_star <- NULL; break }
    eps[bad]    <- rt(sum(bad), df = t_df) * t_scale
    y_star[bad] <- mu[bad] + eps[bad]
  }
}

hatD <- NA_real_
hatM <- NA_real_

if (!is.null(y_star)) {
  y <- if (abs(lambda0) < 1e-8) exp(y_star) else (1 + lambda0 * y_star)^(1 / lambda0)
  if (all(is.finite(y)) && all(y > 0)) {
    hatD <- optimize(objfun_factory(y), lower = lam_lo, upper = lam_hi, tol = 1e-4)$minimum
    hatM <- optimize(nll_factory(y),     lower = lam_lo, upper = lam_hi, tol = 1e-4)$minimum
  }
}

cat(rep_id, hatD, hatM, sep = ",")
cat("\n")
