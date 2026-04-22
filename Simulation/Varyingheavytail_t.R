#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(conquer))
suppressPackageStartupMessages(library(MASS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript sim_lambda_vary_tdf_job.R <t_df> [reps=1000] [B=300] [out_dir=./results_tdf] [lambda0=0.5]")
}
t_df   <- as.numeric(args[1])
reps   <- if (length(args) >= 2) as.integer(args[2]) else 1000L
B      <- if (length(args) >= 3) as.integer(args[3]) else 300L
out_dir <- if (length(args) >= 4) args[4] else "./results_tdf"
lambda0 <- if (length(args) >= 5) as.numeric(args[5]) else 0.5

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Fixed knobs
n      <- 5000L
bG     <- 0.1
pW     <- 9L
taus   <- seq(0.1, 0.9, by = 0.1)
K      <- length(taus)
lam_lo <- -5
lam_hi <-  5

## heavy-tail scale
t_scale <- 0.06

## small domain guard for inverse Box-Cox
max_tries <- 200L

## Fixed regressors & coefficients per (t_df, n, lambda0)
seed_X <- as.integer(abs(1e6 * (t_df + 10) + 1e5 * (lambda0 + 10)) + n) %% .Machine$integer.max
set.seed(seed_X)

PRS <- runif(n)
tgt <- runif(n)
W   <- replicate(pW, runif(n))

b0 <- 0.2
aT <- 0.50
aW <- rep(0.003, pW)

G <- tgt
C <- cbind(PRS, W)
X <- cbind(G, C)        # conquer adds intercept internally

## Precompute things needed for ML estimator (OLS on latent scale)
X_full   <- cbind(1, X)
XtX      <- crossprod(X_full)
XtX_inv  <- tryCatch(solve(XtX), error = function(e) MASS::ginv(XtX))

## Fix D: numerically stable Box-Cox transform
boxcox_stable <- function(y, lam){
  logy <- log(y)
  if (any(!is.finite(logy))) return(rep(NA_real_, length(y)))

  if (abs(lam) < 1e-8) return(logy)

  t <- lam * logy
  t <- pmin(pmax(t, -700), 700)
  ylam <- exp(t)
  (ylam - 1) / lam
}

## (1) D(lambda) objective
objfun_factory <- function(y){
  function(lam){
    ## Box-Cox transform (stable)
    Yst <- boxcox_stable(y, lam)
    if (any(!is.finite(Yst))) return(1e8)

    ## standardize once
    m <- suppressWarnings(median(Yst))
    s <- suppressWarnings(mad(Yst, constant = 1.4826))
    if (!is.finite(s) || s < 1e-8) s <- sd(Yst)
    if (!is.finite(s) || s < 1e-8) s <- 1
    Yst <- (Yst - m)/s

    ## beta(tau) on transformed scale
    beta0 <- numeric(K)
    for (k in seq_len(K)) {
      fitk <- conquer(X = X, Y = Yst, tau = taus[k],
                      ci = "none", tol = 1e-6)
      beta0[k] <- fitk$coeff[2]   # coefficient of G
    }
    if (any(!is.finite(beta0))) return(1e8)

    beta_bar <- mean(beta0)
    if (!is.finite(beta_bar) || abs(beta_bar) < 1e-12) return(1e8)
    r0 <- 1 - beta0 / beta_bar

    ## bootstrap covariance of r(tau)
    set.seed(2025)
    bootR <- matrix(NA_real_, nrow = B, ncol = K)

    for (b in 1:B) {
      idx <- sample.int(n, n, TRUE)
      Yb  <- Yst[idx]
      Xb  <- X[idx, , drop = FALSE]

      bb  <- numeric(K); ok <- TRUE
      for (k in seq_len(K)) {
        fitkb <- try(
          conquer(X = Xb, Y = Yb, tau = taus[k],
                  ci = "none", tol = 1e-6),
          silent = TRUE
        )
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

## (2) Box-Cox profile log-likelihood nll
nll_factory <- function(y){
  ylog <- log(y)
  n    <- length(y)

  function(lam){
    ## Box-Cox transform (stable)
    Yst <- boxcox_stable(y, lam)
    if (any(!is.finite(Yst))) return(1e8)

    beta_hat <- XtX_inv %*% crossprod(X_full, Yst)
    mu_hat   <- as.vector(X_full %*% beta_hat)
    res      <- Yst - mu_hat
    sse      <- sum(res^2)
    if (!is.finite(sse) || sse <= 0) return(1e8)

    (n / 2) * log(sse / n) - (lam - 1) * sum(ylog)
  }
}

## Monte Carlo loop
lambda_hat_D  <- rep(NA_real_, reps)
lambda_hat_ML <- rep(NA_real_, reps)

for (r in 1:reps) {
  set.seed(seed_X + 100000L + r)

  ## latent linear model with heavy-tailed t-errors
  tries <- 0L
  repeat {
    tries <- tries + 1L
    eps <- rt(n, df = t_df) * t_scale
    y_star <- b0 + bG * PRS + aT * tgt + as.vector(W %*% aW) + eps

    ## ensure inverse Box-Cox is defined (no global affine map)
    if (abs(lambda0) < 1e-8) break
    if (all(1 + lambda0 * y_star > 0)) break
    if (tries >= max_tries) { y_star <- NULL; break }
  }
  if (is.null(y_star)) next

  ## generate observed y using the true lambda0
  y <- if (abs(lambda0) < 1e-8) exp(y_star) else (1 + lambda0 * y_star)^(1 / lambda0)
  if (any(!is.finite(y)) || any(y <= 0)) next

  obj_D <- objfun_factory(y)
  sol_D <- optimize(obj_D, lower = lam_lo, upper = lam_hi, tol = 1e-4)
  lambda_hat_D[r] <- sol_D$minimum

  nll    <- nll_factory(y)
  sol_ML <- optimize(nll, lower = lam_lo, upper = lam_hi, tol = 1e-4)
  lambda_hat_ML[r] <- sol_ML$minimum
}

## Summaries & Save
mean_hat_D  <- mean(lambda_hat_D,  na.rm = TRUE)
se_hat_D    <- sd(lambda_hat_D,    na.rm = TRUE)
mean_hat_ML <- mean(lambda_hat_ML, na.rm = TRUE)
se_hat_ML   <- sd(lambda_hat_ML,   na.rm = TRUE)
kept        <- sum(is.finite(lambda_hat_D) & is.finite(lambda_hat_ML))

cat("lambda0:", lambda0, "\n")
cat("t_df   :", t_df, "\n")
cat("n      :", n, "\n")
cat("reps   :", reps, " kept:", kept, "\n\n")

## per-job output file
df_tag  <- gsub("\\.", "p", format(t_df, nsmall = 1, trim = TRUE))
lam_tag <- gsub("\\.", "p", format(lambda0, trim = TRUE))

outfile <- file.path(out_dir,
                     sprintf("sim_lambda_tdf_df%s_n%d_lam%s.csv",
                             df_tag, n, lam_tag))

df_full <- data.frame(
  n             = n,
  lambda0       = lambda0,
  t_df          = t_df,
  rep           = seq_len(reps),
  lambda_hat_D  = lambda_hat_D,
  lambda_hat_ML = lambda_hat_ML,
  mean_hat_D    = mean_hat_D,
  se_hat_D      = se_hat_D,
  mean_hat_ML   = mean_hat_ML,
  se_hat_ML     = se_hat_ML,
  kept          = kept
)
write.csv(df_full, outfile, row.names = FALSE)
cat(sprintf("Wrote %s\n", outfile))

## appended summary
summary_file <- file.path(out_dir, "sim_lambda_tdf_summary.csv")
df_sum <- data.frame(
  n           = n,
  lambda0     = lambda0,
  t_df        = t_df,
  reps        = reps,
  kept        = kept,
  mean_hat_D  = mean_hat_D,
  se_hat_D    = se_hat_D,
  mean_hat_ML = mean_hat_ML,
  se_hat_ML   = se_hat_ML
)

if (!file.exists(summary_file)) {
  write.table(df_sum, summary_file,
              sep = ",", row.names = FALSE, col.names = TRUE)
} else {
  write.table(df_sum, summary_file,
              sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
}
