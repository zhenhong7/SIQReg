#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(conquer))
suppressPackageStartupMessages(library(MASS))

## Usage: Rscript power_heavytail.R <bG_index 1-11> <rep_index 1-1000>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript power_heavytail.R <bG_index 1-11> <rep_index 1-1000>\n",
       call. = FALSE)
}
idx     <- as.integer(args[1])
rep_idx <- as.integer(args[2])

if (is.na(idx) || idx < 1L || idx > 11L) {
  stop("bG_index must be an integer in 1:11", call. = FALSE)
}
if (is.na(rep_idx) || rep_idx < 1L || rep_idx > 1000L) {
  stop("rep_index must be an integer in 1:1000", call. = FALSE)
}

bG_grid <- seq(0, 0.10, by = 0.01)
bG      <- bG_grid[idx]

setTimeLimit(cpu = 3600, elapsed = 3600, transient = TRUE)
out_file <- "results_heavytail.txt"

lambda0   <- 0.5

n         <- 2000L

B         <- 300L
pW        <- 9L
taus      <- seq(0.1, 0.9, by = 0.1)
K         <- length(taus)
lam_lo    <- -5
lam_hi    <- 5

min_good_boot <- K + 20L

## latent mean: Y* = b0 + bG * PRS + aT * tgt + W %*% aW + eps
b0    <- -1.5
aT    <- 1
aW    <- rep(0.3, pW)    ## nonzero covariate effects

## Heavy-tailed errors
eps_sd <- 0.15
rint_transform <- function(x) {
  r <- rank(x, na.last = "keep")
  m <- sum(!is.na(x))
  qnorm((r - 0.5) / m)
}

qr_joint <- function(Y, G, Covs = NULL,
                     quant_levels = seq(0.1, 0.9, by = 0.1),
                     n_boot = 300) {
  Y <- as.numeric(Y)
  G <- as.numeric(G)
  good <- is.finite(Y) & is.finite(G)
  if (!is.null(Covs)) {
    good <- good & apply(Covs, 1, function(r) all(is.finite(r)))
  }
  Y <- Y[good]; G <- G[good]
  if (!is.null(Covs)) Covs <- Covs[good, , drop = FALSE]
  n <- length(Y)
  if (n < 50L) return(list(joint_p = NA_real_))

  Kq <- length(quant_levels)
  ## G is first column so coeff[2] is the PRS effect
  if (!is.null(Covs)) {
    X <- cbind(G = G, Covs)
  } else {
    X <- cbind(G = G)
  }

  beta_vec <- rep(NA_real_, Kq)
  for (i in seq_along(quant_levels)) {
    tau <- quant_levels[i]
    fit <- try(
      conquer(X = X, Y = Y, tau = tau,
              ci = "asymptotic", tol = 1e-5),
      silent = TRUE
    )
    if (inherits(fit, "try-error") ||
        length(fit$coeff) < 2L ||
        is.null(fit$asyCI) ||
        length(fit$asyCI) == 0L ||
        !is.finite(fit$coeff[2])) {
      beta_vec[i] <- NA_real_
    } else {
      beta_vec[i] <- fit$coeff[2]
    }
  }

  if (any(!is.finite(beta_vec))) {
    return(list(joint_p = NA_real_))
  }

  boot_coefs <- matrix(NA_real_, nrow = n_boot, ncol = Kq)
  for (b in seq_len(n_boot)) {
    idxb <- sample.int(n, n, replace = TRUE)
    Xb   <- X[idxb, , drop = FALSE]
    Yb   <- Y[idxb]

    bb <- rep(NA_real_, Kq)
    ok <- TRUE
    for (i in seq_along(quant_levels)) {
      tau <- quant_levels[i]
      fitb <- try(
        conquer(X = Xb, Y = Yb, tau = tau,
                ci = "none", tol = 1e-5),
        silent = TRUE
      )
      if (inherits(fitb, "try-error") ||
          length(fitb$coeff) < 2L ||
          !is.finite(fitb$coeff[2])) {
        ok <- FALSE
        break
      }
      bb[i] <- fitb$coeff[2]
    }
    if (!ok || any(!is.finite(bb))) next
    boot_coefs[b, ] <- bb
  }

  boot_ok <- boot_coefs[stats::complete.cases(boot_coefs), , drop = FALSE]
  if (nrow(boot_ok) < (Kq + 5L)) {
    return(list(joint_p = NA_real_))
  }

  cov_mat <- stats::cov(boot_ok, use = "complete.obs")
  if (any(!is.finite(cov_mat))) {
    return(list(joint_p = NA_real_))
  }

  inv_cov <- tryCatch(solve(cov_mat),
                      error = function(e) MASS::ginv(cov_mat))

  stat    <- as.numeric(t(beta_vec) %*% inv_cov %*% beta_vec)
  joint_p <- stats::pchisq(stat, df = Kq, lower.tail = FALSE)

  list(joint_p = joint_p)
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
      fitk <- conquer(X = X, Y = Yst, tau = taus[k],
                      ci = "none", tol = 1e-5)
      beta0[k] <- fitk$coeff[2]
    }
    if (any(!is.finite(beta0))) return(1e8)

    beta_bar <- mean(beta0)
    if (!is.finite(beta_bar) || abs(beta_bar) < 1e-12) return(1e8)

    r0 <- 1 - beta0 / beta_bar

    set.seed(2025)
    bootR <- matrix(NA_real_, nrow = B, ncol = K)
    good  <- 0L

    for (b in seq_len(B)) {
      idxb <- sample.int(length(Yst), length(Yst), TRUE)
      Yb   <- Yst[idxb]
      Xb   <- X[idxb, , drop = FALSE]

      bb <- numeric(K)
      ok <- TRUE
      for (k in seq_len(K)) {
        fitkb <- try(
          conquer(X = Xb, Y = Yb, tau = taus[k],
                  ci = "none", tol = 1e-5),
          silent = TRUE
        )
        if (inherits(fitkb, "try-error")) { ok <- FALSE; break }
        bb[k] <- fitkb$coeff[2]
      }
      if (!ok || any(!is.finite(bb))) next

      bbm <- mean(bb)
      if (!is.finite(bbm) || abs(bbm) < 1e-12) next

      bootR[b, ] <- 1 - bb / bbm
      good <- good + 1L

      if (good >= min_good_boot) break
    }

    Rmat <- bootR[stats::complete.cases(bootR), , drop = FALSE]
    if (nrow(Rmat) < (K + 5)) return(1e8)

    Sigma_r <- stats::cov(Rmat, use = "complete.obs")
    if (any(!is.finite(Sigma_r))) return(1e8)

    invR <- tryCatch(solve(Sigma_r),
                     error = function(e) NULL)
    if (is.null(invR)) invR <- MASS::ginv(Sigma_r)

    as.numeric(t(r0) %*% invR %*% r0)
  }
}

seed_design <- 24681357L
set.seed(seed_design)
PRS <- runif(n)
tgt <- runif(n)
W   <- replicate(pW, runif(n))

G     <- PRS
X_tgt <- cbind(tgt, G, W)  # first column = tgt for lambda-learning

## covariate matrix for LR and QR adjustment
Covs <- cbind(tgt, W)

lambda_hat <- NA_real_
D_value    <- NA_real_
p_lr_y          <- NA_real_
p_lr_log        <- NA_real_
p_lr_lat        <- NA_real_
p_lr_rint       <- NA_real_
p_lr_oracle     <- NA_real_
joint_p_y       <- NA_real_
joint_p_log     <- NA_real_
joint_p_latent  <- NA_real_
joint_p_rint    <- NA_real_
joint_p_oracle  <- NA_real_
err_flag  <- FALSE

tryCatch({
  seed_eps <- 97531L + idx * 10000L + rep_idx
  set.seed(seed_eps)
  eps <- rt(n, df = 2.5) * eps_sd

  ## latent scale (purely additive)
  y_star <- b0 + bG * PRS + aT * tgt + as.vector(W %*% aW) + eps

  if (lambda0 < 0) {
    wall   <- 1 - 1e-6
    y_star <- pmin(y_star, wall)
  }

  ## forward Box-Cox to generate observed y
  y <- if (abs(lambda0) < 1e-8) {
    exp(y_star)
  } else {
    (1 + lambda0 * y_star)^(1/lambda0)
  }

  ## derived scales
  log_y  <- log(y)
  rint_y <- rint_transform(y)

  ## learn lambda from X_tgt
  objfun <- objfun_factory(y, X_tgt)
  sol    <- optimize(objfun, lower = lam_lo, upper = lam_hi, tol = 1e-4)

  lambda_hat <<- sol$minimum
  D_value    <<- sol$objective

  ## estimated latent scale
  y_star_est <- if (abs(lambda_hat) < 1e-8) {
    log(y)
  } else {
    (y^lambda_hat - 1) / lambda_hat
  }

  fit_lr_y   <- lm(y ~ PRS + tgt + W)
  p_lr_y     <<- summary(fit_lr_y)$coefficients["PRS","Pr(>|t|)"]

  fit_lr_log <- lm(log_y ~ PRS + tgt + W)
  p_lr_log   <<- summary(fit_lr_log)$coefficients["PRS","Pr(>|t|)"]

  fit_lr_lat <- lm(y_star_est ~ PRS + tgt + W)
  p_lr_lat   <<- summary(fit_lr_lat)$coefficients["PRS","Pr(>|t|)"]

  fit_lr_rint <- lm(rint_y ~ PRS + tgt + W)
  p_lr_rint   <<- summary(fit_lr_rint)$coefficients["PRS","Pr(>|t|)"]

  ## ORACLE
  fit_lr_oracle <- lm(y_star ~ PRS + tgt + W)
  p_lr_oracle   <<- summary(fit_lr_oracle)$coefficients["PRS","Pr(>|t|)"]

  qr_y    <- qr_joint(y,          PRS, Covs)
  qr_log  <- qr_joint(log_y,     PRS, Covs)
  qr_lat  <- qr_joint(y_star_est, PRS, Covs)
  qr_rint <- qr_joint(rint_y,    PRS, Covs)
  qr_oracle <- qr_joint(y_star,  PRS, Covs)

  joint_p_y      <<- qr_y$joint_p
  joint_p_log    <<- qr_log$joint_p
  joint_p_latent <<- qr_lat$joint_p
  joint_p_rint   <<- qr_rint$joint_p
  joint_p_oracle <<- qr_oracle$joint_p
},
error = function(e) {
  message("Error / timeout for bG_index = ", idx,
          " (bG = ", bG, "), rep = ", rep_idx, ": ",
          conditionMessage(e))
  err_flag <<- TRUE
})

res <- data.frame(
  error_design    = "heavy_tail",
  lambda0         = lambda0,
  n               = n,
  b0              = b0,
  bG              = bG,
  aT              = aT,
  rep             = rep_idx,
  lambda_hat      = lambda_hat,
  D_value         = D_value,
  p_lr_y          = p_lr_y,
  p_lr_log        = p_lr_log,
  p_lr_latent     = p_lr_lat,
  p_lr_rint       = p_lr_rint,
  p_lr_oracle     = p_lr_oracle,
  joint_p_y       = joint_p_y,
  joint_p_log     = joint_p_log,
  joint_p_latent  = joint_p_latent,
  joint_p_rint    = joint_p_rint,
  joint_p_oracle  = joint_p_oracle
)

write.table(
  res,
  file      = out_file,
  row.names = FALSE,
  col.names = !file.exists(out_file),
  quote     = FALSE,
  append    = TRUE
)

if (err_flag) {
  cat(sprintf(
    "bG_index = %d (bG = %.4f), rep = %d: ERROR/TIMEOUT, lambda_hat = NA\n",
    idx, bG, rep_idx
  ))
} else {
  cat(sprintf("bG_index         = %d (bG = %.4f)\n", idx, bG))
  cat(sprintf("rep_index        = %d\n", rep_idx))
  cat(sprintf("lambda0          = %.4f\n", lambda0))
  cat(sprintf("lambda_hat       = %.4f\n", lambda_hat))
  cat(sprintf("D(lambda_hat)    = %.4f\n", D_value))

  cat("---- Linear regression p-values (PRS ~ covariates) ----\n")
  cat(sprintf("p_lr_y           = %.4g  (original scale)\n", p_lr_y))
  cat(sprintf("p_lr_log         = %.4g  (log(y))\n", p_lr_log))
  cat(sprintf("p_lr_latent      = %.4g  (SIQReg latent: lambda_hat)\n", p_lr_lat))
  cat(sprintf("p_lr_rint        = %.4g  (RINT(y))\n", p_lr_rint))
  cat(sprintf("p_lr_oracle      = %.4g  (oracle: true y*)\n", p_lr_oracle))

  cat("---- QR joint p-values (PRS ~ covariates) ----\n")
  cat(sprintf("joint_p_y        = %.4g  (original scale)\n", joint_p_y))
  cat(sprintf("joint_p_log      = %.4g  (log(y))\n", joint_p_log))
  cat(sprintf("joint_p_latent   = %.4g  (SIQReg latent)\n", joint_p_latent))
  cat(sprintf("joint_p_rint     = %.4g  (RINT(y))\n", joint_p_rint))
  cat(sprintf("joint_p_oracle   = %.4g  (oracle: true y*)\n", joint_p_oracle))
}
