# SIQReg: Scale-Invariant Quantile Regression
#
# The method estimates lambda separately for each covariate (one at a time),
# then combines the per-covariate estimates via meta-analysis.
#
# Usage (as a function):
#   source("SIQReg.R")
#   result <- siqreg(y, C)
#   result$lambda_hat   # pooled lambda
#   result$se           # standard error
#
# Usage (command line):
#   Rscript SIQReg.R --y pheno.txt --c covariates.txt [options]
#
# Required inputs:
#   y   Outcome vector (must be positive)
#   C   Covariate matrix
#
# Output: lambda_hat and se.

suppressPackageStartupMessages({
  library(MASS)
})

siqreg <- function(y, C,
                   taus = seq(0.1, 0.9, by = 0.1),
                   B = 300L,
                   n_chunks = 10L,
                   lambda_lower = -5,
                   lambda_upper = 5,
                   conquer_threshold = 1000L,
                   seed = NULL,
                   verbose = TRUE) {

  y <- as.numeric(y)
  C <- as.matrix(C)
  if (nrow(C) != length(y))
    stop("C must have the same number of rows as length(y)")

  keep <- is.finite(y) & complete.cases(C)
  y <- y[keep]; C <- C[keep, , drop = FALSE]
  n <- length(y)

  if (n < 50) stop("Too few complete observations (n = ", n, ")")
  if (all(y <= 0)) stop("y must contain positive values for Box-Cox transformation")

  n_sub <- ceiling(n / n_chunks)
  K <- ncol(C)
  if (K < 1) stop("C must have at least one column")

  use_conquer <- (n_sub >= conquer_threshold)
  if (use_conquer) {
    if (!requireNamespace("conquer", quietly = TRUE))
      stop("Package 'conquer' is required for large samples. Install it or set conquer_threshold higher.")
  } else {
    if (!requireNamespace("quantreg", quietly = TRUE))
      stop("Package 'quantreg' is required for small samples. Install it with install.packages('quantreg').")
  }

  if (verbose) {
    method <- if (use_conquer) "conquer" else "rq"
    message(sprintf("SIQReg: n=%d, %d covariates, %d chunks, QR method=%s",
                    n, K, n_chunks, method))
  }

  # --- Core functions ---

  boxcox_transform <- function(y, lambda) {
    if (abs(lambda) < 1e-8) {
      log(y)
    } else {
      t <- lambda * log(y)
      t <- pmin(pmax(t, -700), 700)
      (exp(t) - 1) / lambda
    }
  }

  qr_fit <- function(Xmat, Ystar, tau) {
    if (use_conquer) {
      fit <- conquer::conquer(X = Xmat, Y = Ystar, tau = tau,
                              ci = "none", tol = 1e-6)
      as.numeric(fit$coeff[2])
    } else {
      dat <- data.frame(Y = Ystar, Xmat)
      fit <- suppressWarnings(quantreg::rq(Y ~ ., data = dat, tau = tau))
      as.numeric(coef(fit)[2])
    }
  }

  slopes_at_lambda <- function(lambda, y, G, C_rest, taus) {
    Ystar <- boxcox_transform(y, lambda)
    if (any(!is.finite(Ystar))) return(rep(NA_real_, length(taus)))

    m <- median(Ystar, na.rm = TRUE)
    s <- mad(Ystar, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s < 1e-8) s <- sd(Ystar)
    if (!is.finite(s) || s < 1e-8) s <- 1.0
    Ystar <- (Ystar - m) / s

    Xmat <- if (is.null(C_rest)) cbind(G) else cbind(G, C_rest)

    sapply(taus, function(tau) {
      tryCatch(qr_fit(Xmat, Ystar, tau), error = function(e) NA_real_)
    })
  }

  D_ratio <- function(lambda, y, G, C_rest, taus, B, seed) {
    Kt <- length(taus)
    beta0 <- slopes_at_lambda(lambda, y, G, C_rest, taus)
    if (any(is.na(beta0))) return(1e8)

    beta_bar <- mean(beta0)
    if (!is.finite(beta_bar) || abs(beta_bar) < 1e-12) return(1e8)
    r0 <- 1 - beta0 / beta_bar

    if (!is.null(seed)) set.seed(seed)

    nn <- length(y)
    boot_r <- replicate(B, {
      idx <- sample.int(nn, nn, replace = TRUE)
      bvec <- slopes_at_lambda(lambda, y[idx], G[idx],
                               if (is.null(C_rest)) NULL else C_rest[idx, , drop = FALSE],
                               taus)
      if (any(is.na(bvec))) return(rep(NA_real_, Kt))
      bbar <- mean(bvec)
      if (!is.finite(bbar) || abs(bbar) < 1e-12) return(rep(NA_real_, Kt))
      1 - bvec / bbar
    })

    Rmat <- t(boot_r)
    if (!is.matrix(Rmat) || all(!is.finite(Rmat))) return(1e8)

    ok <- complete.cases(Rmat)
    if (sum(ok) < Kt + 5) return(1e8)
    Rmat <- Rmat[ok, , drop = FALSE]

    Sigma_r <- cov(Rmat, use = "complete.obs")
    invR <- tryCatch(solve(Sigma_r), error = function(e) NULL)
    if (is.null(invR)) invR <- MASS::ginv(Sigma_r)

    as.numeric(t(r0) %*% invR %*% r0)
  }

  estimate_lambda_one <- function(y, G, C_rest, taus, B,
                                  lambda_lower, lambda_upper, seed) {
    objfun <- function(lam) {
      if (!is.null(seed)) set.seed(seed)
      D_ratio(lam, y, G, C_rest, taus, B, seed)
    }
    sol <- optimize(objfun, lower = lambda_lower, upper = lambda_upper, tol = 1e-4)
    sol$minimum
  }

  # --- DerSimonian-Laird random-effects meta-analysis ---

  dl_meta <- function(theta, se) {
    w <- 1 / se^2
    theta_fe <- sum(w * theta) / sum(w)
    Q <- sum(w * (theta - theta_fe)^2)
    k <- length(theta)
    tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
    w_re <- 1 / (se^2 + tau2)
    lambda_dl <- sum(w_re * theta) / sum(w_re)
    se_dl <- sqrt(1 / sum(w_re))
    list(lambda_hat = lambda_dl, se = se_dl, k = k)
  }

  # --- Per-covariate lambda estimation via subsampling ---

  if (!is.null(seed)) set.seed(seed)
  perm <- sample.int(n)
  chunk_id <- rep(seq_len(n_chunks), each = n_sub, length.out = n)

  per_cov_lambda <- numeric(K)
  per_cov_se <- numeric(K)
  cov_names <- colnames(C)
  if (is.null(cov_names)) cov_names <- paste0("C", seq_len(K))

  for (j in seq_len(K)) {
    G_full <- C[, j]
    C_rest_full <- if (K == 1) NULL else C[, -j, drop = FALSE]

    chunk_lambdas <- numeric(n_chunks)
    for (ci in seq_len(n_chunks)) {
      idx <- perm[chunk_id == ci]

      seed_chunk <- if (!is.null(seed)) {
        as.integer(abs(seed + 1009L * ci + 17L * j) %% .Machine$integer.max)
      } else NULL

      chunk_lambdas[ci] <- estimate_lambda_one(
        y[idx], G_full[idx],
        if (is.null(C_rest_full)) NULL else C_rest_full[idx, , drop = FALSE],
        taus, B, lambda_lower, lambda_upper, seed_chunk
      )
    }

    per_cov_lambda[j] <- mean(chunk_lambdas)
    per_cov_se[j] <- sd(chunk_lambdas) / sqrt(n_chunks)

    if (verbose) {
      message(sprintf("  Covariate [%s]: lambda=%.4f, SE=%.4f",
                      cov_names[j], per_cov_lambda[j], per_cov_se[j]))
    }
  }

  valid <- is.finite(per_cov_lambda) & is.finite(per_cov_se) & per_cov_se > 0
  if (sum(valid) < 1) stop("No valid per-covariate lambda estimates obtained")

  if (sum(valid) == 1) {
    result <- list(
      lambda_hat = per_cov_lambda[valid],
      se = per_cov_se[valid],
      k = 1L,
      per_covariate = data.frame(
        covariate = cov_names, lambda = per_cov_lambda, se = per_cov_se,
        stringsAsFactors = FALSE
      )
    )
  } else {
    meta <- dl_meta(per_cov_lambda[valid], per_cov_se[valid])
    result <- c(meta, list(
      per_covariate = data.frame(
        covariate = cov_names, lambda = per_cov_lambda, se = per_cov_se,
        stringsAsFactors = FALSE
      )
    ))
  }

  if (verbose) {
    message(sprintf("  Meta-analysis (k=%d): lambda=%.4f, SE=%.4f",
                    result$k, result$lambda_hat, result$se))
  }

  result
}


# --- Command-line interface ---
if (!interactive() && length(sys.frames()) == 0L) {

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    cat("SIQReg: Scale-Invariant Quantile Regression lambda estimator\n\n")
    cat("Usage:\n")
    cat("  Rscript SIQReg.R --y pheno.txt --c covariates.txt [options]\n\n")
    cat("Required:\n")
    cat("  --y FILE        Outcome vector (one value per line, must be positive)\n")
    cat("  --c FILE        Covariate matrix (whitespace-delimited, no header)\n\n")
    cat("Optional:\n")
    cat("  --taus STR      Comma-separated quantile levels (default: 0.1,0.2,...,0.9)\n")
    cat("  --B INT         Bootstrap replicates (default: 300)\n")
    cat("  --n_chunks INT  Number of subsamples per covariate (default: 10)\n")
    cat("  --lambda_lower  Lower bound for lambda search (default: -5)\n")
    cat("  --lambda_upper  Upper bound for lambda search (default: 5)\n")
    cat("  --conquer_threshold INT  Min subsample size to use conquer (default: 2000)\n")
    cat("  --seed INT      Random seed for reproducibility\n")
    cat("  --out FILE      Output CSV file (default: stdout)\n")
    quit(save = "no", status = 0)
  }

  parse_args <- function(args) {
    opts <- list()
    i <- 1
    while (i <= length(args)) {
      key <- sub("^--", "", args[i])
      if (i + 1 <= length(args)) {
        opts[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        stop("Missing value for argument: ", args[i])
      }
    }
    opts
  }

  opts <- parse_args(args)

  if (is.null(opts$y)) stop("--y is required")
  if (is.null(opts$c)) stop("--c is required")

  y_vec <- scan(opts$y, what = numeric(), quiet = TRUE)
  c_mat <- as.matrix(read.table(opts$c, header = FALSE))

  taus <- seq(0.1, 0.9, by = 0.1)
  if (!is.null(opts$taus)) taus <- as.numeric(strsplit(opts$taus, ",")[[1]])

  B <- 300L
  if (!is.null(opts$B)) B <- as.integer(opts$B)

  n_chunks <- 10L
  if (!is.null(opts$n_chunks)) n_chunks <- as.integer(opts$n_chunks)

  lambda_lower <- -5
  if (!is.null(opts$lambda_lower)) lambda_lower <- as.numeric(opts$lambda_lower)

  lambda_upper <- 5
  if (!is.null(opts$lambda_upper)) lambda_upper <- as.numeric(opts$lambda_upper)

  conquer_threshold <- 2000L
  if (!is.null(opts$conquer_threshold)) conquer_threshold <- as.integer(opts$conquer_threshold)

  seed <- NULL
  if (!is.null(opts$seed)) seed <- as.integer(opts$seed)

  result <- siqreg(y_vec, c_mat,
                   taus = taus, B = B, n_chunks = n_chunks,
                   lambda_lower = lambda_lower, lambda_upper = lambda_upper,
                   conquer_threshold = conquer_threshold, seed = seed)

  out_df <- data.frame(lambda_hat = result$lambda_hat, se = result$se)

  if (!is.null(opts$out)) {
    write.csv(out_df, opts$out, row.names = FALSE)
    message("Wrote: ", opts$out)
  } else {
    write.csv(out_df, "", row.names = FALSE)
  }
}
