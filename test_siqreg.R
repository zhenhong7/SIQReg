# Simulation: validate SIQReg
# Usage: Rscript test_siqreg.R

source("SIQReg.R")

lambda_vals <- c(0, 0.5, 1)
n <- 5000
n_reps <- 10
b0 <- -0.1
aT <- 0.50
eps_sd <- 0.06

set.seed(999)
tgt <- runif(n)
C_mat <- matrix(tgt, ncol = 1)
mu <- b0 + aT * tgt

results <- data.frame(lambda0 = numeric(), rep = integer(),
                      lambda_hat = numeric(), se = numeric(),
                      lambda_full = numeric())

for (lam0 in lambda_vals) {
  for (r in 1:n_reps) {
    set.seed(1000 * lam0 + r)
    eps <- rnorm(n, sd = eps_sd)
    y_star <- mu + eps

    if (abs(lam0) >= 1e-8) {
      bound <- -1 / lam0
      for (try in 1:200) {
        bad <- if (lam0 > 0) (y_star <= bound + 1e-6) else (y_star >= bound - 1e-6)
        if (!any(bad)) break
        eps[bad] <- rnorm(sum(bad), sd = eps_sd)
        y_star[bad] <- mu[bad] + eps[bad]
      }
    }

    y <- if (abs(lam0) < 1e-8) exp(y_star) else (1 + lam0 * y_star)^(1 / lam0)
    if (any(!is.finite(y)) || any(y <= 0)) { message("skip"); next }

    res <- siqreg(y, C_mat, B = 300, n_chunks = 10, seed = 42,
                  full_sample = TRUE, verbose = TRUE)

    results <- rbind(results,
      data.frame(lambda0 = lam0, rep = r,
                 lambda_hat = res$lambda_hat, se = res$se,
                 lambda_full = res$lambda_full))

    message(sprintf("lambda0=%.1f rep=%d: chunk=%.4f, full=%.4f, se=%.4f",
                    lam0, r, res$lambda_hat, res$lambda_full, res$se))
  }
}

cat("\n===== Summary =====\n")
for (lam0 in lambda_vals) {
  sub <- results[results$lambda0 == lam0, ]
  cat(sprintf("lambda0=%.1f  chunk_mean:  mean=%.4f, sd=%.4f, bias=%+.4f\n",
              lam0, mean(sub$lambda_hat), sd(sub$lambda_hat),
              mean(sub$lambda_hat) - lam0))
  cat(sprintf("lambda0=%.1f  full_sample: mean=%.4f, sd=%.4f, bias=%+.4f\n",
              lam0, mean(sub$lambda_full), sd(sub$lambda_full),
              mean(sub$lambda_full) - lam0))
}
