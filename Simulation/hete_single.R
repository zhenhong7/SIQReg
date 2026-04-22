#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(conquer))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(wesanderson))
library(ggplot2)

lambda0 <- 1/2
n       <- 5000L
pW      <- 9L
taus    <- seq(0.1, 0.9, by = 0.1)

b0     <- 0
bG     <- 0.5
aT     <- 0.5
aW     <- rep(0, pW)
eps_sd <- 0.02

qr_betas_with_ci <- function(Y, G, taus) {
  Y <- as.numeric(Y); G <- as.numeric(G)
  good <- is.finite(Y) & is.finite(G)
  Y <- Y[good]; G <- G[good]
  X <- cbind(G = G)
  res <- data.frame(tau = taus, beta = NA_real_,
                    ci_lower = NA_real_, ci_upper = NA_real_)
  for (i in seq_along(taus)) {
    fit <- try(conquer(X = X, Y = Y, tau = taus[i], ci = "asymptotic"), silent = TRUE)
    if (inherits(fit, "try-error") || length(fit$coeff) == 0L ||
        is.null(fit$asyCI) || length(fit$asyCI) == 0L) next
    res$beta[i]     <- fit$coeff[2]
    res$ci_lower[i] <- fit$asyCI[2, 1]
    res$ci_upper[i] <- fit$asyCI[2, 2]
  }
  res
}

rint_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

scale_to_unit <- function(z) {
  z <- as.numeric(z)
  mu <- mean(z, na.rm = TRUE)
  sdv <- sd(z, na.rm = TRUE)
  if (!is.finite(sdv) || sdv < 1e-8) sdv <- 1
  (z - mu) / sdv
}

set.seed(123456)
PRS <- rbinom(n, size = 1, prob = 0.5)
tgt <- runif(n)
W   <- replicate(pW, runif(n))
eps <- rnorm(n, 0, eps_sd)

y_star_raw  <- b0 + bG * PRS + aT * tgt + as.vector(W %*% aW) + eps
shift_const <- 0.49 - max(y_star_raw)
y_star      <- y_star_raw + shift_const
sd_y_star   <- sd(y_star)
bG_std_true <- bG / sd_y_star

y      <- (1 + lambda0 * y_star)^(1 / lambda0)
log_y  <- log(y)
rint_y <- rint_transform(y)

y_std      <- scale_to_unit(y)
log_y_std  <- scale_to_unit(log_y)
rint_y_std <- scale_to_unit(rint_y)

## Box-Cox lambda estimate
dat <- data.frame(y = y, PRS = PRS, tgt = tgt)
bc  <- boxcox(y ~ PRS + tgt, data = dat,
              lambda = seq(-4, 2, by = 0.1), plotit = FALSE)
lambda_hat <- bc$x[which.max(bc$y)]

y_star_hat_raw <- if (abs(lambda_hat) < 1e-8) log(y) else (y^lambda_hat - 1) / lambda_hat
y_star_hat     <- scale_to_unit(y_star_hat_raw)

bet_latent <- qr_betas_with_ci(y_star_hat, PRS, taus)
bet_latent$scale <- "latent"
bet_orig <- qr_betas_with_ci(y_std, PRS, taus)
bet_orig$scale <- "original"
bet_log <- qr_betas_with_ci(log_y_std, PRS, taus)
bet_log$scale <- "log"
bet_rint <- qr_betas_with_ci(rint_y_std, PRS, taus)
bet_rint$scale <- "rint"

res_all <- rbind(bet_latent, bet_orig, bet_log, bet_rint)
res_all$scale <- factor(
  res_all$scale,
  levels = c("latent", "original", "log", "rint"),
  labels = c("SIQReg scale", "Original scale", "Log scale", "RINT scale")
)

lr_df <- data.frame(
  scale   = factor(c("SIQReg scale", "Original scale", "Log scale", "RINT scale"),
                   levels = levels(res_all$scale)),
  lr_beta = c(coef(lm(y_star_hat ~ PRS))["PRS"],
              coef(lm(y_std ~ PRS))["PRS"],
              coef(lm(log_y_std ~ PRS))["PRS"],
              coef(lm(rint_y_std ~ PRS))["PRS"])
)

pal <- wes_palette("AsteroidCity1", 4, type = "discrete")
cols <- c(
  "Original scale" = pal[2],
  "Log scale"      = pal[1],
  "SIQReg scale"   = pal[3],
  "RINT scale"     = pal[4]
)

y_range  <- range(c(res_all$beta, bG_std_true), na.rm = TRUE)
y_breaks <- pretty(y_range, n = 5)
if (!any(abs(y_breaks - bG_std_true) < 1e-8)) {
  y_breaks <- sort(c(y_breaks, bG_std_true))
}
y_labels <- sapply(y_breaks, function(v) {
  if (abs(v - bG_std_true) < 1e-8) "oracle" else format(v, digits = 2)
})

p <- ggplot(res_all, aes(x = tau, y = beta,
                         colour = scale, group = scale)) +
  
  # QR lines + points + CI
  geom_line(linewidth = 0.7) +
  geom_point(size = 3.5, shape = 16) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.02, linewidth = 0.4) +
  
  # Oracle true effect (dotted)
  geom_hline(yintercept = bG_std_true,
             linetype = "dotted", linewidth = 0.7) +
  
  scale_x_continuous(breaks = taus) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  
  scale_colour_manual(
    values = cols,
    labels = c(
      "Original scale" = expression("Original scale (" * lambda == 1 * ")"),
      "Log scale"      = expression("Log scale (" * lambda == 0 * ")"),
      "SIQReg scale"   = expression("SIQReg scale (" * lambda == hat(lambda)[0] * ")"),
      "RINT scale"     = expression("RINT scale")
    )
  ) +
  
  labs(
    x      = expression("Quantile level: " * tau),
    y      = "Genetic effect size",
    colour = NULL
  ) +
  
  annotate("text",
           x = max(taus), y = min(y_breaks),
           label = expression(lambda[0] == 0.5),
           hjust = 1.1, vjust = -0.5, size = 6) +
  
  theme_classic(base_size = 10) +
  theme(
    legend.position      = c(0.5, 0.9),
    legend.text          = element_text(size = 14),
    legend.background    = element_rect(fill = "white", colour = "grey70", linewidth = 0.3),
    axis.title           = element_text(size = 18),
    axis.text            = element_text(size = 16, colour = "black"),
    plot.margin          = margin(5, 5, 5, 5)
  )

print(p)


ggsave("fig2f.pdf", p, width=8, height=8)