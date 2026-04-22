library("dplyr")

## INPUTS (unchanged mindset)
PHENO_DIGITS <- commandArgs(trailingOnly = TRUE)[1]
PHENO_NO_DIGITS <- ifelse(grepl("EA4", PHENO_DIGITS), "EA4",
                          ifelse(grepl("FEV1674178", PHENO_DIGITS), "FEV1",
                                 ifelse(grepl("IGF-1674178", PHENO_DIGITS), "IGF-1",
                                        gsub("[0-9]+$", "", PHENO_DIGITS))))
PHENO_LOWER <- tolower(PHENO_NO_DIGITS)

PHENOS <- Sys.getenv("PHENO_DIR", unset = "/path/to/extracted_phenotypes/")
TEST_POP <- "whitebrit"
PGS_COLS <- c("FID","IID","ALLELE_CT","ALLELE_DOSAGE_SUM","SCORE1_AVG","SCORE1_SUM")
THRESH <- c("0.0000000001", "0.00000001", "0.000001", "0.0001", "0.001",
            "0.005", "0.01", "0.05", "0.1", "0.5")
PHENO_FILE <- paste0(PHENOS, PHENO_NO_DIGITS, "/", PHENO_DIGITS, ".pheno")
R2_TEST <- "PEARSON"

# Crash if R2 is not spearman or pearson
if (!(R2_TEST %in% c("SPEARMAN", "PEARSON"))) {
  stop("R2_TEST must be either 'SPEARMAN' or 'PEARSON'")
}

## trait-specific lambda
META_TAB <- Sys.getenv("META_TAB", unset = "/path/to/SI_GWAS/meta_both_summary_all.txt")
meta <- read.table(META_TAB, header=TRUE, sep="\t", stringsAsFactors=FALSE)
ii <- match(PHENO_NO_DIGITS, meta$phenotype)
if (is.na(ii)) stop("Cannot find phenotype in meta table: ", PHENO_NO_DIGITS)
lambda <- meta$lambda_WM[ii]
cat("pheno:", PHENO_LOWER, " lambda_WM:", lambda, "\n")

## OUTPUT PATH
R2_BASE <- Sys.getenv("R2_BASE", unset = "/path/to/SI_GWAS/results/pgs_outputs/")
PGS_BASE <- Sys.getenv("PGS_BASE", unset = "/path/to/SI_GWAS/pgs_outputs/pgs_out/")

if (lambda == 0) {
  if (R2_TEST == "PEARSON") {
    R2 <- paste0(R2_BASE, "r2_out_log/")
  } else {
    R2 <- paste0(R2_BASE, "r2_out_spearman_log/")
  }
  PGS <- paste0(PGS_BASE, "../pgs_out_log/")
} else {
  if (R2_TEST == "PEARSON") {
    R2 <- paste0(R2_BASE, "r2_out/")
  } else {
    R2 <- paste0(R2_BASE, "r2_out_spearman/")
  }
  PGS <- PGS_BASE
}

## SAME functions (unchanged)
R2fxn <- function(x,y,X){
  X <- as.matrix(X)
  y <- resid(lm(y ~ X, na.action = na.exclude))
  if (length(x) == length(y)) {
    x <- resid(lm(x ~ X, na.action = na.exclude))
  } else {
    for (j in 1:ncol(x)) x[,j] <- resid(lm(x[,j] ~ X, na.action = na.exclude))
  }
  summary(lm(y ~ 1 + x))
}

R2fxn_spearman <- function(x,y,X){
  X <- as.matrix(X)
  y <- resid(lm(y ~ X, na.action = na.exclude))
  if (length(x) == length(y)) {
    x <- resid(lm(x ~ X, na.action = na.exclude))
  } else {
    for (j in 1:ncol(x)) x[,j] <- resid(lm(x[,j] ~ X, na.action = na.exclude))
  }
  (cor(x, y, method="spearman", use="complete.obs"))^2
}

## inverse Box-Cox for PRS
inv_boxcox_prs <- function(z, lam) {
  z <- as.numeric(z)
  if (abs(lam) < 1e-12) return(exp(z))
  u <- 1 + lam * z
  out <- rep(NA_real_, length(z))
  ok <- is.finite(u) & (u > 0)
  out[ok] <- u[ok]^(1/lam)
  out
}

## rank metrics (minimal)
resid_on_X <- function(v, X){
  resid(lm(v ~ as.matrix(X), na.action = na.exclude))
}

R2_rank_lin <- function(y, x, X){
  ry <- rank(y, ties.method="average")
  rx <- rank(x, ties.method="average")
  cor(resid_on_X(ry, X), resid_on_X(rx, X), use="complete.obs")^2
}

R2_rank_EC_shared_linear <- function(y, x, X){
  ry <- rank(y, ties.method="average")
  rx <- rank(x, ties.method="average")
  fit <- lm(ry ~ as.matrix(X), na.action = na.exclude)
  gX  <- fitted(fit)
  cor(ry - gX, rx - gX, use="complete.obs")^2
}

## workflow
pheno_table <- read.table(PHENO_FILE, header=TRUE)

dir.create(R2, showWarnings=FALSE, recursive=TRUE)
dir.create(paste0(R2,"all_thresh/"), showWarnings=FALSE, recursive=TRUE)

out_file_rankLin  <- paste0(R2,"processed_",PHENO_LOWER,"_",TEST_POP,"_PGS_rankLin.txt")
out_file_rankEC   <- paste0(R2,"processed_",PHENO_LOWER,"_",TEST_POP,"_PGS_rankEC.txt")
out_file_std      <- paste0(R2,"processed_",PHENO_LOWER,"_",TEST_POP,"_PGS_std.txt")
out_file_spearman <- paste0(R2,"processed_",PHENO_LOWER,"_",TEST_POP,"_PGS_spearman.txt")
all_thresh        <- paste0(R2,"all_thresh/all_thresh_",PHENO_LOWER,"_",TEST_POP,"_PGS.txt")

if (file.exists(out_file_rankLin) && file.exists(out_file_rankEC) &&
    file.exists(out_file_std) && file.exists(out_file_spearman)) {
  quit(save="no", status=0)
}

COVAR_FILE <- Sys.getenv("COVAR_FILE",
  unset = paste0(PHENOS, "covariates_sa40PC/covariates_sa40PC674178.pheno"))
covars <- read.table(COVAR_FILE, header=TRUE)
covs_table <- covars[, 1:14]
rm(covars)

PGS_file <- paste0(PGS, PHENO_LOWER, "/", PHENO_LOWER, "_pgs.", THRESH[1], ".sscore")
PGS_IDs  <- read.table(PGS_file, header=FALSE, col.names=PGS_COLS)[,1:2]

covs  <- covs_table[covs_table$FID %in% PGS_IDs$FID, ]
pheno <- pheno_table[pheno_table$FID %in% PGS_IDs$FID, ]
pheno_cov <- merge(pheno, covs, by=c("FID", "IID"))
colnames(pheno_cov)[3] <- "pheno_code"
rm(PGS_IDs, PGS_file, covs_table, pheno)

y <- pheno_cov$pheno_code
X <- pheno_cov[, !colnames(pheno_cov) %in% c("FID","IID","pheno_code")]

prs.result <- NULL
for (i in THRESH) {
  prs_table <- read.table(paste0(PGS, PHENO_LOWER, "/", PHENO_LOWER, "_pgs.", i, ".sscore"),
                          header=FALSE, col.names=PGS_COLS)[, c(1,6)]
  prs <- merge(pheno_cov, prs_table, by="FID")$SCORE1_SUM

  prs <- inv_boxcox_prs(prs, lambda)

  good <- is.finite(prs) & is.finite(y)
  if (sum(good) == 0L || sd(prs[good]) == 0) {
    r2_std      <- NA_real_
    r2_spearman <- NA_real_
    r2_rankLin  <- NA_real_
    r2_rankEC   <- NA_real_
    w0_std      <- NA_real_
  } else {
    # Keep existing behavior for r2_std (Pearson if R2_TEST="PEARSON", Spearman if set)
    if (R2_TEST == "SPEARMAN") {
      r2_std <- R2fxn_spearman(x=prs, y=y, X=X)
      w0_std <- NA_real_
    } else {
      mod_std <- R2fxn(x=prs, y=y, X=X)
      r2_std  <- mod_std$r.squared
      w0_std  <- mod_std$coef["x","Estimate"]
    }

    # Always compute Spearman^2 on residuals as an extra metric
    r2_spearman <- R2fxn_spearman(x=prs, y=y, X=X)

    r2_rankLin <- R2_rank_lin(y=y, x=prs, X=X)
    r2_rankEC  <- R2_rank_EC_shared_linear(y=y, x=prs, X=X)
  }

  print(c(thresh=i, R2_std=r2_std, R2_spearman=r2_spearman, R2_rankLin=r2_rankLin, R2_rankEC=r2_rankEC))

  prs.result <- rbind(prs.result,
                      data.frame(Threshold=i,
                                 R2_std=r2_std,
                                 R2_spearman=r2_spearman,
                                 R2_rankLin=r2_rankLin,
                                 R2_rankEC=r2_rankEC,
                                 w0_std=w0_std,
                                 type="pgs"))
}

## best threshold under each metric
best_std      <- prs.result[which.max(prs.result$R2_std),      ]
best_spearman <- prs.result[which.max(prs.result$R2_spearman), ]
best_rankLin  <- prs.result[which.max(prs.result$R2_rankLin),  ]
best_rankEC   <- prs.result[which.max(prs.result$R2_rankEC),   ]

print(best_std)
print(best_spearman)
print(best_rankLin)
print(best_rankEC)

write.table(prs.result,    file=all_thresh,        row.names=FALSE, quote=FALSE)
write.table(best_std,      file=out_file_std,      row.names=FALSE, quote=FALSE)
write.table(best_spearman, file=out_file_spearman, row.names=FALSE, quote=FALSE)
write.table(best_rankLin,  file=out_file_rankLin,  row.names=FALSE, quote=FALSE)
write.table(best_rankEC,   file=out_file_rankEC,   row.names=FALSE, quote=FALSE)
