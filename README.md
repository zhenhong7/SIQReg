# SIQReg

Scale-Invariant Quantile Regression (SIQReg) estimates the latent phenotype scale with minimal scale-induced genetic heterogeneity.

## Quick start

```r
source("SIQReg.R")

result <- siqreg(y, C)
result$lambda_hat   # pooled Box-Cox lambda
result$se           # standard error
```

where `y` is a positive numeric vector (the phenotype) and `C` is a covariate matrix (e.g., sex, age, PCs). Each column of `C` is used as a target predictor in turn, and the per-covariate lambda estimates are meta-analyzed.

### Command line

```bash
Rscript SIQReg.R --y pheno.txt --c covariates.txt --out result.csv
```

Input files are headerless. `pheno.txt` has one value per line. `covariates.txt` is whitespace-delimited with one row per sample and one column per covariate. Both files must have the same number of rows. Rows with NA in any column are dropped.

Output is a one-row CSV with columns `lambda_hat` and `se`.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `taus` | 0.1, 0.2, ..., 0.9 | Quantile levels for the QR fits |
| `B` | 300 | Number of bootstrap replicates for the covariance of r(tau) |
| `n_chunks` | 10 | Number of random subsamples per covariate (controls SE estimation) |
| `lambda_lower` | -5 | Lower bound of the lambda search interval |
| `lambda_upper` | 5 | Upper bound of the lambda search interval |
| `conquer_threshold` | 1000 | Minimum subsample size to use `conquer` instead of `rq` |
| `seed` | NULL | Random seed for reproducibility |

## Dependencies

- **R** (>= 4.0)
- **MASS** (included with base R)
- **quantreg** (used when subsample size < `conquer_threshold`)
- **conquer** (used when subsample size >= `conquer_threshold`)

Install the optional packages with:

```r
install.packages(c("quantreg", "conquer"))
```

## Output

The `siqreg()` function returns a list with:

| Field | Description |
|-------|-------------|
| `lambda_hat` | Pooled Box-Cox lambda from DerSimonian-Laird meta-analysis |
| `se` | Standard error of the pooled lambda |
| `k` | Number of covariates used in the meta-analysis |
| `per_covariate` | Data frame with per-covariate lambda and SE |


## Repository structure

```
SIQReg.R                           # Public-facing lambda estimator
config.sh                          # Central path configuration for the analyses below
Simulation/                        # Simulation studies (lambda recovery, power, heavy tails)
Lambda_learning/                   # Lambda estimation pipelines (covariate-based and PRS-based)
Heterogeneity_additive_analysis/   # Downstream heterogeneity tests (GWAS, PRS, TWAS, vQTL, PRSxE)
PRS/                               # PRS prediction intervals and R-squared comparisons
```

