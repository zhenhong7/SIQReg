#!/usr/bin/env bash
#SBATCH -J PI_folds
#SBATCH --partition=tier1q
#SBATCH --time=2:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/folds_%A_%a.out
#SBATCH -e logs/folds_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SCALE=${1:? "Usage: sbatch setup_folds.sh <scale>  (original or bc)"}

if [[ "${SCALE}" != "original" && "${SCALE}" != "bc" ]]; then
    echo "ERROR: scale must be 'original' or 'bc'"
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
TRAIT_STR="${TRAIT_INFO[$TASK_ID]}"
IFS=':' read -r phenoLower phenoFull phenoCode <<< "${TRAIT_STR}"

if [[ "${SCALE}" == "original" ]]; then
    OUTDIR="${OUTBASE}/original_scale/${phenoLower}"
    PHENO_FILE="${PHENO_ORIGINAL_DIR}/${phenoFull}/${phenoFull}${phenoCode}.pheno"
else
    OUTDIR="${OUTBASE}/bc_scale/${phenoLower}"
    PHENO_FILE="${PHENO_BC_DIR}/${phenoFull}${phenoCode}_bc.pheno"
fi

mkdir -p "${OUTDIR}"/{folds,gwas,pgs,logs}

if [[ ! -f "${PHENO_FILE}" ]]; then
    echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
    exit 1
fi

TRAIN_SAMPLES="${OUTBASE}/samples_train.txt"
TEST_SAMPLES="${OUTBASE}/samples_test.txt"

for f in "${TRAIN_SAMPLES}" "${TEST_SAMPLES}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Sample file not found: $f. Run initial_setup.sh first."
        exit 1
    fi
done

module load gcc/12.1.0
module load R/4.4.1

Rscript --vanilla - <<RSCRIPT
library(data.table)

pheno_file <- "${PHENO_FILE}"
train_file <- "${TRAIN_SAMPLES}"
test_file <- "${TEST_SAMPLES}"
outdir <- "${OUTDIR}"
K <- ${K}
scale <- "${SCALE}"

pheno <- fread(pheno_file)
orig_cols <- colnames(pheno)
if (grepl("^#", orig_cols[1])) setnames(pheno, orig_cols[1], sub("^#", "", orig_cols[1]))

if (ncol(pheno) >= 3) {
    colnames(pheno)[1:2] <- c("FID", "IID")
    pheno_col <- colnames(pheno)[3]
    pheno <- pheno[, .(FID, IID, PHENO = get(pheno_col))]
}

pheno[, PHENO := as.numeric(as.character(PHENO))]
pheno <- pheno[!is.na(PHENO)]

train_ids <- fread(train_file, header=FALSE, col.names=c("FID", "IID"))
test_ids <- fread(test_file, header=FALSE, col.names=c("FID", "IID"))

pheno_train_pool <- pheno[IID %in% train_ids[["IID"]]]
pheno_test <- pheno[IID %in% test_ids[["IID"]]]

cat("Train:", nrow(pheno_train_pool), "Test:", nrow(pheno_test), "\n")

pheno_mean <- mean(pheno_train_pool[["PHENO"]], na.rm=TRUE)
pheno_sd <- sd(pheno_train_pool[["PHENO"]], na.rm=TRUE)

stats <- data.table(
    scale = scale,
    mean = pheno_mean,
    sd = pheno_sd,
    n_train_pool = nrow(pheno_train_pool),
    n_test = nrow(pheno_test)
)
fwrite(stats, paste0(outdir, "/pheno_stats.txt"), sep="\t")

set.seed(2024)
pheno_train_pool[, FOLD := sample(rep(1:K, length.out=.N))]
fwrite(pheno_train_pool[, .(FID, IID, FOLD)], paste0(outdir, "/fold_assignments.txt"), sep="\t")

for (k in 1:K) {
    gwas_train_k <- pheno_train_pool[FOLD != k, .(FID, IID)]
    cv_valid_k <- pheno_train_pool[FOLD == k, .(FID, IID)]

    fwrite(gwas_train_k, paste0(outdir, "/folds/gwas_train_fold_", k, ".txt"),
           sep="\t", col.names=FALSE)
    fwrite(cv_valid_k, paste0(outdir, "/folds/cv_valid_fold_", k, ".txt"),
           sep="\t", col.names=FALSE)

    pheno_gwas <- pheno_train_pool[FOLD != k, .(FID, IID, PHENO)]
    fwrite(pheno_gwas, paste0(outdir, "/folds/pheno_gwas_fold_", k, ".txt"), sep="\t")

    cat("Fold", k, ": GWAS train =", nrow(gwas_train_k),
        ", CV valid =", nrow(cv_valid_k), "\n")
}

cv_pheno <- pheno_train_pool[, .(id = IID, y = PHENO, fold = FOLD)]
fwrite(cv_pheno, paste0(outdir, "/cv_phenotypes.txt"), sep="\t")

test_pheno <- pheno_test[, .(id = IID, y = PHENO)]
fwrite(test_pheno, paste0(outdir, "/test_phenotypes.txt"), sep="\t")

fwrite(pheno_train_pool[, .(FID, IID)],
       paste0(outdir, "/train_pool_samples.txt"), sep="\t", col.names=FALSE)

fwrite(pheno_test[, .(FID, IID)],
       paste0(outdir, "/test_samples.txt"), sep="\t", col.names=FALSE)

dir.create(paste0(outdir, "/pgs"), showWarnings=FALSE)
fwrite(pheno_test[, .(FID, IID)],
       paste0(outdir, "/pgs/test.fam"), sep="\t", col.names=FALSE)

RSCRIPT
