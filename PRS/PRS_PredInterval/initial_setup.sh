#!/usr/bin/env bash
#SBATCH -J PI_setup
#SBATCH --partition=tier3q
#SBATCH --time=2:00:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/setup_%j.out
#SBATCH -e logs/setup_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

mkdir -p "${OUTBASE}"/{logs,original_scale,bc_scale}
mkdir -p "${OUTBASE}/PredInterval"

cd "${OUTBASE}"

if [[ ! -f "${OUTBASE}/PredInterval/src/PredInterval.R" ]]; then
    rm -rf "${OUTBASE}/PredInterval"
    git clone https://github.com/xuchang0201/PredInterval.git
fi

FAM_FILE="${GENO_MERGED}.fam"
TRAIN_FILE="${OUTBASE}/samples_train.txt"
TEST_FILE="${OUTBASE}/samples_test.txt"

if [[ ! -f "${FAM_FILE}" ]]; then
    echo "ERROR: FAM file not found: ${FAM_FILE}"
    exit 1
fi

module load gcc/12.1.0
module load R/4.4.1
module load python/3.10.5

Rscript --vanilla - <<RSCRIPT
library(data.table)

fam_file <- "${FAM_FILE}"
train_file <- "${TRAIN_FILE}"
test_file <- "${TEST_FILE}"
train_ratio <- ${TRAIN_RATIO}

fam <- fread(fam_file, header=FALSE)
colnames(fam)[1:2] <- c("FID", "IID")
cat("Total samples in FAM:", nrow(fam), "\n")

set.seed(2024)
n <- nrow(fam)
n_train <- round(n * train_ratio)

train_idx <- sample(1:n, n_train, replace=FALSE)
test_idx <- setdiff(1:n, train_idx)

train_samples <- fam[train_idx, .(FID, IID)]
test_samples <- fam[test_idx, .(FID, IID)]

cat("Train:", nrow(train_samples), "\n")
cat("Test:", nrow(test_samples), "\n")

fwrite(train_samples, train_file, sep="\t", col.names=FALSE)
fwrite(test_samples, test_file, sep="\t", col.names=FALSE)
RSCRIPT

echo "Done. Train: ${TRAIN_FILE}, Test: ${TEST_FILE}"
