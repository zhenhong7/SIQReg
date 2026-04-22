#!/usr/bin/env bash
#SBATCH -J PI_gwas
#SBATCH --partition=tier1q
#SBATCH --time=24:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/gwas_%A_%a.out
#SBATCH -e logs/gwas_%A_%a.err
#SBATCH --array=1-22

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SCALE=${1:? "Usage: sbatch gwas.sh <scale> <trait_idx> <fold>"}
TRAIT_IDX=${2:? "Usage: sbatch gwas.sh <scale> <trait_idx> <fold>"}
FOLD=${3:? "Usage: sbatch gwas.sh <scale> <trait_idx> <fold>"}

CHR=${SLURM_ARRAY_TASK_ID:-1}

TRAIT_STR="${TRAIT_INFO[$TRAIT_IDX]}"
IFS=':' read -r phenoLower phenoFull phenoCode <<< "${TRAIT_STR}"

if [[ "${SCALE}" == "original" ]]; then
    OUTDIR="${OUTBASE}/original_scale/${phenoLower}"
else
    OUTDIR="${OUTBASE}/bc_scale/${phenoLower}"
fi

echo "GWAS: ${phenoLower}, Fold ${FOLD}, Chr ${CHR}"

GWAS_DIR="${OUTDIR}/gwas"
mkdir -p "${GWAS_DIR}"

GENO_FILE="${GENO_CHR_DIR}/${GENO_CHR_PREFIX}${CHR}${GENO_CHR_SUFFIX}"

# UPDATED: Use new naming convention
GWAS_TRAIN_FILE="${OUTDIR}/folds/gwas_train_fold_${FOLD}.txt"
PHENO_FILE="${OUTDIR}/folds/pheno_gwas_fold_${FOLD}.txt"

for f in "${GENO_FILE}.bed" "${GWAS_TRAIN_FILE}" "${PHENO_FILE}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

OUT_PREFIX="${GWAS_DIR}/chr${CHR}_fold${FOLD}"

if [[ -f "${OUT_PREFIX}.PHENO.glm.linear" ]]; then
    echo "Output exists, skipping"
else
    ${PLINK2} \
        --bfile "${GENO_FILE}" \
        --keep "${GWAS_TRAIN_FILE}" \
        --pheno "${PHENO_FILE}" \
        --pheno-name PHENO \
        --covar "${COVAR_FILE}" \
        --covar-col-nums "${COVAR_COLS}" \
        --covar-variance-standardize \
        --glm hide-covar \
        --out "${OUT_PREFIX}"

    echo "Done: ${OUT_PREFIX}.PHENO.glm.linear"
fi
