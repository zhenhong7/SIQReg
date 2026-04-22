#!/usr/bin/env bash
#SBATCH --job-name=TWAS
#SBATCH --time=240:00:00
#SBATCH --partition=tier1q
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

PHENO="${1:?Usage: run_twas_array.sh <PHENO>}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

module load gcc/12.1.0
module load R/4.4.1

# Export R-accessible env vars
export PHENO_DIR="${PHENO_ORIGINAL_DIR}"
export COVAR_FILE="${COVAR_FILE}"
export EXPR_DIR="${EXPR_DIR}"
export KEEP_FILE="${ANCESTRY_DIR}/whitebrit.fid_iid.txt"
export META_BOTH_FILE="${META_BOTH_FILE}"
export TWAS_RESULTS="${TWAS_RESULTS}"
export BIRTH_YEAR_BASELINE="${BIRTH_YEAR_BASELINE}"

TISSUE=$(awk -v p="$PHENO" '$1==p {print $2}' "$PHENO2TISSUE")
if [[ -z "$TISSUE" ]]; then
  echo "No tissue found for phenotype '$PHENO' in $PHENO2TISSUE"; exit 1
fi

MAP="${TWAS_MAPPING_DIR}/${TISSUE}/mapping.txt"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MAP" || true)
if [[ -z "$LINE" ]]; then
  echo "Empty mapping line for task ${SLURM_ARRAY_TASK_ID}"; exit 0
fi
CHR=$(echo "$LINE" | awk '{print $1}')
GENE_INDEX=$(echo "$LINE" | awk '{print $2}')

OUTBASE="${TWAS_RESULTS}/${PHENO}/${TISSUE}"
mkdir -p "${OUTBASE}/logs"

echo "PHENO=${PHENO}  TISSUE=${TISSUE}  CHR=${CHR}  GENE_INDEX=${GENE_INDEX}"

# 1) Transformed scale
Rscript "${SCRIPT_DIR}/TWAS_tr.R"  "$CHR" "$GENE_INDEX" "$PHENO" "$TISSUE"

# 2) Raw scale
Rscript "${SCRIPT_DIR}/TWAS_raw.R" "$CHR" "$GENE_INDEX" "$PHENO" "$TISSUE"
