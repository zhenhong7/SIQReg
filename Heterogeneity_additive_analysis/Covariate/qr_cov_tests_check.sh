#!/usr/bin/env bash
#SBATCH --job-name=qr_cov_meta
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --array=1-1075
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# R env
module load gcc/12.1.0
module load R/4.4.1

# Logs directory
mkdir -p logs

# Temp dir fix (avoids "cannot create 'R_TempDir'")
export TMPDIR="${SLURM_TMPDIR:-/tmp}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"; export TMP="$TMPDIR" TEMP="$TMPDIR"

RS="${SCRIPT_DIR}/qr_cov_tests_check.R"
LIST="${PHENO_DIR}/phenos.txt"       # 25 lines (one per phenotype)
COVARS="${PHENO_DIR}/covariates.txt" # 43 lines (one per covariate)

# Export R-accessible env vars
export PHENO_DIR="${PHENO_ORIGINAL_DIR}"
export COVAR_FILE="${COVAR_FILE}"
export PRS_DIR="${PRS_DIR}"
export TOP1_PRS_CSV="${TOP1_PRS_CSV}"
export META_BOTH_FILE="${META_BOTH_FILE}"
export KEEP_FILE="${ANCESTRY_DIR}/whitebrit.fid_iid.txt"
export OUT_BASE="${RESULTS_ROOT}/qr_cov_tests"
export BIRTH_YEAR_BASELINE="${BIRTH_YEAR_BASELINE}"

# Which lambda to use: mean or median (set via --export=MODE=mean|median; default median)
MODE="${MODE:-median}"

A=${SLURM_ARRAY_TASK_ID}
PHENO_IDX=$(( (A - 1) / 43 + 1 ))   # 1..25
COVAR_IDX=$(( (A - 1) % 43 + 1 ))   # 1..43

# Get names for this task
PHENO=$(sed -n "${PHENO_IDX}p" "$LIST"   | tr -d '\r' || true)
COVAR=$(sed -n "${COVAR_IDX}p" "$COVARS" | tr -d '\r' || true)

if [[ -z "${PHENO:-}" ]]; then
  echo "No phenotype for index ${PHENO_IDX} (array ${A})" >&2
  exit 1
fi
if [[ -z "${COVAR:-}" ]]; then
  echo "No covariate for index ${COVAR_IDX} (array ${A})" >&2
  exit 1
fi

echo "[$(date)] Running for PHENO=${PHENO}  COVAR=${COVAR}  MODE=${MODE}  (A=${A}, i=${PHENO_IDX}, j=${COVAR_IDX})"

# Run: pass MODE as the 3rd argument (mean|median) to select which meta lambda to use
Rscript "$RS" "$PHENO" "$COVAR" "$MODE"
