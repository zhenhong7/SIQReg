#!/usr/bin/env bash
#SBATCH --job-name=qr_prs_check
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-900
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

module load gcc/12.1.0
module load R/4.4.1

mkdir -p logs
export TMPDIR="${SLURM_TMPDIR:-/tmp}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$TMPDIR"; export TMP="$TMPDIR" TEMP="$TMPDIR"

RS="${SCRIPT_DIR}/qr_prs_test_check.R"
PHENO_LIST="${PHENO_DIR}/phenos.txt"     # 25 lines
PRS_LIST_DIR="${PRS_DIR}"                # contains *.pheno files

# Export R-accessible env vars
export PHENO_DIR="${PHENO_ORIGINAL_DIR}"
export COVAR_FILE="${COVAR_FILE}"
export PRS_DIR="${PRS_DIR}"
export KEEP_FILE="${ANCESTRY_DIR}/whitebrit.fid_iid.txt"
export META_BOTH_FILE="${META_BOTH_FILE}"
export OUT_BASE="${RESULTS_ROOT}/qr_prs_tests"
export BIRTH_YEAR_BASELINE="${BIRTH_YEAR_BASELINE}"

# Which lambda to use (mean|median). Default median; override via --export=MODE=median
MODE="${MODE:-median}"

# Enumerate PRS basenames (sorted, without .pheno)
readarray -t PRS_IDS < <(ls -1 "$PRS_LIST_DIR"/*.pheno 2>/dev/null \
                        | xargs -n1 basename \
                        | LC_ALL=C sort \
                        | sed 's/\.pheno$//')

NUM_PRS=${#PRS_IDS[@]}
if [[ $NUM_PRS -eq 0 ]]; then
  echo "No PRS .pheno files found in $PRS_LIST_DIR" >&2
  exit 2
fi
if [[ $NUM_PRS -ne 36 ]]; then
  echo "Warning: expected 36 PRS, found $NUM_PRS in $PRS_LIST_DIR" >&2
fi

A=${SLURM_ARRAY_TASK_ID}
PHENO_IDX=$(( (A - 1) / NUM_PRS + 1 ))   # 1..25
PRS_IDX=$((   (A - 1) % NUM_PRS + 1 ))   # 1..36

PHENO=$(sed -n "${PHENO_IDX}p" "$PHENO_LIST" | tr -d '\r' || true)
PRS_ID=${PRS_IDS[$((PRS_IDX-1))]:-}

if [[ -z "${PHENO:-}" ]]; then
  echo "No phenotype for index ${PHENO_IDX} (array ${A})" >&2; exit 1
fi
if [[ -z "${PRS_ID:-}" ]]; then
  echo "No PRS for index ${PRS_IDX} (array ${A})" >&2; exit 1
fi

echo "[$(date)] PHENO=${PHENO}  PRS_ID=${PRS_ID}  MODE=${MODE}  (A=${A}, i=${PHENO_IDX}, j=${PRS_IDX})"

Rscript "$RS" "$PHENO" "$PRS_ID" "$MODE"
