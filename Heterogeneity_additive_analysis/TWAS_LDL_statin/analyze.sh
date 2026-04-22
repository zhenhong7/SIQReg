#!/bin/bash
#SBATCH --job-name=statin_strat
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=240:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

module load gcc/12.1.0
module load R/4.4.1

mkdir -p logs

# Export R-accessible env vars
export PHENO_DIR="${PHENO_ORIGINAL_DIR}"
export COVAR_FILE="${COVAR_FILE}"
export EXPR_DIR="${EXPR_DIR}"
export KEEP_FILE="${ANCESTRY_DIR}/whitebrit.fid_iid.txt"
export STATIN_FILE="${STATIN_FILE}"
export META_BOTH_FILE="${META_BOTH_FILE}"
export TWAS_RESULTS="${TWAS_RESULTS}"
export BIRTH_YEAR_BASELINE="${BIRTH_YEAR_BASELINE}"

Rscript "${SCRIPT_DIR}/analyze.R"
