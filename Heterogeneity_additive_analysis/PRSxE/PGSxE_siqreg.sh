#!/bin/bash
#SBATCH --job-name=PGSxE_SI
#SBATCH --array=1-900
#SBATCH --time=240:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/slurm_%A.out
#SBATCH --error=logs/slurm_%A.err

ENV=$1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

module load gcc/12.1.0
module load R/4.4.1

mkdir -p logs

# Export R-accessible env vars
export PHENOS_FILE="${PHENO_DIR}/phenos.txt"
export PRS_SCORING_DIR="${PRS_SCORING_DIR}"
export PHENO_ORIGINAL_DIR="${PHENO_ORIGINAL_DIR}"
export COVAR_FULL_FILE="${COVAR_FULL_FILE}"
export KEEP_FILE="${ANCESTRY_DIR}/whitebrit.fid_iid.txt"
export GXE_META_FILE="${GXE_META_FILE}"
export OUTPUT_DIR="${GXE_DIR}"
export ERROR_LOG_FILE="logs/error_log.txt"
export ID_LOG_FILE="logs/id_pheno_prs.txt"

Rscript "${SCRIPT_DIR}/PGSxE_siqreg.R" ${SLURM_ARRAY_TASK_ID} ${ENV}


# Submit for all 7 environments
#for env in sex age EA4 BMI Smoking_status Alcohol_intake_frequency Statins; do
#    sbatch PGSxE_siqreg.sh $env
#done
