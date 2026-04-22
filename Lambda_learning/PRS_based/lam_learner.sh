#!/usr/bin/env bash
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-900
#SBATCH --job-name=lamPRS
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

source ../../config.sh

module load gcc/12.1.0
module load R/4.4.1

# Pull pheno from the first argument
pheno=${1:?Usage: $0 PHENO}
OUTDIR="${RESULTS_DIR:-results}/${pheno}"

# (A) On login node: init + re-sbatch as array
if [ -z "$SLURM_JOB_ID" ]; then
  mkdir -p "${OUTDIR}/logs"
  for f in "${PRS_DIR}"/*.pheno; do
    prs=$(basename "$f" .pheno)
    printf "prs_id,chunk_i,lambda_hat,time_sec\n" \
      > "${OUTDIR}/lambda_${prs}.csv"
  done

  # Re-submit this script as an array
  sbatch "$0" "${pheno}"
  exit 0
fi

# (B) Within array tasks: compute indices
task_id=$SLURM_ARRAY_TASK_ID
n_chunks=25
prs_idx=$(( (task_id - 1) / n_chunks ))
chunk_i=$(( ((task_id - 1) % n_chunks) + 1 ))
PRS_FILES=("${PRS_DIR}"/*.pheno)
THIS_PRS=$(basename "${PRS_FILES[$prs_idx]}" .pheno)

mkdir -p "${OUTDIR}"

# (C) Run R and append
Rscript ./lam_learner.R \
        "${pheno}" "${THIS_PRS}" "${chunk_i}" \
  >> "${OUTDIR}/lambda_${THIS_PRS}.csv"
