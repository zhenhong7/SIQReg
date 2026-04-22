#!/usr/bin/env bash
#SBATCH -J comp_valid_pgs
#SBATCH --partition=tier2q
#SBATCH --time=240:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH -o logs/compile_valid_%A_%a.out
#SBATCH -e logs/compile_valid_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail
module load gcc/12.1.0
module load R/4.4.1

PHENOLIST="${PHENOLIST:-/path/to/SI_GWAS/phenolist_27.txt}"
RSCRIPT_PATH="${RSCRIPT_PATH:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/compile_valid_pgs.R}"

PHENO_DIGITS=$(awk 'NF{print}' "$PHENOLIST" | sed -n "${SLURM_ARRAY_TASK_ID}p")

echo "Date run: $(date)"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "PHENO_DIGITS: ${PHENO_DIGITS}"

if [[ -z "${PHENO_DIGITS:-}" ]]; then
  echo "FATAL: no phenotype found for task id ${SLURM_ARRAY_TASK_ID} in ${PHENOLIST}" >&2
  exit 1
fi

Rscript "$RSCRIPT_PATH" "$PHENO_DIGITS"
