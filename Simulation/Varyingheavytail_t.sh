#!/usr/bin/env bash
#SBATCH --job-name=lambda_tdf
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --mem=1G
#SBATCH --array=1-11
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

mkdir -p logs results_varyingt
module load gcc/12.1.0
module load R/4.4.1

TDFS=(1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5)

task=${SLURM_ARRAY_TASK_ID}
idx=$((task - 1))
t_df=${TDFS[$idx]}

echo "Task ${task}: t_df=${t_df}"

Rscript Varyingheavytail_t.R "${t_df}" 1000 300 "./results_varyingt" 0.5
