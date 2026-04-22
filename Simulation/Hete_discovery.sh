#!/usr/bin/env bash
#SBATCH --job-name=hetero_qr
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --array=1-5000
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail
module load gcc/12.1.0
module load R/4.4.1

OUTDIR="./results_hetero"
LOGDIR="./logs"
mkdir -p "$OUTDIR" "$OUTDIR/rows" "$LOGDIR"

LAMBDA_VALS=(-1 -0.5 0 0.5 1)

TASK_ID=${SLURM_ARRAY_TASK_ID}

# rep_id: 1..1000
rep=$(( (TASK_ID - 1) % 1000 + 1 ))

# lambda index: 0..4 then +1
lambda_idx=$(( (TASK_ID - 1) / 1000 ))
lambda0=${LAMBDA_VALS[$lambda_idx]}

Rscript sim_prs_hetero_job.R "$lambda0" "$rep" "$OUTDIR"
