#!/usr/bin/env bash
#SBATCH --job-name=power_heavytail
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --array=1-11000
#SBATCH --output=logs/power_heavytail_%A.out
#SBATCH --error=logs/power_heavytail_%A.err

set -euo pipefail
module load gcc/12.1.0
module load R/4.4.1

mkdir -p logs

TASK_ID=${SLURM_ARRAY_TASK_ID}

N_BG=11           # number of bG grid points
REPS_PER_BG=1000  # total reps per bG

# Map TASK_ID (1..11000) -> (bG_index 1..11, rep_index 1..1000)
bG_index=$(( (TASK_ID - 1) % N_BG + 1 ))
rep_index=$(( (TASK_ID - 1) / N_BG + 1 ))

# safety check (should always be TRUE for array=1-11000)
if [ "$rep_index" -gt "$REPS_PER_BG" ]; then
  echo "rep_index ${rep_index} > ${REPS_PER_BG}, exiting."
  exit 0
fi

echo "TASK_ID=${TASK_ID}, bG_index=${bG_index}, rep_index=${rep_index}"

Rscript power_heavytail.R "${bG_index}" "${rep_index}"
