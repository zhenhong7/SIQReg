#!/usr/bin/env bash
#SBATCH --job-name=hetero_gamma_power
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --array=1-11000
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail
module load gcc/12.1.0
module load R/4.4.1

OUTDIR="./power"
LOGDIR="./logs"
mkdir -p "$OUTDIR" "$LOGDIR"

OUTFILE="${OUTDIR}/prs_hetero_power.csv"
LOCKFILE="${OUTDIR}/prs_hetero_power.lock"

RSCRIPT="./hete_power.R"
[[ -f "$RSCRIPT" ]] || { echo "FATAL: missing RSCRIPT: $RSCRIPT" >&2; exit 2; }

GAMMAS=(0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0)

TASK_ID=${SLURM_ARRAY_TASK_ID}
rep=$(( (TASK_ID - 1) % 1000 + 1 ))
g_idx=$(( (TASK_ID - 1) / 1000 ))
gamma=${GAMMAS[$g_idx]}

echo "Task ${TASK_ID}: gamma=${gamma}, rep=${rep}"

tmpfile=$(mktemp)
Rscript "$RSCRIPT" "$gamma" "$rep" > "$tmpfile"

# append atomically with a lock; write header once
exec 9>"$LOCKFILE"
flock -x 9
if [ ! -f "$OUTFILE" ]; then
  echo "n,lambda0,gamma,rep,lambda_hat,D_value,eq_p_y,eq_p_log_y,eq_p_rint_y,eq_p_y_star,eq_p_y_star_est" > "$OUTFILE"
fi
cat "$tmpfile" >> "$OUTFILE"
flock -u 9
exec 9>&-

rm -f "$tmpfile"
