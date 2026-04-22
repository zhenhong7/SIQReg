#!/usr/bin/env bash
#SBATCH --job-name=lam_effsize_submit
#SBATCH --partition=tier1q
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -o logs/lam_effsize_submit_%A.out
#SBATCH -e logs/lam_effsize_submit_%A.err

set -euo pipefail
module load gcc/12.1.0
module load R/4.4.1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPS=1000
B=300
ATS=(0.02 0.05 0.10 0.15 0.20 0.30 0.40 0.50)

THROTTLE=500

RSCRIPT="./VaryingEffectSize_t.R"
OUTDIR="./results_t_effsize"
LOGDIR="${OUTDIR}/logs"
RAWDIR="${OUTDIR}/raw"

mkdir -p "$OUTDIR" "$LOGDIR" "$RAWDIR"

for at in "${ATS[@]}"; do
  f="${RAWDIR}/raw_aT${at}.csv"
  if [[ ! -s "$f" ]]; then
    echo "rep,lambda_hat_D,lambda_hat_ML" > "$f"
  fi
done

NAT=${#ATS[@]}
TOTAL=$((NAT * REPS))   # 8 * 1000 = 8000

# array
ARRAY_JOBID=$(sbatch --parsable \
  --job-name=lam1rep_eff \
  --partition=tier1q \
  --time=24:00:00 \
  --mem=1G \
  --cpus-per-task=1 \
  -o "${LOGDIR}/lam1rep_eff_%A.out" \
  -e "${LOGDIR}/lam1rep_eff_%A.err" \
  --array=1-${TOTAL}%${THROTTLE} \
  --wrap "
set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPS=${REPS}
B=${B}
RSCRIPT='${RSCRIPT}'
RAWDIR='${RAWDIR}'
LOGDIR='${LOGDIR}'

ATS=(${ATS[*]})

TASK=\${SLURM_ARRAY_TASK_ID}

rep=\$(( (TASK - 1) % REPS + 1 ))
at_idx=\$(( (TASK - 1) / REPS ))

at=\${ATS[\$at_idx]}

rawfile=\"\${RAWDIR}/raw_aT\${at}.csv\"
lockfile=\"\${rawfile}.lock\"
rerr=\"\${LOGDIR}/R_\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.err\"

Rscript -e 'suppressPackageStartupMessages(library(conquer))' 2> \"\$rerr\"

# args: aT rep_id B
line=\$(Rscript \"\$RSCRIPT\" \"\$at\" \"\$rep\" \"\$B\" 2>> \"\$rerr\")

(
  flock -x 200
  echo \"\$line\" >> \"\$rawfile\"
) 200>\"\$lockfile\"
")

# aggregation
AGG_JOBID=$(sbatch --parsable \
  --dependency=afterok:${ARRAY_JOBID} \
  --job-name=lam_agg_eff \
  --partition=tier1q \
  --time=02:00:00 \
  --mem=4G \
  --cpus-per-task=1 \
  -o "${LOGDIR}/lam_agg_eff_%A.out" \
  -e "${LOGDIR}/lam_agg_eff_%A.err" \
  --wrap "Rscript -e '
outdir <- \"${OUTDIR}\"
rawdir <- \"${RAWDIR}\"
ats <- c(${ATS[*]})

sum_rows <- list()

for (at in ats) {
  at_str <- format(at, nsmall=2, trim=TRUE)
  rawf <- file.path(rawdir, sprintf(\"raw_aT%s.csv\", at_str))
  if (!file.exists(rawf)) next
  dat <- read.csv(rawf)

  mean_hat_D  <- mean(dat\$lambda_hat_D,  na.rm=TRUE)
  se_hat_D    <- sd(dat\$lambda_hat_D,    na.rm=TRUE)
  mean_hat_ML <- mean(dat\$lambda_hat_ML, na.rm=TRUE)
  se_hat_ML   <- sd(dat\$lambda_hat_ML,   na.rm=TRUE)
  kept <- sum(is.finite(dat\$lambda_hat_D) & is.finite(dat\$lambda_hat_ML))

  df_full <- data.frame(
    aT=at, rep=dat\$rep,
    lambda_hat_D=dat\$lambda_hat_D, lambda_hat_ML=dat\$lambda_hat_ML,
    mean_hat_D=mean_hat_D, se_hat_D=se_hat_D,
    mean_hat_ML=mean_hat_ML, se_hat_ML=se_hat_ML
  )
  outf <- file.path(outdir, sprintf(\"sim_lambda_aT%s.csv\", at_str))
  write.csv(df_full, outf, row.names=FALSE)

  sum_rows[[length(sum_rows)+1L]] <- data.frame(
    aT=at, kept=kept,
    mean_hat_D=mean_hat_D, se_hat_D=se_hat_D,
    mean_hat_ML=mean_hat_ML, se_hat_ML=se_hat_ML
  )
}

df_sum <- do.call(rbind, sum_rows)
df_sum <- df_sum[order(df_sum\$aT), ]
write.csv(df_sum, file.path(outdir, \"sim_lambda_effsize_summary.csv\"), row.names=FALSE)
cat(\"Wrote sim_lambda_effsize_summary.csv\\n\")
'")

echo "Submitted array job: ${ARRAY_JOBID} (tasks=${TOTAL})"
echo "Submitted agg job:   ${AGG_JOBID}"
echo "OUTDIR: ${OUTDIR}"
