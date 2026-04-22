#!/usr/bin/env bash
#SBATCH --job-name=lam25000_submit
#SBATCH --partition=tier1q
#SBATCH --time=24:00:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH -o logs/lam25000_submit_%A.out
#SBATCH -e logs/lam25000_submit_%A.err

set -euo pipefail
module load gcc/12.1.0
module load R/4.4.1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPS=1000
B=300
LAMBDAS=(-1 -0.5 0 0.5 1)
NS=(1000 2000 3000 4000 5000)

# throttle for the array
THROTTLE=500

RSCRIPT="./VaryingSampleSize_t.R"
OUTDIR="./results_t"
LOGDIR="${OUTDIR}/logs"
RAWDIR="${OUTDIR}/raw"

mkdir -p "$OUTDIR" "$LOGDIR" "$RAWDIR"

# init 25 raw files ONLY if missing or empty (avoid wiping on rerun)
for lam in "${LAMBDAS[@]}"; do
  for n in "${NS[@]}"; do
    f="${RAWDIR}/raw_n${n}_lam${lam}.csv"
    if [[ ! -s "$f" ]]; then
      echo "rep,lambda_hat_D,lambda_hat_ML" > "$f"
    fi
  done
done

NLAM=${#LAMBDAS[@]}
NN=${#NS[@]}
TOTAL=$((NLAM * NN * REPS))   # 25000

# array: each task runs one rep and appends to its combo raw file
ARRAY_JOBID=$(sbatch --parsable \
  --job-name=lam1rep \
  --partition=tier1q \
  --time=24:00:00 \
  --mem=1G \
  --cpus-per-task=1 \
  -o "${LOGDIR}/lam1rep_%A.out" \
  -e "${LOGDIR}/lam1rep_%A.err" \
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

LAMBDAS=(${LAMBDAS[*]})
NS=(${NS[*]})

NLAM=\${#LAMBDAS[@]}
TASK=\${SLURM_ARRAY_TASK_ID}

rep=\$(( (TASK - 1) % REPS + 1 ))
combo=\$(( (TASK - 1) / REPS ))

lam_idx=\$(( combo % NLAM ))
n_idx=\$(( combo / NLAM ))

lam=\${LAMBDAS[\$lam_idx]}
n=\${NS[\$n_idx]}

rawfile=\"\${RAWDIR}/raw_n\${n}_lam\${lam}.csv\"
lockfile=\"\${rawfile}.lock\"
rerr=\"\${LOGDIR}/R_\${SLURM_ARRAY_JOB_ID}_\${SLURM_ARRAY_TASK_ID}.err\"

Rscript -e 'suppressPackageStartupMessages(library(conquer))' 2> \"\$rerr\"

line=\$(Rscript \"\$RSCRIPT\" \"\$lam\" \"\$n\" \"\$rep\" \"\$B\" 2>> \"\$rerr\")

(
  flock -x 200
  echo \"\$line\" >> \"\$rawfile\"
) 200>\"\$lockfile\"
")

# aggregation job: after array, create final 25 combo outputs + one summary
AGG_JOBID=$(sbatch --parsable \
  --dependency=afterok:${ARRAY_JOBID} \
  --job-name=lam_agg \
  --partition=tier1q \
  --time=02:00:00 \
  --mem=4G \
  --cpus-per-task=1 \
  -o "${LOGDIR}/lam_agg_%A.out" \
  -e "${LOGDIR}/lam_agg_%A.err" \
  --wrap "Rscript -e '
outdir <- \"${OUTDIR}\"
rawdir <- \"${RAWDIR}\"
lambdas <- c(${LAMBDAS[*]})
ns <- c(${NS[*]})

sum_rows <- list()

for (lam0 in lambdas) for (nn in ns) {
  rawf <- file.path(rawdir, sprintf(\"raw_n%d_lam%s.csv\", nn, format(lam0, trim=TRUE)))
  if (!file.exists(rawf)) next
  dat <- read.csv(rawf)

  mean_hat_D  <- mean(dat\$lambda_hat_D,  na.rm=TRUE)
  se_hat_D    <- sd(dat\$lambda_hat_D,    na.rm=TRUE)
  mean_hat_ML <- mean(dat\$lambda_hat_ML, na.rm=TRUE)
  se_hat_ML   <- sd(dat\$lambda_hat_ML,   na.rm=TRUE)
  kept <- sum(is.finite(dat\$lambda_hat_D) & is.finite(dat\$lambda_hat_ML))

  df_full <- data.frame(
    n=nn,
    lambda0=lam0,
    rep=dat\$rep,
    lambda_hat_D=dat\$lambda_hat_D,
    lambda_hat_ML=dat\$lambda_hat_ML,
    mean_hat_D=mean_hat_D,
    se_hat_D=se_hat_D,
    mean_hat_ML=mean_hat_ML,
    se_hat_ML=se_hat_ML
  )

  outf <- file.path(outdir, sprintf(\"sim_lambda_n%d_lam%s.csv\", nn, format(lam0, trim=TRUE)))
  write.csv(df_full, outf, row.names=FALSE)

  sum_rows[[length(sum_rows)+1L]] <- data.frame(
    n=nn, lambda0=lam0, kept=kept,
    mean_hat_D=mean_hat_D, se_hat_D=se_hat_D,
    mean_hat_ML=mean_hat_ML, se_hat_ML=se_hat_ML
  )
}

df_sum <- do.call(rbind, sum_rows)
df_sum <- df_sum[order(df_sum\$n, df_sum\$lambda0), ]
write.csv(df_sum, file.path(outdir, \"sim_lambda_summary.csv\"), row.names=FALSE)
cat(\"Wrote sim_lambda_summary.csv\\n\")
'")

echo "Submitted array job: ${ARRAY_JOBID} (tasks=${TOTAL})"
echo "Submitted agg job:   ${AGG_JOBID}"
echo "OUTDIR: ${OUTDIR}"
