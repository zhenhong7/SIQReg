#!/usr/bin/env bash
#SBATCH --partition=tier1q
#SBATCH --time=240:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-215
#SBATCH --job-name=lamcov
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.err

set -euo pipefail

source ../../config.sh

module load gcc/12.1.0
module load R/4.4.1

PHENO=${1:?Usage: $0 PHENO_ID [ANCESTRY]}

RSCRIPT="./lambda_covariate.R"

RAW_ANCESTRY="${2:-whitebrit}"
ANCESTRY="$(printf "%s" "$RAW_ANCESTRY" | tr -d '\r' | xargs)"
export ANCESTRY

# chunks by ancestry
case "$ANCESTRY" in
  whitebrit)  N_CHUNKS=25 ;;
  white_euro) N_CHUNKS=20 ;;
  afr)        N_CHUNKS=5  ;;
  asn)        N_CHUNKS=5  ;;
  *) echo "Unknown ANCESTRY=<$ANCESTRY> (raw=<$RAW_ANCESTRY>)" >&2; exit 1 ;;
esac
export N_CHUNKS
echo "ANCESTRY=<${ANCESTRY}>  N_CHUNKS=${N_CHUNKS}"

# Results root
OUTROOT="${OUTROOT:-results_covariates}"
export OUTROOT
OUTDIR="${OUTROOT}/${ANCESTRY}/${PHENO}"
mkdir -p "${OUTDIR}"

# Covariate list (43): sex, age, age2, PC1..PC40
COVARS=(sex age age2 $(printf "PC%d " {1..40}))

chunks=${N_CHUNKS}
task_id=${SLURM_ARRAY_TASK_ID}
cov_idx=$(( (task_id - 1) / chunks + 1 ))   # 1..43
chunk_i=$(( (task_id - 1) % chunks + 1 ))   # 1..N_CHUNKS

if (( cov_idx < 1 || cov_idx > ${#COVARS[@]} )); then
  echo "cov_idx out of range: ${cov_idx}" >&2
  exit 1
fi
cov="${COVARS[$((cov_idx-1))]}"

BASE_TMP="${SLURM_TMPDIR:-/scratch/$USER/tmp}"
mkdir -p "$BASE_TMP" || true
TMPDIR="$(mktemp -d -p "$BASE_TMP" "Rtmp_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}_XXXXXXXX")" || {
  mkdir -p "$HOME/tmp"
  TMPDIR="$(mktemp -d -p "$HOME/tmp" "Rtmp_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}_XXXXXXXX")"
}
export TMPDIR TMP="$TMPDIR" TEMP="$TMPDIR"
chmod 700 "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT
echo "Node: $(hostname)"
echo "TMPDIR=$TMPDIR"
df -h "$TMPDIR" || true

echo "PHENO=${PHENO}  COV=${cov}  CHUNK=${chunk_i}"
Rscript "$RSCRIPT" "$PHENO" "$chunk_i" "$cov"
