#!/usr/bin/env bash
#SBATCH -J HLMM
#SBATCH --time=240:00:00
#SBATCH --partition=tier3q
#SBATCH --cpus-per-task=1
#SBATCH --mem=55G
#SBATCH --array=1-22
#SBATCH -o logs/hlmm_%x_chr%a_%A.out
#SBATCH -e logs/hlmm_%x_chr%a_%A.err

set -euo pipefail

PHENO_NAME=${1:? "Usage: sbatch hlmm_original.sh <PHENO_DIR_NAME>   (e.g. LDL, Height, BMI)"}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# 0) Clean module env
module purge
module load gcc/12.1.0

# 1) Make conda activation robust (avoid libmamba plugin issues)
unset LD_PRELOAD
unset LD_LIBRARY_PATH
export CONDA_NO_PLUGINS=true
export CONDA_SOLVER=classic

source "$HOME/miniconda/etc/profile.d/conda.sh"
conda activate hlmm_py27

# 2) Fix pysnptools .so runtime: use GCC12 libstdc++/libgcc (NOT conda's)
LIBSTDCXX="$(g++ --print-file-name=libstdc++.so.6)"
LIBGCC="$(gcc --print-file-name=libgcc_s.so.1)"
export LD_PRELOAD="${LIBSTDCXX}:${LIBGCC}"
export LD_LIBRARY_PATH="$(dirname "$LIBSTDCXX"):${LD_LIBRARY_PATH:-}"

# 3) Paths
GENO_DIR="${GENO_DIR}/${TRAIN_POP}"
PHENO_DIR="${HLMM_PHENO_DIR}"
PHENO_FILE=$(ls -1 "${PHENO_DIR}/${PHENO_NAME}"*.pheno | head -n 1)
COV="${HLMM_COV}"
HLMM_PY_BIN="${HLMM_PY}"

OUTDIR="${HLMM_RESULTS}/${PHENO_NAME}"
mkdir -p "${OUTDIR}"

CHR=${SLURM_ARRAY_TASK_ID}
BFILE=${GENO_DIR}/ukb_chr${CHR}_${TRAIN_POP}_QC
BED="${BFILE}.bed"
BIM="${BFILE}.bim"

# 4) Chunked HLMM to avoid MemoryError
M=$(wc -l < "${BIM}")
echo "PHENO=${PHENO_NAME}  CHR=${CHR}  M=${M} SNPs"
echo "PHENO_FILE=${PHENO_FILE}"

BLOCK=${BLOCK:-10000}   # override: sbatch --export=BLOCK=20000 hlmm_original.sh LDL
echo "Using BLOCK=${BLOCK} SNPs per chunk"

OUTPREFIX="${OUTDIR}/${PHENO_NAME}.chr${CHR}"

start=0
chunk=0
while [ "${start}" -lt "${M}" ]; do
  end=$(( start + BLOCK ))
  if [ "${end}" -gt "${M}" ]; then end="${M}"; fi

  echo "Running chunk ${chunk}: SNP [${start}, ${end})"

  if [ "${start}" -eq 0 ]; then
    python "${HLMM_PY_BIN}" "${BED}" "${start}" "${end}" "${PHENO_FILE}" "${OUTPREFIX}" \
      --phen_index 1 \
      --mean_covar "${COV}" \
      --var_covar  "${COV}"
  else
    python "${HLMM_PY_BIN}" "${BED}" "${start}" "${end}" "${PHENO_FILE}" "${OUTPREFIX}" \
      --phen_index 1 \
      --mean_covar "${COV}" \
      --var_covar  "${COV}" \
      --append
  fi

  start="${end}"
  chunk=$((chunk + 1))
done

echo "DONE ${PHENO_NAME} chr${CHR}"
