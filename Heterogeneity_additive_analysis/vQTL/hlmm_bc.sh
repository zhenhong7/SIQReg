#!/usr/bin/env bash
#SBATCH -J HLMM_BC
#SBATCH --time=240:00:00
#SBATCH --partition=tier3q
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --array=1-22
#SBATCH -o logs/hlmm_%x_%A.out
#SBATCH -e logs/hlmm_%x_%A.err

set -euo pipefail

pheno=${1:? "Usage: sbatch hlmm_bc.sh <PHENO>  (e.g. BMI, Height, LDL)"}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

# 0) Clean module env
module purge
module load gcc/12.1.0

# 1) Conda activation robust (avoid libmamba plugin issues)
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
LOCAL_GENO_DIR="${GENO_DIR}/${TRAIN_POP}"
COV="${HLMM_COV}"
HLMM_PY_BIN="${HLMM_FAST_PY}"

PHENO_DIR="${HLMM_PHENO_BC_DIR}"
PHENO_FILE=$(ls -1 "${PHENO_DIR}/${pheno}"*_bc.pheno 2>/dev/null | head -n 1)
if [ -z "${PHENO_FILE}" ]; then
  echo "ERROR: cannot find phenotype file: ${PHENO_DIR}/${pheno}*_bc.pheno" >&2
  exit 1
fi

OUTDIR="${HLMM_RESULTS}/${pheno}_bc"
mkdir -p "${OUTDIR}"

CHR=${SLURM_ARRAY_TASK_ID}
BFILE=${LOCAL_GENO_DIR}/ukb_chr${CHR}_${TRAIN_POP}_QC
BED="${BFILE}.bed"
BIM="${BFILE}.bim"

M=$(wc -l < "${BIM}")
BLOCK=${BLOCK:-20000}

OUTPREFIX="${OUTDIR}/${pheno}.chr${CHR}"

# 4) Chunked HLMM
start=0
while [ "${start}" -lt "${M}" ]; do
  end=$(( start + BLOCK ))
  if [ "${end}" -gt "${M}" ]; then end="${M}"; fi

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
done

echo "DONE ${pheno} chr${CHR}"
