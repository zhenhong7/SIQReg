#!/usr/bin/env bash
#SBATCH -J join_gwas
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --partition=tier2q
#SBATCH -o logs/join_gwas_%A.out
#SBATCH -e logs/join_gwas_%A.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

echo "Date run: $(date)"

pheno=${1:? "Usage: sbatch join_gwas.sh <phenoLower>  (e.g. bmi)"}

GWAS="${GWAS_RESULTS}/gwas_results"
JOINED_GWAS="${GWAS_RESULTS}/joined_gwas_results"

input_dir="${GWAS}/${pheno}"
output_dir="${JOINED_GWAS}/${pheno}"
mkdir -p "$output_dir"

# Find chr1 file (glm.linear) to get the exact suffix like ".21001.glm.linear"
chr1_file=$(ls -1 "${input_dir}/chr1_${TRAIN_POP}_${pheno}."*.glm.linear 2>/dev/null | head -n 1 || true)
if [[ -z "$chr1_file" ]]; then
  echo "FATAL: cannot find chr1 GWAS file in ${input_dir} for pheno=${pheno}" >&2
  echo "Expected something like: ${input_dir}/chr1_${TRAIN_POP}_${pheno}.<code>.glm.linear" >&2
  exit 1
fi

suffix="${chr1_file#${input_dir}/chr1_${TRAIN_POP}_${pheno}.}"   # e.g. "21001.glm.linear"
output_file="${output_dir}/chr1-22_${TRAIN_POP}_${pheno}.${suffix}"

tmp=$(mktemp "${output_file}.tmp.XXXXXX")

echo "Merging into: $output_file"
echo "Using suffix:  $suffix"

# header from chr1
head -n 1 "$chr1_file" > "$tmp"

# append chr1..22 data (no headers)
tail -n +2 -q "${input_dir}"/chr{1..22}_"${TRAIN_POP}"_"${pheno}"."${suffix}" >> "$tmp"

mv "$tmp" "$output_file"
echo "Done."
