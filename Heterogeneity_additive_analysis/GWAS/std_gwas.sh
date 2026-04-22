#!/usr/bin/env bash

#SBATCH -J std_gwas
#SBATCH --mem=10GB
#SBATCH --time=24:00:00
#SBATCH --partition=tier2q
#SBATCH --array=1-22
#SBATCH -o logs/std_gwas_%a_%A.out
#SBATCH -e logs/std_gwas_%a_%A.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

now=$(date)
start_time=$(date +%s)

echo "Date run: $now"

pheno=${1:? "Usage: sbatch std_gwas.sh <pheno_name>"}
echo "Phenotype argument: $pheno"

GENOS="${GENO_DIR}/"
PHENOS="${PHENO_ORIGINAL_DIR}/"
GWAS="${GWAS_RESULTS}/gwas_results/"
GWAS_LOG="${GWAS_RESULTS}/gwas_results_log/"
LOCAL_PLINK2="${PLINK2}"
INTERMEDIATE_DIR="${GWAS_RESULTS}/intermediate_files"

chr=${SLURM_ARRAY_TASK_ID:-21}
pop_path="${GENOS}${TRAIN_POP}/"
pop_file="ukb_chr${chr}_${TRAIN_POP}_QC"
keep_file="${GENO_DIR}/${TRAIN_POP}/ukb_chr1-22_${TRAIN_POP}_QC_train"

META_TAB="${META_BOTH_FILE}"

echo "Processing pheno ${pheno} and chr ${chr}"

PROCESS_PHENO_NAME() {
    local pheno_name=$1
    local phenoNoDigits=${pheno_name%%[0-9]*}
    local phenoLower
    phenoLower=$(echo "$phenoNoDigits" | tr '[:upper:]' '[:lower:]')

    case $pheno_name in
        *"FEV1674178"*) phenoNoDigits="FEV1"; phenoLower="fev1" ;;
        *"IGF-1674178"*) phenoNoDigits="IGF-1"; phenoLower="igf-1" ;;
        *"HbA1c674178"*) phenoNoDigits="HbA1c"; phenoLower="hba1c" ;;
        *"EA4"*) phenoNoDigits="EA4"; phenoLower="ea4" ;;
        *"A1c674178"*) phenoNoDigits="A1c"; phenoLower="a1c" ;;
    esac

    echo "$phenoNoDigits $phenoLower"
}

parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
phenoNoDigits=${parse_pheno%% *}
phenoLower=${parse_pheno##* }

if [[ "$phenoNoDigits" == "FEV" ]]; then
    phenoNoDigits="FEV1"
    phenoLower="fev1"
fi

echo "Processed phenotype name: phenoNoDigits=${phenoNoDigits}, phenoLower=${phenoLower}"

lambda=$(awk -v p="$phenoNoDigits" 'NR>1 && $1==p {print $11; exit}' "$META_TAB")
if [[ -z "$lambda" ]]; then
    echo "ERROR: Could not find lambda for phenotype '${phenoNoDigits}' in ${META_TAB}" >&2
    exit 1
fi
echo "Lambda: $lambda"

pheno_file="${PHENOS}${phenoNoDigits}/${pheno}.pheno"
echo "Phenotype file: $pheno_file"

pheno_code=$(head -1 "$pheno_file" | awk '{ print $3 }')
echo "Phenotype code: $pheno_code"

echo "Applying Box-Cox transform with lambda = $lambda"
mkdir -p "$INTERMEDIATE_DIR"
pheno_file_bc="${INTERMEDIATE_DIR}/${pheno}_bc.pheno"

lock="${pheno_file_bc}.lock"
exec 9>"$lock"
flock -x 9

if [[ ! -s "$pheno_file_bc" ]]; then
    tmp=$(mktemp "${pheno_file_bc}.tmp.XXXXXX")

    awk -v lam="$lambda" '
    BEGIN { OFS="\t" }
    NR == 1 { print $0; next }  # Keep header row
    {
        # Require numeric, positive phenotype
        if (($3+0)==$3 && $3 > 0) {
            y = $3
            if (lam == 0) {
                $3 = log(y)
            } else {
                $3 = (exp(lam * log(y)) - 1.0) / lam
            }
        } else {
            $3 = "NA"
        }
        print $0
    }' "$pheno_file" > "$tmp"

    # atomic replace (prevents PLINK from ever seeing a half-written file)
    mv "$tmp" "$pheno_file_bc"
fi

flock -u 9
exec 9>&-

pheno_file="$pheno_file_bc"
echo "Transformed phenotype file saved to: $pheno_file"

if (( $(echo "$lambda == 0" | bc -l) )); then
    gwas_path="${GWAS_LOG}${phenoLower}/"
else
    gwas_path="${GWAS}${phenoLower}/"
fi
mkdir -p "$gwas_path"

covar_file="${COVAR_FILE}"
ncovar='3-14'

out_file="${gwas_path}chr${chr}_${TRAIN_POP}_${phenoLower}"

if [[ ! -f "$out_file.$pheno_code.glm.linear" ]]; then
    ${LOCAL_PLINK2} --bfile ${pop_path}${pop_file} \
        --keep ${keep_file} \
        --pheno "$pheno_file" \
        --pheno-name "$pheno_code" \
        --no-input-missing-phenotype \
        --glm hide-covar \
        --covar-variance-standardize 34-0.0 \
        --covar "$covar_file" \
        --covar-col-nums "$ncovar" \
        --out "$out_file"
else
    echo "Output file already exists: $out_file.$pheno_code.glm.linear"
fi

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds"
