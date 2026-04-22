#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config.sh"

PHENOS="${PHENO_ORIGINAL_DIR}/"
META_TAB="${META_BOTH_FILE}"

GWAS_SCRIPT="${SCRIPT_DIR}/std_gwas.sh"

while read phenotype k_DL lambda_DL rest; do
    if [[ "$phenotype" == "phenotype" ]]; then
        continue
    fi

    pheno_dir="${PHENOS}${phenotype}"

    if [[ ! -d "$pheno_dir" ]]; then
        echo "WARNING: phenotype directory not found: $pheno_dir" >&2
        continue
    fi

    pheno_file=$(ls "${pheno_dir}"/*.pheno 2>/dev/null | head -n 1)

    if [[ -z "$pheno_file" ]]; then
        echo "WARNING: no .pheno file found in $pheno_dir" >&2
        continue
    fi

    pheno_id=$(basename "$pheno_file" .pheno)

    echo "Submitting GWAS for phenotype=${phenotype}, pheno_id=${pheno_id}"
    sbatch "$GWAS_SCRIPT" "$pheno_id"

done < "$META_TAB"
