#!/usr/bin/env bash
#SBATCH -J PI_pgs
#SBATCH --partition=tier3q
#SBATCH --time=24:00:00
#SBATCH --mem=75GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/pgs_%A_%a.out
#SBATCH -e logs/pgs_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SCALE=${1:? "Usage: sbatch merge_clump_pgs.sh <scale>"}

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
TRAIT_STR="${TRAIT_INFO[$TASK_ID]}"
IFS=':' read -r phenoLower phenoFull phenoCode <<< "${TRAIT_STR}"

if [[ "${SCALE}" == "original" ]]; then
    OUTDIR="${OUTBASE}/original_scale/${phenoLower}"
else
    OUTDIR="${OUTBASE}/bc_scale/${phenoLower}"
fi

echo "${phenoLower} (${SCALE})"

GWAS_DIR="${OUTDIR}/gwas"
PGS_DIR="${OUTDIR}/pgs"
mkdir -p "${PGS_DIR}"

module load gcc/12.1.0
module load R/4.4.1

CV_VALID_PREFIX="${OUTDIR}/folds/cv_valid_fold"
TEST_SAMPLES="${OUTDIR}/test_samples.txt"

if [[ ! -f "${TEST_SAMPLES}" ]]; then
    echo "ERROR: Test samples file not found: ${TEST_SAMPLES}"
    exit 1
fi

for fold in $(seq 1 ${K}); do
    echo "fold ${fold}"

    MERGED="${GWAS_DIR}/gwas_fold_${fold}.txt"
    CHR1="${GWAS_DIR}/chr1_fold${fold}.PHENO.glm.linear"

    if [[ ! -f "${CHR1}" ]]; then
        echo "ERROR: Chr1 GWAS not found: ${CHR1}"
        continue
    fi

    head -n 1 "${CHR1}" > "${MERGED}"
    for chr in $(seq 1 22); do
        CHR_FILE="${GWAS_DIR}/chr${chr}_fold${fold}.PHENO.glm.linear"
        if [[ -f "${CHR_FILE}" ]]; then
            tail -n +2 "${CHR_FILE}" >> "${MERGED}"
        fi
    done
    CLUMP_OUT="${GWAS_DIR}/clumped_fold_${fold}"
    ${PLINK2} --bfile "${GENO_MERGED}" --clump "${MERGED}" --clump-p1 1 \
        --clump-r2 ${CLUMP_R2} --clump-kb ${CLUMP_KB} --out "${CLUMP_OUT}" \
        > "${OUTDIR}/logs/clump_fold_${fold}.log" 2>&1 || true

    CLUMP_FILE="${CLUMP_OUT}.clumps"
    if [[ -f "${CLUMP_FILE}" ]]; then
        awk 'NR>1 && NF>0 {print $3}' "${CLUMP_FILE}" > "${GWAS_DIR}/clumped_snps_${fold}.txt"
    else
        awk 'NR>1 {print $3}' "${MERGED}" > "${GWAS_DIR}/clumped_snps_${fold}.txt"
    fi

    SCORE_FILE="${GWAS_DIR}/score_fold_${fold}.txt"

    Rscript --vanilla - <<RSCRIPT
library(data.table)
gwas <- fread("${MERGED}")

id_col <- intersect(c("ID", "SNP"), colnames(gwas))[1]
a1_col <- intersect(c("A1", "ALT"), colnames(gwas))[1]
beta_col <- intersect(c("BETA", "B"), colnames(gwas))[1]
p_col <- intersect(c("P", "P_VALUE"), colnames(gwas))[1]

setnames(gwas, c(id_col, a1_col, beta_col, p_col),
         c("SNP", "A1", "BETA", "P"), skip_absent=TRUE)

clumped <- fread("${GWAS_DIR}/clumped_snps_${fold}.txt", header=FALSE)[["V1"]]
gwas <- gwas[SNP %in% clumped]

gwas <- gwas[!is.na(P) & P < ${P_THRESHOLD}]

cat("SNPs passing P <", ${P_THRESHOLD}, ":", nrow(gwas), "\n")

score <- gwas[, .(SNP, A1, BETA)]
fwrite(score, "${SCORE_FILE}", sep="\t", col.names=TRUE)
RSCRIPT

    CV_VALID_FILE="${CV_VALID_PREFIX}_${fold}.txt"
    ${PLINK2} --bfile "${GENO_MERGED}" \
        --keep "${CV_VALID_FILE}" \
        --score "${SCORE_FILE}" 1 2 3 header cols=+scoresums \
        --out "${PGS_DIR}/cv_valid_tmp_${fold}" \
        > "${OUTDIR}/logs/pgs_cv_valid_${fold}.log" 2>&1

    ${PLINK2} --bfile "${GENO_MERGED}" \
        --keep "${TEST_SAMPLES}" \
        --score "${SCORE_FILE}" 1 2 3 header cols=+scoresums \
        --out "${PGS_DIR}/test_tmp_${fold}" \
        > "${OUTDIR}/logs/pgs_test_${fold}.log" 2>&1

    for pgs_type in cv_valid test; do
        SSCORE="${PGS_DIR}/${pgs_type}_tmp_${fold}.sscore"
        PROFILE="${PGS_DIR}/${pgs_type}_PGS_fold_${fold}.profile"

        if [[ -f "${SSCORE}" ]]; then
            Rscript --vanilla - <<RSCRIPT
library(data.table)
sscore <- fread("${SSCORE}")
if ("#FID" %in% colnames(sscore)) setnames(sscore, "#FID", "FID")

score_col <- grep("SCORE.*SUM", colnames(sscore), value=TRUE)[1]

profile <- data.table(
    FID = sscore[["FID"]],
    IID = sscore[["IID"]],
    PHENO = NA,
    CNT = sscore[["ALLELE_CT"]],
    CNT2 = NA,
    SCORESUM = sscore[[score_col]]
)
fwrite(profile, "${PROFILE}", sep="\t")
cat("  Created: ${PROFILE} (n =", nrow(profile), ")\n")
RSCRIPT
            rm -f "${SSCORE}"
        else
            echo "  WARNING: ${SSCORE} not found"
        fi
    done
done

rm -f "${PGS_DIR}"/*_tmp_* 2>/dev/null || true
