#!/usr/bin/env bash
#SBATCH -J CA_pgs
#SBATCH --partition=tier3q
#SBATCH --time=12:00:00
#SBATCH --mem=30GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/pgs_%A_%a.out
#SBATCH -e logs/pgs_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config_cross_ancestry.sh"

ANCESTRY=${1:? "Usage: sbatch cross_score_pgs.sh <ancestry> <scale>"}
SCALE=${2:? "Usage: sbatch cross_score_pgs.sh <ancestry> <scale>"}

if [[ -z "${ANCESTRY_IDS[$ANCESTRY]+x}" ]]; then
    echo "ERROR: Unknown ancestry '${ANCESTRY}'"
    exit 1
fi
if [[ "${SCALE}" != "original" && "${SCALE}" != "bc" ]]; then
    echo "ERROR: scale must be 'original' or 'bc'"
    exit 1
fi

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
TRAIT_STR="${TRAIT_INFO[$TASK_ID]}"
IFS=':' read -r phenoLower phenoFull phenoCode <<< "${TRAIT_STR}"

echo "${ANCESTRY} ${phenoLower} ${SCALE}"

OUTDIR="${OUTBASE}/${ANCESTRY}/${SCALE}_scale/${phenoLower}"
PGS_DIR="${OUTDIR}/pgs"
mkdir -p "${PGS_DIR}"

WB_GWAS_DIR="${WB_RESULTS}/${SCALE}_scale/${phenoLower}/gwas"
GENO_DIR="${GENO_DIRS[$ANCESTRY]}"
GENO_TAG="${GENO_SUFFIX[$ANCESTRY]}"
ALL_SAMPLES="${OUTDIR}/all_samples.txt"

if [[ ! -f "${ALL_SAMPLES}" ]]; then
    echo "ERROR: ${ALL_SAMPLES} not found. Run cross_setup.sh first."
    exit 1
fi

module load gcc/12.1.0
module load R/4.4.1

for fold in $(seq 1 ${K}); do
    SCORE_FILE="${WB_GWAS_DIR}/score_fold_${fold}.txt"

    if [[ ! -f "${SCORE_FILE}" ]]; then
        echo "WARNING: Score file not found: ${SCORE_FILE}, skipping fold ${fold}"
        continue
    fi

    echo ""
    echo "Fold ${fold}"
    echo "  Score file: ${SCORE_FILE}"

    N_SNPS=$(tail -n +2 "${SCORE_FILE}" | wc -l)
    echo "  SNPs in score file: ${N_SNPS}"

    for chr in $(seq 1 22); do
        GENO_FILE="${GENO_DIR}/ukb_chr${chr}_${GENO_TAG}_QC"

        if [[ ! -f "${GENO_FILE}.bed" ]]; then
            echo "  WARNING: ${GENO_FILE}.bed not found, skipping chr${chr}"
            continue
        fi

        CHR_OUT="${PGS_DIR}/tmp_fold${fold}_chr${chr}"

        ${PLINK2} --bfile "${GENO_FILE}" \
            --keep "${ALL_SAMPLES}" \
            --score "${SCORE_FILE}" 1 2 3 header cols=+scoresums \
            --out "${CHR_OUT}" \
            > "${OUTDIR}/logs/pgs_fold${fold}_chr${chr}.log" 2>&1 || true
    done

    echo "  Merging PGS across chromosomes for fold ${fold}..."

    Rscript --vanilla - <<RSCRIPT
library(data.table)

pgs_dir <- "${PGS_DIR}"
fold <- ${fold}

chr_scores <- list()
for (chr in 1:22) {
    sscore_file <- paste0(pgs_dir, "/tmp_fold", fold, "_chr", chr, ".sscore")
    if (file.exists(sscore_file)) {
        dt <- fread(sscore_file)
        if ("#FID" %in% colnames(dt)) setnames(dt, "#FID", "FID")
        score_col <- grep("SCORE.*SUM", colnames(dt), value = TRUE)[1]
        if (!is.na(score_col)) {
            chr_scores[[chr]] <- dt[, .(IID, score = get(score_col))]
        }
    }
}

if (length(chr_scores) == 0) {
    cat("ERROR: No chromosome scores found for fold", fold, "\n")
    quit(status = 1)
}

merged <- rbindlist(chr_scores)
pgs_total <- merged[, .(SCORESUM = sum(score, na.rm = TRUE)), by = IID]

out_file <- paste0(pgs_dir, "/pgs_fold_", fold, ".txt")
fwrite(pgs_total, out_file, sep = "\t")
cat("  Fold", fold, ": n =", nrow(pgs_total), "individuals scored\n")
RSCRIPT

    rm -f "${PGS_DIR}"/tmp_fold${fold}_chr*.sscore 2>/dev/null || true
    rm -f "${PGS_DIR}"/tmp_fold${fold}_chr*.log 2>/dev/null || true
done

echo ""
echo "Averaging PGS across ${K} folds..."

Rscript --vanilla - <<RSCRIPT
library(data.table)

pgs_dir <- "${PGS_DIR}"
K <- ${K}

pgs_list <- list()
for (k in 1:K) {
    f <- paste0(pgs_dir, "/pgs_fold_", k, ".txt")
    if (file.exists(f)) {
        dt <- fread(f)
        if (nrow(dt) > 0 && "SCORESUM" %in% colnames(dt)) {
            setnames(dt, "SCORESUM", paste0("PGS_fold_", k))
            pgs_list[[paste0("fold_", k)]] <- dt
        } else {
            cat("  WARNING: Fold", k, "file exists but is empty or malformed, skipping\n")
        }
    } else {
        cat("  WARNING: Fold", k, "PGS file not found (WB score file likely missing), skipping\n")
    }
}

if (length(pgs_list) == 0) {
    cat("ERROR: No fold PGS files found\n")
    quit(status = 1)
}

cat("  Averaging across", length(pgs_list), "available folds\n")

pgs_merged <- pgs_list[[1]]
if (length(pgs_list) > 1) {
    for (k in 2:length(pgs_list)) {
        pgs_merged <- merge(pgs_merged, pgs_list[[k]], by = "IID", all = TRUE)
    }
}

pgs_cols <- grep("^PGS_fold_", colnames(pgs_merged), value = TRUE)
pgs_merged[, PGS := rowMeans(.SD, na.rm = TRUE), .SDcols = pgs_cols]

out_file <- paste0(pgs_dir, "/pgs_averaged.txt")
fwrite(pgs_merged[, .(IID, PGS)], out_file, sep = "\t")
cat("Final averaged PGS saved: n =", nrow(pgs_merged), "\n")
cat("  Mean PGS:", round(mean(pgs_merged\$PGS, na.rm = TRUE), 6), "\n")
cat("  SD PGS:", round(sd(pgs_merged\$PGS, na.rm = TRUE), 6), "\n")
RSCRIPT

rm -f "${PGS_DIR}"/tmp_* 2>/dev/null || true

echo ""
echo "PGS scoring complete: ${phenoLower} (${ANCESTRY}, ${SCALE})"
echo "  Output: ${PGS_DIR}/pgs_averaged.txt"
