#!/usr/bin/env bash
#SBATCH -J PI_multiconf
#SBATCH --partition=tier2q
#SBATCH --time=24:00:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/multiconf_%A_%a.out
#SBATCH -e logs/multiconf_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SCALE=${1:? "Usage: sbatch multi_conf.sh <scale>  (original or bc)"}

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
TRAIT_STR="${TRAIT_INFO[$TASK_ID]}"
IFS=':' read -r phenoLower phenoFull phenoCode <<< "${TRAIT_STR}"

if [[ "${SCALE}" == "original" ]]; then
    OUTDIR="${OUTBASE}/original_scale/${phenoLower}"
elif [[ "${SCALE}" == "bc" ]]; then
    OUTDIR="${OUTBASE}/bc_scale/${phenoLower}"
else
    echo "ERROR: SCALE must be 'original' or 'bc'"
    exit 1
fi

echo "${phenoLower} (${SCALE})"

PGS_DIR="${OUTDIR}/pgs"
PI_SCRIPT="${OUTBASE}/PredInterval/src/PredInterval.R"
MULTICONF_DIR="${OUTDIR}/multiconf"
mkdir -p "${MULTICONF_DIR}"

CONF_LEVELS="0.75 0.80 0.85 0.90 0.95"

CV_PHENO="${OUTDIR}/cv_phenotypes.txt"
CV_PHENO_2COL="${OUTDIR}/cv_phenotypes_2col.txt"
TEST_PHENO="${OUTDIR}/test_phenotypes.txt"
TEST_FAM="${OUTDIR}/pgs/test.fam"

for f in "${PI_SCRIPT}" "${CV_PHENO}" "${TEST_PHENO}" "${TEST_FAM}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing $f"
        exit 1
    fi
done

for k in $(seq 1 ${K}); do
    for prefix in cv_valid_PGS test_PGS; do
        if [[ ! -f "${PGS_DIR}/${prefix}_fold_${k}.profile" ]]; then
            echo "ERROR: Missing ${PGS_DIR}/${prefix}_fold_${k}.profile"
            exit 1
        fi
    done
done

module load gcc/12.1.0
module load R/4.4.1

awk -v OFS='\t' '{print $1, $2}' "${CV_PHENO}" > "${CV_PHENO_2COL}"

for k in $(seq 1 ${K}); do
    ln -sf "${PGS_DIR}/cv_valid_PGS_fold_${k}.profile" \
           "${PGS_DIR}/cv_valid_PGS_subset_${k}.profile"
    ln -sf "${PGS_DIR}/test_PGS_fold_${k}.profile" \
           "${PGS_DIR}/test_PGS_subset_${k}.profile"
done

for CONF in ${CONF_LEVELS}; do
    echo "conf=${CONF}"

    OUTPUT_PI="${MULTICONF_DIR}/${phenoLower}_conf${CONF}_prediction_intervals.txt"
    OUTPUT_FULL="${MULTICONF_DIR}/${phenoLower}_conf${CONF}_full_results.txt"

    if [[ -f "${OUTPUT_FULL}" ]]; then
        continue
    fi

    Rscript "${PI_SCRIPT}" \
        "${CV_PHENO_2COL}" \
        "${PGS_DIR}/cv_valid_PGS_subset" \
        "${TEST_FAM}" \
        "${PGS_DIR}/test_PGS_subset" \
        "${K}" \
        "${OUTPUT_PI}" \
        "${CONF}"

    if [[ ! -f "${OUTPUT_PI}" ]]; then
        echo "  WARNING: PredInterval output not created for conf=${CONF}"
        continue
    fi

    Rscript --vanilla - "${OUTPUT_PI}" "${TEST_PHENO}" "${OUTPUT_FULL}" \
        "${PGS_DIR}" "${K}" "${CONF}" <<'RSCRIPT'
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
pi_file    <- args[1]
test_file  <- args[2]
out_full   <- args[3]
pgs_dir    <- args[4]
K          <- as.integer(args[5])
conf_level <- as.numeric(args[6])

pi_out <- fread(pi_file)
setnames(pi_out, names(pi_out)[1:3], c("id", "lower_bound", "upper_bound"))

test_pheno <- fread(test_file)

test_pgs_list <- vector("list", K)
for (k in 1:K) {
    f <- file.path(pgs_dir, paste0("test_PGS_subset_", k, ".profile"))
    pgs_k <- fread(f)
    test_pgs_list[[k]] <- pgs_k[, .(IID, value = SCORESUM)]
    setnames(test_pgs_list[[k]], "value", paste0("PGS_fold_", k))
}
test_pgs <- Reduce(function(a, b) merge(a, b, by = "IID", all = FALSE), test_pgs_list)
test_pgs[, pgs := rowMeans(.SD), .SDcols = paste0("PGS_fold_", 1:K)]

full_data <- merge(test_pheno, pi_out, by = "id", all = FALSE)
full_data <- merge(full_data, test_pgs[, .(IID, pgs)], by.x = "id", by.y = "IID", all = FALSE)

full_data[, width := upper_bound - lower_bound]
full_data[, covered := (y >= lower_bound) & (y <= upper_bound)]
full_data[, conf_level := conf_level]

fwrite(full_data[, .(id, y, pgs, lower_bound, upper_bound, width, covered, conf_level)],
       out_full, sep = "\t")

cat("  Conf:", conf_level,
    "Coverage:", round(mean(full_data$covered, na.rm = TRUE), 4),
    "Width:", round(mean(full_data$width, na.rm = TRUE), 4), "\n")
RSCRIPT

done

