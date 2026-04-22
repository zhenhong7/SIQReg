#!/usr/bin/env bash
#SBATCH -J PI_run
#SBATCH --partition=tier3q
#SBATCH --time=24:00:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/predint_%A_%a.out
#SBATCH -e logs/predint_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

SCALE=${1:? "Usage: sbatch run_predinterval.sh <scale>"}

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

CV_PHENO="${OUTDIR}/cv_phenotypes.txt"
CV_PHENO_2COL="${OUTDIR}/cv_phenotypes_2col.txt"
TEST_PHENO="${OUTDIR}/test_phenotypes.txt"
TEST_FAM="${OUTDIR}/pgs/test.fam"

OUTPUT_PI="${OUTDIR}/${phenoLower}_prediction_intervals.txt"
OUTPUT_SUMMARY="${OUTDIR}/${phenoLower}_summary.txt"
OUTPUT_FULL="${OUTDIR}/${phenoLower}_full_results.txt"

module load gcc/12.1.0
module load R/4.4.1

if [[ ! -f "${PI_SCRIPT}" ]]; then
    echo "ERROR: PredInterval.R not found at ${PI_SCRIPT}"
    exit 1
fi

if [[ ! -f "${CV_PHENO}" ]]; then
    echo "ERROR: Missing ${CV_PHENO}"
    exit 1
fi

if [[ ! -f "${TEST_PHENO}" ]]; then
    echo "ERROR: Missing ${TEST_PHENO}"
    exit 1
fi

if [[ ! -f "${TEST_FAM}" ]]; then
    echo "ERROR: Missing ${TEST_FAM}"
    exit 1
fi

for k in $(seq 1 ${K}); do
    if [[ ! -f "${PGS_DIR}/cv_valid_PGS_fold_${k}.profile" ]]; then
        echo "ERROR: Missing ${PGS_DIR}/cv_valid_PGS_fold_${k}.profile"
        exit 1
    fi
    if [[ ! -f "${PGS_DIR}/test_PGS_fold_${k}.profile" ]]; then
        echo "ERROR: Missing ${PGS_DIR}/test_PGS_fold_${k}.profile"
        exit 1
    fi
done

awk -v OFS='\t' '{print $1, $2}' "${CV_PHENO}" > "${CV_PHENO_2COL}"

for k in $(seq 1 ${K}); do
    ln -sf "${PGS_DIR}/cv_valid_PGS_fold_${k}.profile" "${PGS_DIR}/cv_valid_PGS_subset_${k}.profile"
    ln -sf "${PGS_DIR}/test_PGS_fold_${k}.profile" "${PGS_DIR}/test_PGS_subset_${k}.profile"
done


Rscript "${PI_SCRIPT}" \
    "${CV_PHENO_2COL}" \
    "${PGS_DIR}/cv_valid_PGS_subset" \
    "${TEST_FAM}" \
    "${PGS_DIR}/test_PGS_subset" \
    "${K}" \
    "${OUTPUT_PI}" \
    "${CONF_LEVEL}"


if [[ ! -f "${OUTPUT_PI}" ]]; then
    echo "ERROR: PredInterval output not created: ${OUTPUT_PI}"
    exit 1
fi


Rscript --vanilla - <<RSCRIPT
library(data.table)

pi_file    <- "${OUTPUT_PI}"
test_file  <- "${TEST_PHENO}"
out_full   <- "${OUTPUT_FULL}"
out_sum    <- "${OUTPUT_SUMMARY}"
pgs_dir    <- "${PGS_DIR}"
trait_name <- "${phenoLower}"
scale_name <- "${SCALE}"
K          <- as.integer("${K}")
conf_level <- as.numeric("${CONF_LEVEL}")
cv_file    <- "${CV_PHENO_2COL}"

pi_out <- fread(pi_file)

if (ncol(pi_out) < 3) {
    stop("PredInterval output has fewer than 3 columns: ", pi_file)
}
setnames(pi_out, names(pi_out)[1:3], c("id", "lower_bound", "upper_bound"))

test_pheno <- fread(test_file)

if (!all(c("id", "y") %in% names(test_pheno))) {
    stop("test_phenotypes.txt must contain columns named 'id' and 'y'")
}

test_pgs_list <- vector("list", K)

for (k in 1:K) {
    f <- file.path(pgs_dir, paste0("test_PGS_subset_", k, ".profile"))
    if (!file.exists(f)) stop("Missing PGS profile: ", f)

    pgs_k <- fread(f)

    if (!all(c("IID", "SCORESUM") %in% names(pgs_k))) {
        stop("Profile file missing IID or SCORESUM: ", f)
    }

    test_pgs_list[[k]] <- pgs_k[, .(IID, value = SCORESUM)]
    setnames(test_pgs_list[[k]], "value", paste0("PGS_fold_", k))
}

test_pgs <- test_pgs_list[[1]]
if (K > 1) {
    for (k in 2:K) {
        test_pgs <- merge(test_pgs, test_pgs_list[[k]], by = "IID", all = FALSE)
    }
}
test_pgs[, pgs := rowMeans(.SD), .SDcols = paste0("PGS_fold_", 1:K)]

full_data <- merge(test_pheno, pi_out, by = "id", all = FALSE)
full_data <- merge(full_data, test_pgs[, .(IID, pgs)], by.x = "id", by.y = "IID", all = FALSE)

full_data[, width := upper_bound - lower_bound]
full_data[, covered := (y >= lower_bound) & (y <= upper_bound)]

fwrite(
    full_data[, .(id, y, pgs, lower_bound, upper_bound, width, covered)],
    out_full,
    sep = "\t"
)

coverage   <- mean(full_data[["covered"]], na.rm = TRUE)
mean_width <- mean(full_data[["width"]], na.rm = TRUE)

summary_dt <- data.table(
    trait = trait_name,
    scale = scale_name,
    n_cv = nrow(fread(cv_file)),
    n_test = nrow(full_data),
    K = K,
    target_coverage = conf_level,
    observed_coverage = coverage,
    mean_width = mean_width,
    sd_width = sd(full_data[["width"]], na.rm = TRUE),
    mean_pgs = mean(full_data[["pgs"]], na.rm = TRUE),
    sd_pgs = sd(full_data[["pgs"]], na.rm = TRUE)
)

fwrite(summary_dt, out_sum, sep = "\t")

cat("Target coverage:", conf_level, "\\n")
cat("Observed coverage:", round(coverage, 4), "\\n")
cat("Mean interval width:", round(mean_width, 4), "\\n")
RSCRIPT

