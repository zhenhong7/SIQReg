#!/usr/bin/env bash
#SBATCH -J CA_pi
#SBATCH --partition=tier3q
#SBATCH --time=4:00:00
#SBATCH --mem=10GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/pi_%A_%a.out
#SBATCH -e logs/pi_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config_cross_ancestry.sh"

ANCESTRY=${1:? "Usage: sbatch cross_predinterval.sh <ancestry> <scale>"}
SCALE=${2:? "Usage: sbatch cross_predinterval.sh <ancestry> <scale>"}

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
PGS_FILE="${PGS_DIR}/pgs_averaged.txt"
CALIB_PHENO="${OUTDIR}/calib_phenotypes.txt"
TEST_PHENO="${OUTDIR}/test_phenotypes.txt"

# OUTBASE points to .../cross_ancestry, so PredInterval is one directory up
PI_SCRIPT="${OUTBASE}/../PredInterval/src/PredInterval.R"

CALIB_PHENO_2COL="${OUTDIR}/calib_phenotypes_2col.txt"
TEST_FAM="${OUTDIR}/test.fam"

CALIB_PROFILE="${PGS_DIR}/calib_PGS_subset_1.profile"
TEST_PROFILE="${PGS_DIR}/test_PGS_subset_1.profile"

OUTPUT_PI="${OUTDIR}/${phenoLower}_prediction_intervals.txt"
OUTPUT_SUMMARY="${OUTDIR}/${phenoLower}_summary.txt"
OUTPUT_FULL="${OUTDIR}/${phenoLower}_full_results.txt"

module load gcc/12.1.0
module load R/4.4.1

for f in "${CALIB_PHENO}" "${TEST_PHENO}" "${PGS_FILE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Missing ${f}"
        exit 1
    fi
done

if [[ ! -f "${PI_SCRIPT}" ]]; then
    echo "ERROR: PredInterval.R not found at ${PI_SCRIPT}"
    exit 1
fi

mkdir -p "${PGS_DIR}"


Rscript --vanilla - <<RSCRIPT
library(data.table)

pgs_file   <- "${PGS_FILE}"
calib_file <- "${CALIB_PHENO}"
test_file  <- "${TEST_PHENO}"

out_calib_pheno <- "${CALIB_PHENO_2COL}"
out_test_fam    <- "${TEST_FAM}"
out_calib_prof  <- "${CALIB_PROFILE}"
out_test_prof   <- "${TEST_PROFILE}"

pgs <- fread(pgs_file)
calib_pheno <- fread(calib_file)
test_pheno  <- fread(test_file)

if (!("IID" %in% names(pgs))) {
    stop("pgs_averaged.txt must contain column 'IID'")
}

pgs_col <- NULL
if ("PGS" %in% names(pgs)) {
    pgs_col <- "PGS"
} else if ("pgs" %in% names(pgs)) {
    pgs_col <- "pgs"
} else {
    stop("pgs_averaged.txt must contain either column 'PGS' or 'pgs'")
}

if (!all(c("id","y") %in% names(calib_pheno))) {
    stop("calib_phenotypes.txt must contain columns 'id' and 'y'")
}
if (!all(c("id","y") %in% names(test_pheno))) {
    stop("test_phenotypes.txt must contain columns 'id' and 'y'")
}

setnames(pgs, pgs_col, "PGS")

calib_data <- merge(calib_pheno[, .(id, y)],
                    pgs[, .(IID, PGS)],
                    by.x = "id", by.y = "IID",
                    all = FALSE)

if (nrow(calib_data) == 0) {
    stop("No overlapping IDs between calibration phenotypes and averaged PGS")
}

shift <- mean(calib_data[["y"]], na.rm = TRUE) - mean(calib_data[["PGS"]], na.rm = TRUE)
cat("Calculated PGS shift for ${ANCESTRY}: ", round(shift, 6), "\\n", sep = "")

pgs[, shifted_PGS := PGS + shift]

fwrite(calib_pheno[, .(id, y)], out_calib_pheno, sep = "\t")

# official PredInterval uses the 2nd column as IID
test_fam <- test_pheno[, .(FID = id, IID = id)]
fwrite(test_fam, out_test_fam, sep = "\t", col.names = FALSE)

calib_prof <- pgs[IID %in% calib_pheno[["id"]], .(IID, SCORESUM = shifted_PGS)]
test_prof  <- pgs[IID %in% test_pheno[["id"]],  .(IID, SCORESUM = shifted_PGS)]

if (nrow(calib_prof) == 0) stop("Calibration profile has 0 rows")
if (nrow(test_prof)  == 0) stop("Test profile has 0 rows")

fwrite(calib_prof, out_calib_prof, sep = "\t")
fwrite(test_prof,  out_test_prof,  sep = "\t")
RSCRIPT

for f in "${CALIB_PHENO_2COL}" "${TEST_FAM}" "${CALIB_PROFILE}" "${TEST_PROFILE}"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: Expected generated file not found: ${f}"
        exit 1
    fi
done

Rscript "${PI_SCRIPT}" \
    "${CALIB_PHENO_2COL}" \
    "${PGS_DIR}/calib_PGS_subset" \
    "${TEST_FAM}" \
    "${PGS_DIR}/test_PGS_subset" \
    1 \
    "${OUTPUT_PI}" \
    "${CONF_LEVEL}"

if [[ ! -f "${OUTPUT_PI}" ]]; then
    echo "ERROR: PredInterval output not created: ${OUTPUT_PI}"
    exit 1
fi

Rscript --vanilla - <<RSCRIPT
library(data.table)

pi_file       <- "${OUTPUT_PI}"
test_file     <- "${TEST_PHENO}"
test_prof     <- "${TEST_PROFILE}"
out_full      <- "${OUTPUT_FULL}"
out_sum       <- "${OUTPUT_SUMMARY}"
trait_name    <- "${phenoLower}"
scale_name    <- "${SCALE}"
ancestry_name <- "${ANCESTRY}"
conf_level    <- as.numeric("${CONF_LEVEL}")
calib_file    <- "${CALIB_PHENO_2COL}"

pi_out <- fread(pi_file)
if (ncol(pi_out) < 3) {
    stop("PredInterval output has fewer than 3 columns")
}
setnames(pi_out, names(pi_out)[1:3], c("id", "lower_bound", "upper_bound"))

test_pheno <- fread(test_file)
if (!all(c("id","y") %in% names(test_pheno))) {
    stop("test_phenotypes.txt must contain columns 'id' and 'y'")
}

test_pgs <- fread(test_prof)
if (!all(c("IID","SCORESUM") %in% names(test_pgs))) {
    stop("test profile must contain columns 'IID' and 'SCORESUM'")
}
setnames(test_pgs, "SCORESUM", "pgs")

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

cor_y_pgs <- NA_real_
if (nrow(full_data) > 1) {
    cor_y_pgs <- suppressWarnings(cor(full_data[["y"]], full_data[["pgs"]], use = "complete.obs"))
}

summary_dt <- data.table(
    trait = trait_name,
    scale = scale_name,
    ancestry = ancestry_name,
    n_calib = nrow(fread(calib_file)),
    n_test = nrow(full_data),
    target_coverage = conf_level,
    observed_coverage = coverage,
    mean_width = mean_width,
    sd_width = sd(full_data[["width"]], na.rm = TRUE),
    mean_pgs = mean(full_data[["pgs"]], na.rm = TRUE),
    sd_pgs = sd(full_data[["pgs"]], na.rm = TRUE),
    cor_y_pgs = cor_y_pgs
)

fwrite(summary_dt, out_sum, sep = "\t")

cat("Target coverage:", conf_level, "\\n")
cat("Observed coverage:", round(coverage, 4), "\\n")
cat("Mean interval width:", round(mean_width, 4), "\\n")
RSCRIPT

