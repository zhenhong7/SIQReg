#!/usr/bin/env bash
#SBATCH -J CA_setup
#SBATCH --partition=tier3q
#SBATCH --time=2:00:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/setup_%A_%a.out
#SBATCH -e logs/setup_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config_cross_ancestry.sh"

ANCESTRY=${1:? "Usage: sbatch cross_setup.sh <ancestry> <scale>"}
SCALE=${2:? "Usage: sbatch cross_setup.sh <ancestry> <scale>"}

if [[ -z "${ANCESTRY_IDS[$ANCESTRY]+x}" ]]; then
    echo "ERROR: Unknown ancestry '${ANCESTRY}'. Use: asn, afr, or white_euro"
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
mkdir -p "${OUTDIR}/logs"

ANCESTRY_ID_FILE="${ANCESTRY_IDS[$ANCESTRY]}"
CALIB_FRAC="${CALIB_RATIO[$ANCESTRY]}"

if [[ "${SCALE}" == "original" ]]; then
    PHENO_FILE="${PHENO_ORIGINAL_DIR}/${phenoFull}/${phenoFull}${phenoCode}.pheno"
else
    PHENO_FILE="${PHENO_BC_DIRS[$ANCESTRY]}/${phenoFull}${phenoCode}_bc.pheno"
fi

for f in "${ANCESTRY_ID_FILE}" "${PHENO_FILE}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: File not found: $f"
        exit 1
    fi
done

module load gcc/12.1.0
module load R/4.4.1

Rscript --vanilla - <<RSCRIPT
library(data.table)

ancestry_id_file <- "${ANCESTRY_ID_FILE}"
pheno_file <- "${PHENO_FILE}"
outdir <- "${OUTDIR}"
calib_frac <- ${CALIB_FRAC}
ancestry <- "${ANCESTRY}"
scale <- "${SCALE}"
phenoLower <- "${phenoLower}"

anc_ids <- fread(ancestry_id_file, header = FALSE, col.names = c("FID", "IID"))
cat("Ancestry IDs loaded:", nrow(anc_ids), "\n")

pheno <- fread(pheno_file)
orig_cols <- colnames(pheno)
if (grepl("^#", orig_cols[1])) setnames(pheno, orig_cols[1], sub("^#", "", orig_cols[1]))

if (ncol(pheno) >= 3) {
    colnames(pheno)[1:2] <- c("FID", "IID")
    pheno_col <- colnames(pheno)[3]
    pheno <- pheno[, .(FID, IID, PHENO = get(pheno_col))]
}

pheno[, PHENO := as.numeric(as.character(PHENO))]
pheno <- pheno[!is.na(PHENO)]

pheno_anc <- pheno[IID %in% anc_ids[["IID"]]]
cat("Individuals with phenotype in", ancestry, ":", nrow(pheno_anc), "\n")

if (nrow(pheno_anc) < 50) {
    cat("ERROR: Too few individuals (", nrow(pheno_anc), ") for", phenoLower, "\n")
    quit(status = 1)
}

set.seed(2024)
n <- nrow(pheno_anc)
n_calib <- round(n * calib_frac)

calib_idx <- sample(1:n, n_calib, replace = FALSE)
test_idx <- setdiff(1:n, calib_idx)

pheno_calib <- pheno_anc[calib_idx]
pheno_test <- pheno_anc[test_idx]

cat("Calibration set:", nrow(pheno_calib), "\n")
cat("Test set:", nrow(pheno_test), "\n")

fwrite(pheno_calib[, .(FID, IID)], paste0(outdir, "/calib_samples.txt"),
       sep = "\t", col.names = FALSE)
fwrite(pheno_test[, .(FID, IID)], paste0(outdir, "/test_samples.txt"),
       sep = "\t", col.names = FALSE)
fwrite(pheno_anc[, .(FID, IID)], paste0(outdir, "/all_samples.txt"),
       sep = "\t", col.names = FALSE)

fwrite(pheno_calib[, .(id = IID, y = PHENO)], paste0(outdir, "/calib_phenotypes.txt"),
       sep = "\t")
fwrite(pheno_test[, .(id = IID, y = PHENO)], paste0(outdir, "/test_phenotypes.txt"),
       sep = "\t")

stats <- data.table(
    ancestry = ancestry, scale = scale, trait = phenoLower,
    n_ancestry_total = nrow(anc_ids), n_with_pheno = nrow(pheno_anc),
    n_calib = nrow(pheno_calib), n_test = nrow(pheno_test),
    calib_frac = calib_frac,
    pheno_mean = mean(pheno_anc\$PHENO), pheno_sd = sd(pheno_anc\$PHENO)
)
fwrite(stats, paste0(outdir, "/pheno_stats.txt"), sep = "\t")

RSCRIPT
