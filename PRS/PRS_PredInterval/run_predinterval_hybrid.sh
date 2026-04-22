#!/usr/bin/env bash
#SBATCH -J PI_hybrid
#SBATCH --partition=tier3q
#SBATCH --time=240:00:00
#SBATCH --mem=20GB
#SBATCH --cpus-per-task=1
#SBATCH -o logs/hybrid_%A_%a.out
#SBATCH -e logs/hybrid_%A_%a.err
#SBATCH --array=1-25


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
TRAIT_STR="${TRAIT_INFO[$TASK_ID]}"
IFS=':' read -r phenoLower phenoFull phenoCode <<< "${TRAIT_STR}"

echo "${phenoLower}: hybrid decomposition"

ORIG_DIR="${OUTBASE}/original_scale/${phenoLower}"
BC_DIR="${OUTBASE}/bc_scale/${phenoLower}"
HYBRID_B_DIR="${OUTBASE}/hybrid_bc2orig/${phenoLower}"
HYBRID_D_DIR="${OUTBASE}/hybrid_orig2bc/${phenoLower}"

mkdir -p "${HYBRID_B_DIR}"/{pgs,logs}
mkdir -p "${HYBRID_D_DIR}"/{pgs,logs}

PI_SCRIPT="${OUTBASE}/PredInterval/src/PredInterval.R"

LAMBDA_META="${LAMBDA_META:-${PROJECT_DIR}/SI_GWAS/meta_both_summary_all.txt}"

if [[ ! -f "${LAMBDA_META}" ]]; then
    echo "ERROR: Lambda meta file not found: ${LAMBDA_META}"
    exit 1
fi

LAMBDA_VAL=$(awk -F'\t' -v trait="${phenoFull}" 'NR>1 && $1==trait {print $3}' "${LAMBDA_META}")

if [[ -z "${LAMBDA_VAL}" ]]; then
    echo "ERROR: Could not find lambda for ${phenoFull} in ${LAMBDA_META}"
    exit 1
fi

echo "Lambda for ${phenoFull}: ${LAMBDA_VAL}"

for k in $(seq 1 ${K}); do
    for f in "${BC_DIR}/pgs/cv_valid_PGS_fold_${k}.profile" \
             "${BC_DIR}/pgs/test_PGS_fold_${k}.profile" \
             "${ORIG_DIR}/pgs/cv_valid_PGS_fold_${k}.profile" \
             "${ORIG_DIR}/pgs/test_PGS_fold_${k}.profile" ; do
        if [[ ! -f "$f" ]]; then
            echo "ERROR: Missing PGS file: $f"
            exit 1
        fi
    done
done

for f in "${ORIG_DIR}/cv_phenotypes.txt" "${BC_DIR}/cv_phenotypes.txt"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing phenotype file: $f"
        exit 1
    fi
done

module load gcc/12.1.0
module load R/4.4.1

TRANSFORM_SCRIPT=$(mktemp /tmp/transform_pgs_XXXXXX.R)
cat > "${TRANSFORM_SCRIPT}" <<'RSCRIPT'
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
src_pgs_dir    <- args[1]
dst_pgs_dir    <- args[2]
lambda         <- as.numeric(args[3])
K              <- as.integer(args[4])
direction      <- args[5]  # "inv" or "fwd"

cat("Lambda:", lambda, ", Direction:", direction, "\n")

fwd_boxcox <- function(z, lam) {
    if (abs(lam) < 1e-10) {
        return(log(pmax(z, 1e-10)))
    } else {
        return((pmax(z, 1e-10)^lam - 1) / lam)
    }
}

inv_boxcox <- function(z, lam) {
    if (abs(lam) < 1e-10) {
        return(exp(z))
    } else {
        base <- pmax(lam * z + 1, 1e-10)
        return(base^(1 / lam))
    }
}

transform_fn <- if (direction == "inv") inv_boxcox else fwd_boxcox

for (pgs_type in c("cv_valid", "test")) {
    for (k in 1:K) {
        infile  <- file.path(src_pgs_dir, paste0(pgs_type, "_PGS_fold_", k, ".profile"))
        outfile <- file.path(dst_pgs_dir, paste0(pgs_type, "_PGS_fold_", k, ".profile"))

        if (!file.exists(infile)) {
            cat("WARNING: Missing", infile, "\n")
            next
        }

        pgs <- fread(infile)
        pgs[, SCORESUM := transform_fn(SCORESUM, lambda)]
        fwrite(pgs, outfile, sep = "\t")
        cat("  Created:", outfile, "(n =", nrow(pgs), ")\n")
    }
}
cat("Transform complete.\n")
RSCRIPT

POSTPROC_SCRIPT=$(mktemp /tmp/postproc_pgs_XXXXXX.R)
cat > "${POSTPROC_SCRIPT}" <<'RSCRIPT'
library(data.table)
args <- commandArgs(trailingOnly = TRUE)
pi_file <- args[1]; test_file <- args[2]; out_full <- args[3]; out_sum <- args[4]
pgs_dir <- args[5]; trait_name <- args[6]; K <- as.integer(args[7])
conf_level <- as.numeric(args[8]); scale_label <- args[9]

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

fwrite(full_data[, .(id, y, pgs, lower_bound, upper_bound, width, covered)], out_full, sep = "\t")

coverage <- mean(full_data$covered, na.rm = TRUE)
mean_width <- mean(full_data$width, na.rm = TRUE)
summary_dt <- data.table(trait = trait_name, scale = scale_label,
    n_test = nrow(full_data), K = K, target_coverage = conf_level,
    observed_coverage = coverage, mean_width = mean_width,
    sd_width = sd(full_data$width, na.rm = TRUE),
    mean_pgs = mean(full_data$pgs, na.rm = TRUE),
    sd_pgs = sd(full_data$pgs, na.rm = TRUE))
fwrite(summary_dt, out_sum, sep = "\t")
cat("Coverage:", round(coverage, 4), ", Width:", round(mean_width, 4), "\n")
RSCRIPT

# Cell B: g^{-1}(PGS_bc) + Y_orig
Rscript --vanilla "${TRANSFORM_SCRIPT}" \
    "${BC_DIR}/pgs" "${HYBRID_B_DIR}/pgs" "${LAMBDA_VAL}" "${K}" "inv"

awk -v OFS='\t' '{print $1, $2}' "${ORIG_DIR}/cv_phenotypes.txt" \
    > "${HYBRID_B_DIR}/cv_phenotypes_2col.txt"

for k in $(seq 1 ${K}); do
    ln -sf "${HYBRID_B_DIR}/pgs/cv_valid_PGS_fold_${k}.profile" \
           "${HYBRID_B_DIR}/pgs/cv_valid_PGS_subset_${k}.profile"
    ln -sf "${HYBRID_B_DIR}/pgs/test_PGS_fold_${k}.profile" \
           "${HYBRID_B_DIR}/pgs/test_PGS_subset_${k}.profile"
done

Rscript "${PI_SCRIPT}" \
    "${HYBRID_B_DIR}/cv_phenotypes_2col.txt" \
    "${HYBRID_B_DIR}/pgs/cv_valid_PGS_subset" \
    "${ORIG_DIR}/pgs/test.fam" \
    "${HYBRID_B_DIR}/pgs/test_PGS_subset" \
    "${K}" \
    "${HYBRID_B_DIR}/${phenoLower}_prediction_intervals.txt" \
    "${CONF_LEVEL}"

if [[ ! -f "${HYBRID_B_DIR}/${phenoLower}_prediction_intervals.txt" ]]; then
    echo "ERROR: Cell B PredInterval output not created"; exit 1
fi

Rscript --vanilla "${POSTPROC_SCRIPT}" \
    "${HYBRID_B_DIR}/${phenoLower}_prediction_intervals.txt" \
    "${ORIG_DIR}/test_phenotypes.txt" \
    "${HYBRID_B_DIR}/${phenoLower}_full_results.txt" \
    "${HYBRID_B_DIR}/${phenoLower}_summary.txt" \
    "${HYBRID_B_DIR}/pgs" "${phenoLower}" "${K}" "${CONF_LEVEL}" "hybrid_bc2orig"

# Cell D: g(PGS_orig) + Y_bc
Rscript --vanilla "${TRANSFORM_SCRIPT}" \
    "${ORIG_DIR}/pgs" "${HYBRID_D_DIR}/pgs" "${LAMBDA_VAL}" "${K}" "fwd"

awk -v OFS='\t' '{print $1, $2}' "${BC_DIR}/cv_phenotypes.txt" \
    > "${HYBRID_D_DIR}/cv_phenotypes_2col.txt"

for k in $(seq 1 ${K}); do
    ln -sf "${HYBRID_D_DIR}/pgs/cv_valid_PGS_fold_${k}.profile" \
           "${HYBRID_D_DIR}/pgs/cv_valid_PGS_subset_${k}.profile"
    ln -sf "${HYBRID_D_DIR}/pgs/test_PGS_fold_${k}.profile" \
           "${HYBRID_D_DIR}/pgs/test_PGS_subset_${k}.profile"
done

Rscript "${PI_SCRIPT}" \
    "${HYBRID_D_DIR}/cv_phenotypes_2col.txt" \
    "${HYBRID_D_DIR}/pgs/cv_valid_PGS_subset" \
    "${BC_DIR}/pgs/test.fam" \
    "${HYBRID_D_DIR}/pgs/test_PGS_subset" \
    "${K}" \
    "${HYBRID_D_DIR}/${phenoLower}_prediction_intervals.txt" \
    "${CONF_LEVEL}"

if [[ ! -f "${HYBRID_D_DIR}/${phenoLower}_prediction_intervals.txt" ]]; then
    echo "ERROR: Cell D PredInterval output not created"; exit 1
fi

Rscript --vanilla "${POSTPROC_SCRIPT}" \
    "${HYBRID_D_DIR}/${phenoLower}_prediction_intervals.txt" \
    "${BC_DIR}/test_phenotypes.txt" \
    "${HYBRID_D_DIR}/${phenoLower}_full_results.txt" \
    "${HYBRID_D_DIR}/${phenoLower}_summary.txt" \
    "${HYBRID_D_DIR}/pgs" "${phenoLower}" "${K}" "${CONF_LEVEL}" "hybrid_orig2bc"

rm -f "${TRANSFORM_SCRIPT}" "${POSTPROC_SCRIPT}" 2>/dev/null || true
