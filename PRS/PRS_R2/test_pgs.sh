#!/usr/bin/env bash

#SBATCH -J PGS
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=75GB
#SBATCH --partition=tier2q
#SBATCH -o logs/PGS_PGSC_%A_%a.out
#SBATCH -e logs/PGS_PGSC_%A_%a.err
#SBATCH --array=1-25

now=$(date)
start_time=$(date +%s)
echo "Date run: $now"

lambda=1
GENOS="${GENO_BASE:-/path/to/genotypes/pop_genos}/"
PHENOS="${PHENO_DIR:-/path/to/extracted_phenotypes}/"
LOCAL_PLINK2="${PLINK2:-plink2}"
TRAIN_POP="whitebrit"
TEST_POP="whitebrit"
JOINED_GWAS="${JOINED_GWAS_DIR:-/path/to/SI_GWAS/joined_gwas_results}/"
pop_path="${GENOS}${TEST_POP}/"
genos="${GENOS}${TEST_POP}/"
pop_file="ukb_chr1-22_${TEST_POP}"
PGS="${PGS_OUT_DIR:-/path/to/SI_GWAS/pgs_outputs/pgs_out}/"
CLUMPED_GENOS="${CLUMPED_GENOS_DIR:-/path/to/SI_GWAS/pgs_output/clumped_genos}/"

RANGE_LIST="${RANGE_LIST:-/path/to/SI_GWAS/std_range_list}"
if [[ ! -s "$RANGE_LIST" ]]; then
  echo "FATAL: range list not found or empty: $RANGE_LIST" >&2
  exit 1
fi

THRESH=(0.0000000001 0.00000001 0.000001 0.0001 0.001 0.005 0.01 0.05 0.1 0.5)
PHENOARRAY=(
  "arm_fat-free_mass_left"
  "arm_fat-free_mass_right"
  "birth_weight"
  "bmi"
  "c-reactive_protein"
  "calcium"
  "creatinine"
  "diastolicbp_auto"
  "fev1"
  "hba1c"
  "hdl"
  "height"
  "hip_to_waist"
  "igf-1"
  "ldl"
  "leukocyte_count"
  "platelet_count"
  "pulse_rate"
  "rbc"
  "shbg"
  "systolicbp_auto"
  "testosterone"
  "urate"
  "urea"
  "vitamin_d"
)

if (( $(echo "$lambda == 0" | bc -l) )); then
  JOINED_GWAS="${JOINED_GWAS_LOG_DIR:-/path/to/SI_GWAS/joined_gwas_results_log}/"
  PGS="${PGS_OUT_LOG_DIR:-/path/to/SI_GWAS/pgs_outputs/pgs_out_log}/"
  CLUMPED_GENOS="${CLUMPED_GENOS_LOG_DIR:-/path/to/SI_GWAS/pgs_output/clumped_genos_log}/"
fi

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "PHENOARRAY length: ${#PHENOARRAY[@]}"

pheno_num=$((SLURM_ARRAY_TASK_ID-1))
echo "Calculated pheno_num: $pheno_num"

selected_pheno="${PHENOARRAY[pheno_num]}"
echo "$selected_pheno"
if [[ ! $selected_pheno ]]; then
  echo "no pheno present"
  exit 1
fi

pgs_path="${PGS}${selected_pheno}/"

mkdir -p "${CLUMPED_GENOS}"
mkdir -p "${pgs_path}"
mkdir -p "${JOINED_GWAS}snp_pval/"

echo "Building std ${selected_pheno} PGS"

# clump snps
output_file="${CLUMPED_GENOS}${pop_file}_${selected_pheno}"

gwas_file=$(ls -1 "${JOINED_GWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}."*.glm.linear 2>/dev/null | head -n 1 || true)
if [[ -z "$gwas_file" ]]; then
  echo "FATAL: joined GWAS file not found for ${selected_pheno}" >&2
  exit 1
fi
echo "Using GWAS file: $gwas_file"

valid_snp_file="${CLUMPED_GENOS}${pop_file}_${selected_pheno}.valid.snp"

if [[ ! -f "${output_file}.clumps" ]]; then
  ${LOCAL_PLINK2} \
    --bfile "${pop_path}${pop_file}_QC" \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump "$gwas_file" \
    --out "$output_file"
  echo "$selected_pheno clump done"
fi

# (re)generate valid SNP list from clumps
awk 'NR>1{print $3}' "${output_file}.clumps" > "$valid_snp_file"

# build pgs
pgs_out_file="${pgs_path}${selected_pheno}_pgs"
last_thresh="${THRESH[${#THRESH[@]}-1]}"

if [[ ! -f "${pgs_out_file}.${last_thresh}.sscore" ]]; then
  ${LOCAL_PLINK2} \
    --bfile "${pop_path}${pop_file}_QC" \
    --keep "${genos}${pop_file}_QC_test" \
    --score "$gwas_file" 3 7 12 header cols=+scoresums \
    --q-score-range "$RANGE_LIST" "$gwas_file" 3 15 header \
    --extract "$valid_snp_file" \
    --out "$pgs_out_file"
  echo "Std PGS complete"
else
  echo "std pgs already exists"
fi

echo "test PGS done"

end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds"
