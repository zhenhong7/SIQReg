#!/usr/bin/env bash
#SBATCH -J val_PGS_rank
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --partition=tier2q
#SBATCH -o logs/val_PGS_rank_%A_%a.out
#SBATCH -e logs/val_PGS_rank_%A_%a.err
#SBATCH --array=1-25

set -euo pipefail

echo "Date run: $(date)"
start_time=$(date +%s)

# INPUTS
GENOS="${GENO_BASE:-/path/to/genotypes/pop_genos}/"
LOCAL_PLINK2="${PLINK2:-plink2}"

TRAIN_POP="whitebrit"   # GWAS source pop
TEST_POP="whitebrit"    # where you chose best thresh + clumped SNPs

JOINED_GWAS="${JOINED_GWAS_DIR:-/path/to/SI_GWAS/joined_gwas_results}/"
CLUMPED_GENOS="${CLUMPED_GENOS_DIR:-/path/to/SI_GWAS/pgs_output/clumped_genos}/"

META_TAB="${META_TAB:-/path/to/SI_GWAS/meta_both_summary_all.txt}"

R2_TEST="PEARSON"
R2_BASE="${R2_BASE:-/path/to/SI_GWAS/results/pgs_outputs}/"

# requested output base
PGS_OUT_BASE="${PGS_OUT_BASE:-${R2_BASE}pgs_out}/"

VALID_POPS=("white_euro" "afr" "asn")

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

# Pick phenotype by array index
pheno_num=$((SLURM_ARRAY_TASK_ID-1))
selected_pheno="${PHENOARRAY[$pheno_num]:-}"

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Selected phenotype: $selected_pheno"

if [[ -z "$selected_pheno" ]]; then
  echo "FATAL: no phenotype for task id $SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# Decide R2_DIR using lambda_WM (matches compile_test_pgs.R logic)
if [[ ! -s "$META_TAB" ]]; then
  echo "FATAL: meta table missing/empty: $META_TAB" >&2
  exit 1
fi

lambda_val=$(
  awk -F'\t' -v ph="$selected_pheno" '
    function norm(s){ s=tolower(s); gsub(/-/, "_", s); return s }
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="phenotype") p=i;
        if($i=="lambda_WM") l=i;
      }
      next
    }
    (norm($p) == norm(ph)){ print $l; exit }
  ' "$META_TAB"
)

if [[ -z "$lambda_val" ]]; then
  echo "FATAL: could not find lambda_WM for phenotype=${selected_pheno} in $META_TAB" >&2
  exit 1
fi
echo "lambda_WM for ${selected_pheno}: ${lambda_val}"

is_log=$(
  awk -v x="$lambda_val" 'BEGIN{
    if ((x+0)==0) print 1; else print 0
  }'
)

if [[ "$is_log" -eq 1 ]]; then
  if [[ "$R2_TEST" == "PEARSON" ]]; then
    R2_DIR="${R2_BASE}r2_out_log/"
  else
    R2_DIR="${R2_BASE}r2_out_spearman_log/"
  fi
else
  if [[ "$R2_TEST" == "PEARSON" ]]; then
    R2_DIR="${R2_BASE}r2_out/"
  else
    R2_DIR="${R2_BASE}r2_out_spearman/"
  fi
fi
echo "Using R2_DIR: $R2_DIR"

# Find joined GWAS file (*.glm.linear)
gwas_file=$(ls -1 "${JOINED_GWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}."*.glm.linear 2>/dev/null | head -n 1 || true)
if [[ -z "$gwas_file" ]]; then
  echo "FATAL: joined GWAS file not found for ${selected_pheno} in ${JOINED_GWAS}${selected_pheno}/" >&2
  exit 1
fi
echo "Using GWAS file: $gwas_file"

# Clumped SNP list (from test step)
valid_snp_file="${CLUMPED_GENOS}ukb_chr1-22_${TEST_POP}_${selected_pheno}.valid.snp"
if [[ ! -s "$valid_snp_file" ]]; then
  echo "FATAL: clumped SNP list missing/empty: $valid_snp_file" >&2
  exit 1
fi
echo "Using SNP list: $valid_snp_file"

# For each metric: std, spearman, rankLin, rankEC
for metric in std spearman rankLin rankEC; do
  processed_pgs="${R2_DIR}processed_${selected_pheno}_${TEST_POP}_PGS_${metric}.txt"
  if [[ ! -s "$processed_pgs" ]]; then
    echo "FATAL: missing processed PGS file: $processed_pgs" >&2
    exit 1
  fi

  best_pgs_thresh=$(awk 'NR==2{print $1; exit}' "$processed_pgs")
  if [[ -z "$best_pgs_thresh" ]]; then
    echo "FATAL: could not parse best threshold from: $processed_pgs" >&2
    exit 1
  fi
  echo "[${metric}] Best threshold: $best_pgs_thresh"

  # Single-threshold range list for plink2
  best_thresh_dir="${PGS_OUT_BASE}best_thresh_${metric}/"
  pgs_path="${PGS_OUT_BASE}valid_${metric}/"
  mkdir -p "$best_thresh_dir" "$pgs_path"

  pgs_valid_thresh="${best_thresh_dir}${selected_pheno}_thresh.txt"
  printf "%s 0 %s\n" "$best_pgs_thresh" "$best_pgs_thresh" > "$pgs_valid_thresh"
  echo "[${metric}] Wrote threshold file: $pgs_valid_thresh"

  # Score each validation population
  for valid_pop in "${VALID_POPS[@]}"; do
    pop_path="${GENOS}${valid_pop}/"
    pop_file="ukb_chr1-22_${valid_pop}"

    # choose bfile prefix robustly
    if [[ -f "${pop_path}QC/${pop_file}_QC.bed" ]]; then
      bfile="${pop_path}QC/${pop_file}_QC"
    elif [[ -f "${pop_path}${pop_file}_QC.bed" ]]; then
      bfile="${pop_path}${pop_file}_QC"
    else
      bfile="${pop_path}${pop_file}"
    fi

    if [[ ! -f "${bfile}.bed" ]]; then
      echo "FATAL: missing genotype bed for ${valid_pop}: ${bfile}.bed" >&2
      exit 1
    fi

    echo "[${metric}] Scoring ${selected_pheno} in ${valid_pop} using bfile=${bfile}"

    pgs_out_file="${pgs_path}${selected_pheno}_${valid_pop}_pgs_valid_${metric}"

    if [[ ! -f "${pgs_out_file}.${best_pgs_thresh}.sscore" ]]; then
      "${LOCAL_PLINK2}" \
        --bfile "${bfile}" \
        --keep "${bfile}.fam" \
        --score "${gwas_file}" 3 7 12 header cols=+scoresums \
        --q-score-range "${pgs_valid_thresh}" "${gwas_file}" 3 15 header \
        --extract "${valid_snp_file}" \
        --out "${pgs_out_file}"
    else
      echo "[${metric}] Exists: ${pgs_out_file}.${best_pgs_thresh}.sscore"
    fi
  done
done


echo "valid PGS done"
end_time=$(date +%s)
echo "Job Runtime: $((end_time - start_time)) seconds"
