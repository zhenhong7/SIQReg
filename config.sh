#!/usr/bin/env bash
# Central configuration: update paths below to match your environment.

PROJECT_ROOT="/path/to/SIQReg"

GENO_DIR="${PROJECT_ROOT}/genotypes"
GENO_CHR_DIR="${GENO_DIR}/per_chr"               # per-chromosome PLINK files
GENO_MERGED="${GENO_DIR}/merged/ukb_merged_QC"   # merged genotype prefix

PHENO_DIR="${PROJECT_ROOT}/phenotypes"
PHENO_ORIGINAL_DIR="${PHENO_DIR}/original"
PHENO_BC_DIR="${PHENO_DIR}/box_cox"
COVAR_FILE="${PHENO_DIR}/covariates_sa40PC.pheno"

ANCESTRY_DIR="${PROJECT_ROOT}/ancestry_ids"
# Expected files: whitebrit.fid_iid.txt, white_euro.fid_iid.txt,
#                 afr.fid_iid.txt, asn.fid_iid.txt

PRS_DIR="${PROJECT_ROOT}/PRS/ALL_PRS"
TOP1_PRS_CSV="${PROJECT_ROOT}/PRS/R2/top1_prs_by_pheno.csv"

META_BOTH_FILE="${PROJECT_ROOT}/Lambda_learning/meta_both_summary_all.txt"

EXPR_DIR="${PROJECT_ROOT}/TWAS/imputed_expr"
PHENO2TISSUE="${PROJECT_ROOT}/TWAS/pheno2tissue.txt"
TWAS_MAPPING_DIR="${PROJECT_ROOT}/TWAS/mapping"

PRS_SCORING_DIR="${PROJECT_ROOT}/PRS/polygenic_scores"
COVAR_FULL_FILE="${PHENO_DIR}/covar_full/covar_full_age2.pheno"
STATIN_FILE="${PHENO_DIR}/Statins/Statins.pheno"

BIRTH_YEAR_BASELINE=2025

TRAIN_POP="whitebrit"

HLMM_BASE="${PROJECT_ROOT}/HLMM"
HLMM_COV="${HLMM_BASE}/inputs/covar_3_14.fam"
HLMM_PY="${HLMM_BASE}/code/hlmm/bin/hlmm_chr.py"
HLMM_FAST_PY="${HLMM_BASE}/code/hlmm_fast/bin/hlmm_chr.py"
HLMM_PHENO_DIR="${HLMM_BASE}/pheno"
HLMM_PHENO_BC_DIR="${HLMM_BASE}/pheno_bc"

GXE_DIR="${RESULTS_ROOT}/GxE"
GXE_META_FILE="${GXE_DIR}/meta_both_summary_all.txt"

PLINK2="plink2"   # path to PLINK2 binary

RESULTS_ROOT="${PROJECT_ROOT}/results"
GWAS_RESULTS="${RESULTS_ROOT}/GWAS"
TWAS_RESULTS="${RESULTS_ROOT}/TWAS"
HLMM_RESULTS="${RESULTS_ROOT}/HLMM"
PREDINTERVAL_OUTBASE="${RESULTS_ROOT}/PredInterval"

SCRATCH_DIR="${SCRATCH:-/tmp}/${USER}/SIQReg"

get_pop_geno_dir() {
  local pop="$1"
  echo "${GENO_DIR}/${pop}"
}

get_pop_geno_prefix() {
  local pop="$1" chr="$2"
  echo "${GENO_DIR}/${pop}/ukb_chr${chr}_${pop}_QC"
}

PHENOS=(
  Arm_fat-free_mass_left
  Arm_fat-free_mass_right
  BMI
  Birth_weight
  C-reactive_protein
  Calcium
  Creatinine
  DiastolicBP_auto
  FEV1
  HDL
  HbA1c
  Height
  Hip_to_waist
  IGF-1
  LDL
  Leukocyte_count
  Platelet_count
  Pulse_rate
  RBC
  SHBG
  SystolicBP_auto
  Testosterone
  Urate
  Urea
  Vitamin_D
)

COVARS=(sex age age2 $(printf "PC%d " {1..40}))
