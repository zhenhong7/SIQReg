#!/usr/bin/env bash
# PredInterval pipeline config. Update paths before running.

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"

GENO_MERGED="${GENO_MERGED:-${PROJECT_DIR}/genotypes/pop_genos/whitebrit/ukb_chr1-22_whitebrit}"

GENO_CHR_DIR="${GENO_CHR_DIR:-${PROJECT_DIR}/genotypes/whitebrit}"
GENO_CHR_PREFIX="ukb_chr"
GENO_CHR_SUFFIX="_whitebrit_QC"

PHENO_ORIGINAL_DIR="${PHENO_ORIGINAL_DIR:-${PROJECT_DIR}/extracted_phenotypes}"
PHENO_BC_DIR="${PHENO_BC_DIR:-${PROJECT_DIR}/SI_GWAS/intermediate_files}"

COVAR_FILE="${COVAR_FILE:-${PHENO_ORIGINAL_DIR}/covariates_sa40PC/covariates_sa40PC674178.pheno}"
COVAR_COLS="3-14"

OUTBASE="${OUTBASE:-${PROJECT_DIR}/SI_PredInterval}"

PLINK2="${PLINK2:-plink2}"

K=5
TRAIN_RATIO=0.8
P_THRESHOLD=0.05
CLUMP_R2=0.1
CLUMP_KB=250
CONF_LEVEL=0.95

declare -A TRAIT_INFO
TRAIT_INFO[1]="arm_fat-free_mass_left:Arm_fat-free_mass_left:674178"
TRAIT_INFO[2]="arm_fat-free_mass_right:Arm_fat-free_mass_right:674178"
TRAIT_INFO[3]="bmi:BMI:674178"
TRAIT_INFO[4]="birth_weight:Birth_weight:674178"
TRAIT_INFO[5]="c-reactive_protein:C-reactive_protein:674178"
TRAIT_INFO[6]="calcium:Calcium:674178"
TRAIT_INFO[7]="creatinine:Creatinine:674178"
TRAIT_INFO[8]="diastolicbp_auto:DiastolicBP_auto:674178"
TRAIT_INFO[9]="fev1:FEV1:674206"
TRAIT_INFO[10]="hdl:HDL:674178"
TRAIT_INFO[11]="hba1c:HbA1c:674178"
TRAIT_INFO[12]="height:Height:674178"
TRAIT_INFO[13]="hip_to_waist:Hip_to_waist:674178"
TRAIT_INFO[14]="igf-1:IGF-1:674178"
TRAIT_INFO[15]="ldl:LDL:674178"
TRAIT_INFO[16]="leukocyte_count:Leukocyte_count:674178"
TRAIT_INFO[17]="platelet_count:Platelet_count:674178"
TRAIT_INFO[18]="pulse_rate:Pulse_rate:674178"
TRAIT_INFO[19]="rbc:RBC:674178"
TRAIT_INFO[20]="shbg:SHBG:674178"
TRAIT_INFO[21]="systolicbp_auto:SystolicBP_auto:674178"
TRAIT_INFO[22]="testosterone:Testosterone:674178"
TRAIT_INFO[23]="urate:Urate:674178"
TRAIT_INFO[24]="urea:Urea:674178"
TRAIT_INFO[25]="vitamin_d:Vitamin_D:674178"

get_trait_info() {
    local idx=$1
    echo "${TRAIT_INFO[$idx]}"
}

export GENO_MERGED GENO_CHR_DIR GENO_CHR_PREFIX GENO_CHR_SUFFIX
export PHENO_ORIGINAL_DIR PHENO_BC_DIR COVAR_FILE COVAR_COLS
export OUTBASE PLINK2
export K TRAIN_RATIO P_THRESHOLD CLUMP_R2 CLUMP_KB CONF_LEVEL
