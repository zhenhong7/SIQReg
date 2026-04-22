#!/usr/bin/env bash
# Cross-ancestry PredInterval config. Update paths before running.

PROJECT_DIR="${PROJECT_DIR:-/path/to/project}"
WB_RESULTS="${WB_RESULTS:-${PROJECT_DIR}/SI_PredInterval}"

declare -A ANCESTRY_LABEL
ANCESTRY_LABEL[asn]="East Asian"
ANCESTRY_LABEL[afr]="African"
ANCESTRY_LABEL[white_euro]="White European (non-British)"

declare -A ANCESTRY_IDS
ANCESTRY_IDS[asn]="${ANCESTRY_IDS_DIR:-${PROJECT_DIR}/ancestry_ids}/asn.fid_iid.txt"
ANCESTRY_IDS[afr]="${ANCESTRY_IDS_DIR:-${PROJECT_DIR}/ancestry_ids}/afr.fid_iid.txt"
ANCESTRY_IDS[white_euro]="${ANCESTRY_IDS_DIR:-${PROJECT_DIR}/ancestry_ids}/white_euro.fid_iid.txt"

declare -A GENO_DIRS
GENO_DIRS[asn]="${GENO_BASE:-${PROJECT_DIR}/genotypes}/asn"
GENO_DIRS[afr]="${GENO_BASE:-${PROJECT_DIR}/genotypes}/afr"
GENO_DIRS[white_euro]="${GENO_BASE:-${PROJECT_DIR}/genotypes}/white_euro"

declare -A GENO_SUFFIX
GENO_SUFFIX[asn]="asn"
GENO_SUFFIX[afr]="afr"
GENO_SUFFIX[white_euro]="white_euro"

declare -A PHENO_BC_DIRS
PHENO_BC_DIRS[asn]="${PHENO_BC_BASE:-${PROJECT_DIR}/SI_GWAS}/asn/intermediate_files"
PHENO_BC_DIRS[afr]="${PHENO_BC_BASE:-${PROJECT_DIR}/SI_GWAS}/afr/intermediate_files"
PHENO_BC_DIRS[white_euro]="${PHENO_BC_BASE:-${PROJECT_DIR}/SI_GWAS}/white_euro/intermediate_files"

PHENO_ORIGINAL_DIR="${PHENO_ORIGINAL_DIR:-${PROJECT_DIR}/extracted_phenotypes}"

declare -A CALIB_RATIO
CALIB_RATIO[asn]=0.5
CALIB_RATIO[afr]=0.7
CALIB_RATIO[white_euro]=0.7

OUTBASE="${OUTBASE:-${PROJECT_DIR}/SI_PredInterval/cross_ancestry}"
PLINK2="${PLINK2:-plink2}"
K=5
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

export WB_RESULTS PHENO_ORIGINAL_DIR OUTBASE PLINK2 K CONF_LEVEL
