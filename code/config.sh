#!/usr/bin/env bash
# =============================================================================
# config.sh  —  Shared variables for all pipeline steps
# Source this file at the top of every step script:
#   source /home/ec2-user/code/config.sh
# =============================================================================

# --- Paths -------------------------------------------------------------------
export BASE_DIR="/home/ec2-user"
export CODE_DIR="${BASE_DIR}/code"
export REF_DIR="${BASE_DIR}/ref/mm10"
export TMP_DIR="/tmp/neoantig_pipeline"

# --- Reference ---------------------------------------------------------------
export REF="${REF_DIR}/mm10.fa"
export REF_DICT="${REF_DIR}/mm10.dict"
export REF_FAI="${REF_DIR}/mm10.fa.fai"
export EXON_BED="${REF_DIR}/mm10_exons.bed"
export KNOWN_SNPS="${REF_DIR}/mgp_snps.vcf.gz"
export SNPEFF_DB="GRCm38.86"

# --- Run ID (timestamped) ----------------------------------------------------
# Set once by run_pipeline.sh and exported; individual steps inherit it.
# If a step is run standalone, a new timestamp is generated.
export RUN_ID="${RUN_ID:-res_$(date '+%Y%m%d_%H%M%S')}"
export RESULTS_DIR="${BASE_DIR}/results/${RUN_ID}"

# --- S3 ----------------------------------------------------------------------
export S3_ROOT="s3://neoantigen2026-rerun"
export S3_BAM="${S3_ROOT}/data/bam/wes"
export S3_REF="${S3_ROOT}/data/reference/wes"
export S3_RESULTS="${S3_ROOT}/results/${RUN_ID}"

# --- Tools -------------------------------------------------------------------
export CONDA_BIN="/home/ec2-user/miniforge3/bin"
export SAMTOOLS="${CONDA_BIN}/samtools"
export BCFTOOLS="${CONDA_BIN}/bcftools"
export GATK_LOCAL_JAR="${CONDA_BIN}/../share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"
export GATK="${CONDA_BIN}/gatk"
export PYTHON="${CONDA_BIN}/python3"
export PICARD="${CONDA_BIN}/picard"
export SNPEFF="${CONDA_BIN}/snpEff"
export NETMHCPAN="${CONDA_BIN}/netMHCpan"

# --- Samples -----------------------------------------------------------------
export NORMAL_SAMPLE="423_D0_old"
export NORMAL_SM="D0_old"
export PON_NORMALS=("423_D0_old" "424_D0__old")

# Tumour samples — ordered by timepoint
export TUMOR_SAMPLES=("443_D21_new" "428_D20_new" "34_D52_old" "36_D99_new" "38_D99_new" "42_D122_old")
declare -A TUMOR_SM_TAGS
TUMOR_SM_TAGS["443_D21_new"]="D21_new"
TUMOR_SM_TAGS["428_D20_new"]="D20_new"
TUMOR_SM_TAGS["34_D52_old"]="D52_old"
TUMOR_SM_TAGS["36_D99_new"]="D99_new"
TUMOR_SM_TAGS["38_D99_new"]="D99_new"
TUMOR_SM_TAGS["42_D122_old"]="D122_old"
export TUMOR_SM_TAGS

# MHC alleles (C57BL/6 H-2b haplotype)
export MHC_ALLELES="H-2-Kb,H-2-Db"
export PEPTIDE_LENGTHS="8,9,10"

# --- Runtime -----------------------------------------------------------------
export THREADS=4
export JAVA_OPTS="-Xmx12g"

# --- Utility functions -------------------------------------------------------
log()  { echo "[$(date '+%H:%M:%S')] $*"; }
skip() { log "[SKIP] $1 already exists"; }
s3up() {
    local src="$1" dst="${S3_RESULTS}/$2"
    [[ -f "${src}" ]] && aws s3 cp "${src}" "${dst}" && log "Uploaded: ${dst}"
}
s3rm_local() { [[ -f "$1" ]] && rm -f "$1" && log "Removed local: $1"; }

mkdir -p "${TMP_DIR}" "${RESULTS_DIR}"
log "Run ID: ${RUN_ID}  →  local: ${RESULTS_DIR}  |  S3: ${S3_RESULTS}"
