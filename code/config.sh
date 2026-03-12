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
# Use /scratch if available (500GB EBS scratch volume — see setup_scratch_volume.sh)
# Fall back to /tmp (15GB tmpfs) — safe for sequential runs, too small for parallel
export TMP_DIR="${TMP_DIR:-$( [[ -d /scratch ]] && echo /scratch/neoantig_pipeline || echo /home/ec2-user/tmp/neoantig_pipeline )}"

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
export S3_BAM_PREPROC="${S3_ROOT}/data/bam/wes/preprocessed"   # permanent preprocessed BAMs (reused across runs)
export S3_REF="${S3_ROOT}/data/reference/wes"
export S3_RESULTS="${S3_ROOT}/results/${RUN_ID}"

# --- Tools -------------------------------------------------------------------
export CONDA_BIN="/home/ec2-user/miniforge3/bin"
export SAMTOOLS="${CONDA_BIN}/samtools"
export BCFTOOLS="${CONDA_BIN}/bcftools"
export GATK_LOCAL_JAR="${CONDA_BIN}/../share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"
export GATK="${CONDA_BIN}/gatk"
export GATK_SPARK="GATK_LOCAL_JAR=${GATK_LOCAL_JAR} ${CONDA_BIN}/gatk"   # enables Spark tools
export PYTHON="${CONDA_BIN}/python3"
export PICARD="${CONDA_BIN}/picard"
export SNPEFF="${CONDA_BIN}/snpEff"
export NETMHCPAN="${CONDA_BIN}/netMHCpan"

# --- Samples (loaded from manifest) -----------------------------------------
# Edit code/samples.tsv to add/remove/rename samples.
# Columns: sample_id  role  sm_tag  timepoint  library
# Roles (comma-separated): normal | pon | tumor
export SAMPLES_META="${CODE_DIR}/samples.tsv"

TUMOR_SAMPLES=()
PON_NORMALS=()
declare -A TUMOR_SM_TAGS
declare -A PON_SM_TAGS

while IFS=$'\t' read -r _sid _role _sm _tp _lib; do
    # skip comment and blank lines
    [[ "${_sid}" =~ ^[[:space:]]*# || -z "${_sid}" ]] && continue
    if [[ "${_role}" == *"normal"* ]]; then
        export NORMAL_SAMPLE="${_sid}"
        export NORMAL_SM="${_sm}"
    fi
    if [[ "${_role}" == *"pon"* ]]; then
        PON_NORMALS+=("${_sid}")
        PON_SM_TAGS["${_sid}"]="${_sm}"
    fi
    if [[ "${_role}" == *"tumor"* ]]; then
        TUMOR_SAMPLES+=("${_sid}")
        TUMOR_SM_TAGS["${_sid}"]="${_sm}"
    fi
done < "${SAMPLES_META}"

export PON_NORMALS PON_SM_TAGS TUMOR_SAMPLES TUMOR_SM_TAGS

# MHC alleles (C57BL/6 H-2b haplotype)
export MHC_ALLELES="H-2-Kb,H-2-Db"
export PEPTIDE_LENGTHS="8,9,10"

# --- Runtime -----------------------------------------------------------------
# Instance sizing guide:
#   r5.xlarge  (4 vCPU,  32GB):  THREADS=3,  JAVA_OPTS="-Xmx12g"
#   r5.2xlarge (8 vCPU,  64GB):  THREADS=6,  JAVA_OPTS="-Xmx24g"
#   r5.4xlarge (16 vCPU, 128GB): THREADS=12, JAVA_OPTS="-Xmx40g"
#   Current:   16 vCPU, 30GB EBS root + /scratch EBS volume
export THREADS=12
export JAVA_OPTS="-Xmx20g"

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
