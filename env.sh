#!/usr/bin/env bash
# =============================================================================
# env.sh — NeoAntigen 2026 WES Pipeline Environment
#
# Copy this file and source it before running the pipeline:
#   cp env.sh ~/pipeline_env.sh
#   source ~/pipeline_env.sh
#
# All variables here can override the defaults in code/config.sh.
# Only edit what differs from your setup.
# =============================================================================

# -----------------------------------------------------------------------------
# Directory layout
# Adjust BASE_DIR if your home directory or mount point differs.
# -----------------------------------------------------------------------------
export BASE_DIR="/home/ec2-user"
export CODE_DIR="${BASE_DIR}/code"
export REF_DIR="${BASE_DIR}/ref/mm10"
export TMP_DIR="/tmp/neoantig_pipeline"

# -----------------------------------------------------------------------------
# Conda / tool binaries
# Point to your Miniforge/Miniconda/Anaconda installation.
# -----------------------------------------------------------------------------
export CONDA_BIN="${BASE_DIR}/miniforge3/bin"
export SAMTOOLS="${CONDA_BIN}/samtools"
export BCFTOOLS="${CONDA_BIN}/bcftools"
export GATK="${CONDA_BIN}/gatk"
export GATK_LOCAL_JAR="${CONDA_BIN}/../share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"
export PYTHON="${CONDA_BIN}/python3"
export PICARD="${CONDA_BIN}/picard"
export SNPEFF="${CONDA_BIN}/snpEff"
export NETMHCPAN="${CONDA_BIN}/netMHCpan"

# Add conda bin to PATH
export PATH="${CONDA_BIN}:${PATH}"

# -----------------------------------------------------------------------------
# S3 storage
# Update S3_ROOT if using a different bucket.
# -----------------------------------------------------------------------------
export S3_ROOT="s3://neoantigen2026-rerun"
export S3_BAM="${S3_ROOT}/data/bam/wes"
export S3_REF="${S3_ROOT}/data/reference/wes"
# S3_RESULTS is set automatically per run — do not set here

# -----------------------------------------------------------------------------
# AWS region (set if not already configured via aws configure)
# -----------------------------------------------------------------------------
export AWS_DEFAULT_REGION="us-east-1"

# -----------------------------------------------------------------------------
# Runtime resources
# Adjust THREADS to match your instance vCPU count.
# Adjust JAVA_OPTS heap based on available RAM (leave ~4 GB for OS).
#   r5.large  (2 vCPU,  16 GB):  THREADS=2, JAVA_OPTS="-Xmx10g"
#   r5.xlarge (4 vCPU,  32 GB):  THREADS=4, JAVA_OPTS="-Xmx26g"
#   r5.2xlarge(8 vCPU,  64 GB):  THREADS=8, JAVA_OPTS="-Xmx56g"
# -----------------------------------------------------------------------------
export THREADS=4
export JAVA_OPTS="-Xmx12g"

# -----------------------------------------------------------------------------
# Reference genome
# -----------------------------------------------------------------------------
export REF="${REF_DIR}/mm10.fa"
export REF_DICT="${REF_DIR}/mm10.dict"
export REF_FAI="${REF_DIR}/mm10.fa.fai"
export EXON_BED="${REF_DIR}/mm10_exons.bed"
export KNOWN_SNPS="${REF_DIR}/mgp_snps.vcf.gz"
export SNPEFF_DB="GRCm38.86"

# -----------------------------------------------------------------------------
# Samples
# Modify if running a different set of samples.
# -----------------------------------------------------------------------------
export NORMAL_SAMPLE="423_D0_old"
export NORMAL_SM="D0_old"
export PON_NORMALS=("423_D0_old" "424_D0__old")

export TUMOR_SAMPLES=("428_D20_new" "34_D52_old" "36_D99_new" "38_D99_new" "42_D122_old")

declare -A TUMOR_SM_TAGS
TUMOR_SM_TAGS["428_D20_new"]="D20_new"
TUMOR_SM_TAGS["34_D52_old"]="D52_old"
TUMOR_SM_TAGS["36_D99_new"]="D99_new"
TUMOR_SM_TAGS["38_D99_new"]="D99_new"
TUMOR_SM_TAGS["42_D122_old"]="D122_old"
export TUMOR_SM_TAGS

# -----------------------------------------------------------------------------
# MHC alleles — C57BL/6 H-2b haplotype
# Change if using a different mouse strain.
# -----------------------------------------------------------------------------
export MHC_ALLELES="H-2-Kb,H-2-Db"
export PEPTIDE_LENGTHS="8,9,10"

# -----------------------------------------------------------------------------
echo "[env.sh] Environment loaded."
echo "  BASE_DIR   : ${BASE_DIR}"
echo "  CONDA_BIN  : ${CONDA_BIN}"
echo "  S3_ROOT    : ${S3_ROOT}"
echo "  THREADS    : ${THREADS}"
echo "  JAVA_OPTS  : ${JAVA_OPTS}"
echo "  NORMAL     : ${NORMAL_SAMPLE}"
echo "  TUMORS     : ${TUMOR_SAMPLES[*]}"
