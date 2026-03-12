#!/usr/bin/env bash
# =============================================================================
# config_rnaseq.sh  —  RNASeq-specific configuration
# Sources the main config.sh first, then adds RNASeq-specific variables.
# Source this at the top of every RNASeq step script:
#   source /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh
# =============================================================================

source /home/ec2-user/NeoAntigen2026-aws-rerun/code/config.sh

# --- RNASeq directories ------------------------------------------------------
export RNASEQ_DIR="${CODE_DIR}/rnaseq"
export RNASEQ_STEPS_DIR="${RNASEQ_DIR}/steps"

# --- Reference files ---------------------------------------------------------
export STAR_INDEX="/home/ec2-user/ref/mm10/star_index"
export GTF="/home/ec2-user/ref/mm10/mm10.gtf"

# --- S3 paths ----------------------------------------------------------------
export S3_FASTQ="${S3_ROOT}/data/fastq/rnaseq"
export S3_RNASEQ_OUT="${S3_RESULTS}/rnaseq"

# --- Threading ---------------------------------------------------------------
export RNASEQ_THREADS=16
# Number of STAR alignment jobs to run concurrently (2 × 16-thread = 32 threads
# in flight; safe on a 16-vCPU machine because STAR I/O-waits frequently).
export STAR_ALIGN_PARALLEL=2

# --- Local paths -------------------------------------------------------------
export RNASEQ_RESULTS="${RESULTS_DIR}/rnaseq"
export RNASEQ_TMP="${TMP_DIR}/rnaseq"

# --- Manifest ----------------------------------------------------------------
export RNASEQ_META="${RNASEQ_DIR}/samples_rnaseq.tsv"

# --- Parse sample manifest into arrays ---------------------------------------
# Arrays populated from samples_rnaseq.tsv (skips comment/blank lines).
#
#   RNASEQ_SAMPLES   — ordered list of sample_id values
#   RNASEQ_WES_IDS   — associative: sample_id → wes_id   (NA for extras)
#   RNASEQ_BATCHES   — associative: sample_id → batch
#   RNASEQ_ROLES     — associative: sample_id → role
#   RNASEQ_PATTERNS  — associative: sample_id → fastq_pattern (regex)

RNASEQ_SAMPLES=()
declare -A RNASEQ_WES_IDS
declare -A RNASEQ_BATCHES
declare -A RNASEQ_ROLES
declare -A RNASEQ_STAGES
declare -A RNASEQ_PATTERNS

while IFS=$'\t' read -r _sid _wes _tp _batch _role _stage _pat; do
    # Skip comment lines and blank lines
    [[ "${_sid}" =~ ^[[:space:]]*# || -z "${_sid}" ]] && continue
    RNASEQ_SAMPLES+=("${_sid}")
    RNASEQ_WES_IDS["${_sid}"]="${_wes}"
    RNASEQ_BATCHES["${_sid}"]="${_batch}"
    RNASEQ_ROLES["${_sid}"]="${_role}"
    RNASEQ_STAGES["${_sid}"]="${_stage}"
    RNASEQ_PATTERNS["${_sid}"]="${_pat}"
done < "${RNASEQ_META}"

export RNASEQ_SAMPLES RNASEQ_WES_IDS RNASEQ_BATCHES RNASEQ_ROLES RNASEQ_STAGES RNASEQ_PATTERNS

# --- Ensure local output dirs exist ------------------------------------------
mkdir -p "${RNASEQ_RESULTS}" "${RNASEQ_TMP}"

log "RNASeq config loaded — ${#RNASEQ_SAMPLES[@]} samples | STAR_INDEX=${STAR_INDEX}"
