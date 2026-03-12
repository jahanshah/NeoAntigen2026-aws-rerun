#!/usr/bin/env bash
# =============================================================================
# run_rnaseq_pipeline.sh — RNASeq Pipeline (standalone)
#
# Independent of the WES pipeline. Runs STAR alignment, featureCounts
# quantification, DESeq2 normalisation, DEA, and limma batch correction.
# Outputs expression_matrix.tsv to ${RESULTS_DIR}/rnaseq/ and S3, where
# the WES step 12 (10_filter_results.R) picks it up for final ranking.
#
# Usage:
#   bash run_rnaseq_pipeline.sh           # full pipeline R1→R4
#   bash run_rnaseq_pipeline.sh R2        # start from alignment (STAR index exists)
#   bash run_rnaseq_pipeline.sh R2 R3     # alignment + quantification only
#   RUN_ID=res_20260311_225555 bash run_rnaseq_pipeline.sh   # tie to WES run ID
#
# Logs: ${RESULTS_DIR}/rnaseq_pipeline.log
# =============================================================================

set -euo pipefail

export RUN_ID="${RUN_ID:-res_$(date '+%Y%m%d_%H%M%S')}"
source /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh

# Normalise start/end arguments: accept "R2" or "2"
_normalise_step() {
    local s="${1#R}"  # strip leading "R" if present
    echo "${s}"
}

START_STEP="${1:-R1}"
END_STEP="${2:-R4}"
START_NUM=$(_normalise_step "${START_STEP}")
END_NUM=$(_normalise_step "${END_STEP}")

# Pipeline log file (append mode so re-runs accumulate)
PIPELINE_LOG="${RESULTS_DIR}/rnaseq_pipeline.log"
mkdir -p "${RESULTS_DIR}"
exec > >(tee -a "${PIPELINE_LOG}") 2>&1

log "============================================================"
log "RNASeq Pipeline  RUN_ID=${RUN_ID}"
log "Steps  : R${START_NUM} → R${END_NUM}"
log "Samples: ${RNASEQ_SAMPLES[*]}"
log "Local  : ${RNASEQ_RESULTS}"
log "S3     : ${S3_RNASEQ_OUT}"
log "Log    : ${PIPELINE_LOG}"
log "============================================================"

# =============================================================================
# Helper: run_step
# =============================================================================
run_step() {
    local STEP_TAG="$1"   # e.g. "R1", "R2"
    local SCRIPT="$2"     # full path to script
    local DESC="$3"

    local STEP_NUM
    STEP_NUM=$(_normalise_step "${STEP_TAG}")

    # Skip if outside requested range
    if [[ ${STEP_NUM} -lt ${START_NUM} || ${STEP_NUM} -gt ${END_NUM} ]]; then
        return 0
    fi

    log ""
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Step ${STEP_TAG}: ${DESC}"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    local T0=$SECONDS

    local EXT="${SCRIPT##*.}"
    case "${EXT}" in
        sh)  bash    "${SCRIPT}" ;;
        R)   Rscript "${SCRIPT}" ;;
        py)  python3 "${SCRIPT}" ;;
        *)   log "[ERROR] Unknown extension: ${EXT}"; exit 1 ;;
    esac

    log "Step ${STEP_TAG} done in $(( SECONDS - T0 ))s."
}

# =============================================================================
# RNASeq Steps R1 → R4
# =============================================================================
run_step R1 "${RNASEQ_STEPS_DIR}/R01_setup_ref.sh"   "Setup: install tools, GTF, STAR index"
run_step R2 "${RNASEQ_STEPS_DIR}/R02_align.sh"        "Align FASTQs with STAR (2-pass)"
run_step R3 "${RNASEQ_STEPS_DIR}/R03_quantify.sh"     "Quantify with featureCounts"
run_step R4 "${RNASEQ_STEPS_DIR}/R04_expression.R"    "DESeq2 normalisation + DEA + limma batch correction"

log ""
log "============================================================"
log "RNASeq pipeline complete."
log "  Expression matrix : ${RNASEQ_RESULTS}/expression_matrix.tsv"
log "  DEA results       : ${RNASEQ_RESULTS}/dea_results.tsv"
log "  S3                : ${S3_RNASEQ_OUT}"
log "  Log               : ${PIPELINE_LOG}"
log ""
log "  To run final neoantigen ranking once WES is also done:"
log "    RUN_ID=${RUN_ID} Rscript /home/ec2-user/code/steps/10_filter_results.R"
log "============================================================"
