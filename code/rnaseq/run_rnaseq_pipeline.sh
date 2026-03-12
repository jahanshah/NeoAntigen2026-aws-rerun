#!/usr/bin/env bash
# =============================================================================
# run_rnaseq_pipeline.sh — RNASeq Pipeline Orchestrator
#
# Runs in parallel with the WES pipeline (run_pipeline.sh).
# Steps R1–R4 run sequentially; after R4 completes, the script polls for
# the WES netMHCpan output and then runs the final integration (step 12).
#
# Usage:
#   bash run_rnaseq_pipeline.sh           # full RNA pipeline (R1→R4 + integration)
#   bash run_rnaseq_pipeline.sh R2        # start from alignment (skip R1)
#   bash run_rnaseq_pipeline.sh R2 R3     # run only R2 and R3
#   RUN_ID=res_20260312_120000 bash run_rnaseq_pipeline.sh  # use specific run ID
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
run_step R4 "${RNASEQ_STEPS_DIR}/R04_expression.R"    "DESeq2 + VST + limma batch correction"

# =============================================================================
# Integration: poll for WES netMHCpan output, then run step 12
# (Only when the full pipeline was requested, i.e. END_NUM >= 4)
# =============================================================================
if [[ ${END_NUM} -ge 4 && ${START_NUM} -le 4 ]]; then
    log ""
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Waiting for WES netMHCpan output before running Step 12..."
    log "  Polling: ${S3_RESULTS}/netmhcpan/all_predictions.tsv"
    log "  Timeout: 8 hours | Interval: 10 minutes"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    S3_NETMHC="${S3_RESULTS}/netmhcpan/all_predictions.tsv"
    EXPR_MATRIX="${RNASEQ_RESULTS}/expression_matrix.tsv"
    MAX_WAIT=28800    # 8 hours in seconds
    POLL_INTERVAL=600 # 10 minutes
    WAITED=0

    netmhc_ready=0
    while [[ ${WAITED} -lt ${MAX_WAIT} ]]; do
        if aws s3 ls "${S3_NETMHC}" &>/dev/null; then
            log "netMHCpan output found on S3: ${S3_NETMHC}"
            netmhc_ready=1
            break
        fi
        log "netMHCpan not ready yet — waiting ${POLL_INTERVAL}s (${WAITED}/${MAX_WAIT}s elapsed)..."
        sleep "${POLL_INTERVAL}"
        WAITED=$(( WAITED + POLL_INTERVAL ))
    done

    if [[ ${netmhc_ready} -eq 0 ]]; then
        log "[WARN] Timeout waiting for netMHCpan output after 8 hours."
        log "[WARN] Step 12 will run with whatever is available."
    fi

    # Verify expression_matrix.tsv exists (from R04 above)
    if [[ ! -f "${EXPR_MATRIX}" ]]; then
        log "[WARN] expression_matrix.tsv not found locally — attempting S3 fetch..."
        aws s3 cp "${S3_RNASEQ_OUT}/expression_matrix.tsv" "${EXPR_MATRIX}" || \
            log "[WARN] expression_matrix.tsv not on S3 either — step 12 will run without expression data"
    fi

    log ""
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Step 12: Neoantigen candidate ranking (10_filter_results.R)"
    log "  Expression data available: $([[ -f "${EXPR_MATRIX}" ]] && echo YES || echo NO)"
    log "  netMHCpan data available : $([[ ${netmhc_ready} -eq 1 ]] && echo YES || echo NO)"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    T0=$SECONDS
    Rscript /home/ec2-user/NeoAntigen2026-aws-rerun/code/steps/10_filter_results.R
    log "Step 12 done in $(( SECONDS - T0 ))s."
fi

log ""
log "============================================================"
log "RNASeq pipeline complete."
log "  Local  : ${RNASEQ_RESULTS}"
log "  S3     : ${S3_RNASEQ_OUT}"
log "  Log    : ${PIPELINE_LOG}"
log "============================================================"
