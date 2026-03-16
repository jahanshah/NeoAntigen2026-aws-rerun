#!/usr/bin/env bash
# Wait for:
#   1. Current Mutect2 run (PID 7058) to finish
#   2. Kill main pipeline (PID 7017) before it proceeds to steps 3-12 with
#      incomplete data (only 443/428) — PyClone requires all 6 tumours.
#   3. All 3 missing preprocessed BAMs to appear on S3
#   4. Re-run Mutect2 with boosted resources for the 4 remaining tumours
#   5. Run steps 3-12 over all 6 tumours

set -euo pipefail

CURRENT_MUTECT2_PID=7058
MAIN_PIPELINE_PID=7017
MISSING_SAMPLES=(36_D99_new 38_D99_new 42_D122_old)
S3_PREPROC="s3://neoantigen2026-rerun/data/bam/wes/preprocessed"
LOG="/home/ec2-user/results/res_20260311_225555/mutect2_rerun.log"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "${LOG}"; }

mkdir -p "$(dirname "${LOG}")"
log "=== Mutect2 re-run watcher started ==="
log "  Waiting for PID ${CURRENT_MUTECT2_PID} + BAMs: ${MISSING_SAMPLES[*]}"

# ── 1. Wait for current Mutect2 to finish ───────────────────────────────────
log "Waiting for current Mutect2 (PID ${CURRENT_MUTECT2_PID}) to finish..."
while kill -0 "${CURRENT_MUTECT2_PID}" 2>/dev/null; do
    sleep 60
done
log "Current Mutect2 (443/428) done."

# ── 2. Kill main pipeline before it runs steps 3-12 with incomplete data ────
# PyClone (step 6) is a joint multi-sample analysis — must have all 6 tumours.
if kill -0 "${MAIN_PIPELINE_PID}" 2>/dev/null; then
    log "Killing main pipeline (PID ${MAIN_PIPELINE_PID}) to prevent premature steps 3-12..."
    kill "${MAIN_PIPELINE_PID}" 2>/dev/null || true
    sleep 5
    kill -9 "${MAIN_PIPELINE_PID}" 2>/dev/null || true
    log "Main pipeline stopped."
else
    log "[WARN] Main pipeline (PID ${MAIN_PIPELINE_PID}) already exited — steps 3-12 may have started prematurely."
fi

# ── 3. Wait for all 3 missing BAMs on S3 ────────────────────────────────────
log "Waiting for preprocessed BAMs on S3..."
while true; do
    ALL_READY=1
    for SAMPLE in "${MISSING_SAMPLES[@]}"; do
        if ! aws s3 ls "${S3_PREPROC}/${SAMPLE}.preproc.bam" &>/dev/null; then
            log "  [WAIT] ${SAMPLE} not yet on S3..."
            ALL_READY=0
        else
            log "  [OK]   ${SAMPLE} ready on S3"
        fi
    done
    [[ "${ALL_READY}" -eq 1 ]] && break
    sleep 300  # check every 5 min
done
log "All 3 BAMs on S3. Launching Mutect2 with boosted resources."

# ── 4. Run Mutect2 with extra resources for 4 remaining tumours ─────────────
export PATH=/home/ec2-user/miniforge3/bin:$PATH
export RUN_ID=res_20260311_225555
export THREADS=16
export JAVA_OPTS="-Xmx100g"
export MUTECT2_PARALLEL=2

log "Config: THREADS=${THREADS}, JAVA_OPTS=${JAVA_OPTS}, MUTECT2_PARALLEL=${MUTECT2_PARALLEL}"
log "(443/428 already on S3 — will be auto-skipped)"

bash /home/ec2-user/code/steps/02_mutect2.sh 2>&1 | tee -a "${LOG}"
log "=== Mutect2 re-run complete. All 6 tumours done. ==="

# ── 5. Run steps 3-12 over all 6 tumours ────────────────────────────────────
log "Starting steps 3-12 (full pipeline with all 6 tumours)..."
cd /home/ec2-user
bash /home/ec2-user/NeoAntigen2026-aws-rerun/code/run_pipeline.sh 3 12 \
    2>&1 | tee -a "${LOG}"

log "=== WES pipeline complete ==="
