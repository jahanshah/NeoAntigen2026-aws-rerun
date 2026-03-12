#!/usr/bin/env bash
# =============================================================================
# instance-resize.sh — Stop, resize, restart EC2 and re-run pipelines
#
# Run this from your LOCAL Mac, not from the EC2 instance.
#
# Usage:
#   bash instance-resize.sh                   # resize to default (r5.4xlarge)
#   bash instance-resize.sh r5.8xlarge        # specify instance type
#
# What it does:
#   1. Stops the EC2 instance
#   2. Resizes to the target instance type
#   3. Starts the instance and waits for SSH
#   4. Updates pipeline configs (THREADS, JAVA_OPTS, parallelism)
#   5. Kills any leftover pipeline processes
#   6. Launches WES (step 2 onwards) + RNAseq in separate screen sessions
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION — edit these if needed
# =============================================================================
INSTANCE_ID="i-01a99e6e02659e1e7"
REGION="ca-central-1"
SSH_KEY="${SSH_KEY:-~/.ssh/id_rsa}"          # override with: SSH_KEY=~/.ssh/mykey.pem bash instance-resize.sh
SSH_USER="ec2-user"
TARGET_TYPE="${1:-r5.4xlarge}"
RUN_ID="res_20260311_225555"                 # existing run to resume
# =============================================================================

# Per-instance-type tuning
# Format: THREADS  JAVA_XMX  MUTECT2_PARALLEL  STAR_ALIGN_PARALLEL
# Mutect2 peak RAM = MUTECT2_PARALLEL × JAVA_XMX
# RNAseq peak RAM  = STAR_ALIGN_PARALLEL × ~28GB (genome per process)
declare -A CFG_THREADS       CFG_JAVA_XMX       CFG_MUTECT_PAR    CFG_STAR_PAR
CFG_THREADS[r5.2xlarge]=4;   CFG_JAVA_XMX[r5.2xlarge]=12;  CFG_MUTECT_PAR[r5.2xlarge]=3;  CFG_STAR_PAR[r5.2xlarge]=1
CFG_THREADS[r5.4xlarge]=4;   CFG_JAVA_XMX[r5.4xlarge]=16;  CFG_MUTECT_PAR[r5.4xlarge]=6;  CFG_STAR_PAR[r5.4xlarge]=2
CFG_THREADS[r5.8xlarge]=8;   CFG_JAVA_XMX[r5.8xlarge]=20;  CFG_MUTECT_PAR[r5.8xlarge]=6;  CFG_STAR_PAR[r5.8xlarge]=4
CFG_THREADS[r5.16xlarge]=10; CFG_JAVA_XMX[r5.16xlarge]=28; CFG_MUTECT_PAR[r5.16xlarge]=6; CFG_STAR_PAR[r5.16xlarge]=6
CFG_THREADS[r4.8xlarge]=8;   CFG_JAVA_XMX[r4.8xlarge]=20;  CFG_MUTECT_PAR[r4.8xlarge]=6;  CFG_STAR_PAR[r4.8xlarge]=4

if [[ -z "${CFG_THREADS[$TARGET_TYPE]+x}" ]]; then
    echo "ERROR: Unknown instance type '${TARGET_TYPE}'. Add it to the CFG_ maps in this script."
    exit 1
fi

THREADS="${CFG_THREADS[$TARGET_TYPE]}"
JAVA_XMX="${CFG_JAVA_XMX[$TARGET_TYPE]}"
MUTECT_PAR="${CFG_MUTECT_PAR[$TARGET_TYPE]}"
STAR_PAR="${CFG_STAR_PAR[$TARGET_TYPE]}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

log "============================================================"
log "Target instance : ${TARGET_TYPE}"
log "Mutect2         : ${MUTECT_PAR} parallel × -Xmx${JAVA_XMX}g  (${THREADS} threads each)"
log "STAR alignment  : ${STAR_PAR} parallel"
log "Resume RUN_ID   : ${RUN_ID}"
log "============================================================"

# =============================================================================
# 1. Stop instance
# =============================================================================
log "Stopping instance ${INSTANCE_ID}..."
aws ec2 stop-instances --instance-ids "${INSTANCE_ID}" --region "${REGION}" --output text > /dev/null

log "Waiting for instance to stop..."
aws ec2 wait instance-stopped --instance-ids "${INSTANCE_ID}" --region "${REGION}"
log "Instance stopped."

# =============================================================================
# 2. Resize
# =============================================================================
log "Resizing to ${TARGET_TYPE}..."
aws ec2 modify-instance-attribute \
    --instance-id "${INSTANCE_ID}" \
    --instance-type "{\"Value\": \"${TARGET_TYPE}\"}" \
    --region "${REGION}"
log "Instance type updated."

# =============================================================================
# 3. Start and wait for SSH
# =============================================================================
log "Starting instance..."
aws ec2 start-instances --instance-ids "${INSTANCE_ID}" --region "${REGION}" --output text > /dev/null

log "Waiting for instance to be running..."
aws ec2 wait instance-running --instance-ids "${INSTANCE_ID}" --region "${REGION}"

PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids "${INSTANCE_ID}" \
    --region "${REGION}" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' \
    --output text)
log "Instance running at ${PUBLIC_IP}"

log "Waiting for SSH to be ready..."
for i in $(seq 1 30); do
    if ssh -i "${SSH_KEY}" -o StrictHostKeyChecking=no -o ConnectTimeout=5 \
            "${SSH_USER}@${PUBLIC_IP}" "echo ok" &>/dev/null; then
        log "SSH is ready."
        break
    fi
    [[ $i -eq 30 ]] && { log "ERROR: SSH not ready after 5 min. Check security groups."; exit 1; }
    sleep 10
done

# =============================================================================
# 4. Remote setup: kill leftovers, update configs, launch pipelines
# =============================================================================
log "Connecting to ${PUBLIC_IP} to configure and launch pipelines..."

ssh -i "${SSH_KEY}" -o StrictHostKeyChecking=no "${SSH_USER}@${PUBLIC_IP}" \
    THREADS="${THREADS}" \
    JAVA_XMX="${JAVA_XMX}" \
    MUTECT_PAR="${MUTECT_PAR}" \
    STAR_PAR="${STAR_PAR}" \
    RUN_ID="${RUN_ID}" \
    'bash -s' << 'REMOTE'

set -euo pipefail
log() { echo "[$(date '+%H:%M:%S')] $*"; }

# --------------------------------------------------------------------------
# Kill any leftover pipeline processes from the previous run
# --------------------------------------------------------------------------
log "Killing any leftover pipeline processes..."
pkill -f "02_mutect2.sh"    2>/dev/null || true
pkill -f "run_rnaseq"       2>/dev/null || true
pkill -f "gatk.*Mutect2"    2>/dev/null || true
pkill -f "STAR.*genomeDir"  2>/dev/null || true
pkill -f "featureCounts"    2>/dev/null || true
sleep 2
log "Done."

# --------------------------------------------------------------------------
# Update config.sh
# --------------------------------------------------------------------------
CONFIG="/home/ec2-user/NeoAntigen2026-aws-rerun/code/config.sh"
log "Updating ${CONFIG}: THREADS=${THREADS}, JAVA_OPTS=-Xmx${JAVA_XMX}g"
sed -i "s/^export THREADS=.*/export THREADS=${THREADS}/" "${CONFIG}"
sed -i "s/^export JAVA_OPTS=.*/export JAVA_OPTS=\"-Xmx${JAVA_XMX}g\"/" "${CONFIG}"

# --------------------------------------------------------------------------
# Update config_rnaseq.sh
# --------------------------------------------------------------------------
RNA_CONFIG="/home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh"
log "Updating ${RNA_CONFIG}: RNASEQ_THREADS=${THREADS}, STAR_ALIGN_PARALLEL=${STAR_PAR}"
sed -i "s/^export RNASEQ_THREADS=.*/export RNASEQ_THREADS=${THREADS}/" "${RNA_CONFIG}"
sed -i "s/^export STAR_ALIGN_PARALLEL=.*/export STAR_ALIGN_PARALLEL=${STAR_PAR}/" "${RNA_CONFIG}"

# --------------------------------------------------------------------------
# Clear incomplete STAR index so R01 rebuilds cleanly
# --------------------------------------------------------------------------
STAR_INDEX="/home/ec2-user/ref/mm10/star_index"
if [[ ! -f "${STAR_INDEX}/SAindex" ]]; then
    log "Clearing incomplete STAR index..."
    rm -f "${STAR_INDEX}"/*
fi

# --------------------------------------------------------------------------
# Launch WES: step 2 onwards (PoN + preprocessed BAMs already on S3)
# --------------------------------------------------------------------------
log "Launching WES pipeline (step 2) in screen session 'wes'..."
screen -dmS wes bash -c "
    export PATH=/home/ec2-user/miniforge3/bin:\$PATH
    export RUN_ID=${RUN_ID}
    export MUTECT2_PARALLEL=${MUTECT_PAR}
    cd /home/ec2-user
    bash /home/ec2-user/NeoAntigen2026-aws-rerun/code/run_pipeline.sh 2 12 \
        2>&1 | tee -a /home/ec2-user/results/${RUN_ID}/wes_pipeline.log
    echo '[WES] Pipeline exited with code \$?' >> /home/ec2-user/results/${RUN_ID}/wes_pipeline.log
"

# --------------------------------------------------------------------------
# Launch RNAseq: full pipeline R1→R4 (S3 will cache the STAR index after first build)
# --------------------------------------------------------------------------
log "Launching RNAseq pipeline (R1) in screen session 'rnaseq'..."
screen -dmS rnaseq bash -c "
    export PATH=/home/ec2-user/miniforge3/bin:\$PATH
    export RUN_ID=${RUN_ID}
    cd /home/ec2-user
    bash /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/run_rnaseq_pipeline.sh R1 \
        2>&1 | tee -a /home/ec2-user/results/${RUN_ID}/rnaseq_pipeline.log
    echo '[RNASEQ] Pipeline exited with code \$?' >> /home/ec2-user/results/${RUN_ID}/rnaseq_pipeline.log
"

log "============================================================"
log "Both pipelines launched. Monitor with:"
log "  screen -r wes        # WES pipeline"
log "  screen -r rnaseq     # RNAseq pipeline"
log ""
log "Or tail logs directly:"
log "  tail -f ~/results/${RUN_ID}/wes_pipeline.log"
log "  tail -f ~/results/${RUN_ID}/rnaseq_pipeline.log"
log "============================================================"

REMOTE

log "============================================================"
log "Done! Instance ${INSTANCE_ID} is now a ${TARGET_TYPE} at ${PUBLIC_IP}"
log "SSH:  ssh -i ${SSH_KEY} ${SSH_USER}@${PUBLIC_IP}"
log "============================================================"
