#!/usr/bin/env bash
# =============================================================================
# bam_preprocessing.sh — BAM sort, MarkDuplicates, BQSR, ValidateSamFile
# Wrapper for Step 1 of the NeoAntigen pipeline.
#
# Usage:
#   bash /home/ec2-user/code/bam_preprocessing.sh [SAMPLE]
#   bash /home/ec2-user/code/bam_preprocessing.sh          # all 7 samples
#   bash /home/ec2-user/code/bam_preprocessing.sh 34_D52_old
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

LOG="${RESULTS_DIR}/bam_preprocessing.log"
exec > >(tee -a "${LOG}") 2>&1

log "bam_preprocessing.sh started: $(date)"
log "Samples: ${NORMAL_SAMPLE} ${PON_NORMALS[1]} ${TUMOR_SAMPLES[*]}"
log "Log: ${LOG}"

bash /home/ec2-user/code/steps/01_bam_preprocess.sh "$@"

log "bam_preprocessing.sh complete: $(date)"
