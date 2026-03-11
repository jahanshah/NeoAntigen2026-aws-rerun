#!/bin/bash
# =============================================================================
# index_mm10_and_run.sh
# Wait for mm10.fa and all BAMs to finish downloading, then:
#   1. Index mm10.fa  (samtools faidx, GATK dict, BWA index)
#   2. Check SnpEff GRCm38.86 database
#   3. Launch neoantigen_pipeline.sh
#
# Run in background after starting downloads:
#   nohup bash scripts/index_mm10_and_run.sh > logs/pipeline_launch.log 2>&1 &
# =============================================================================

set -euo pipefail
WORK_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="${WORK_DIR}/data"
REF_DB="${WORK_DIR}/reference_db"
LOG_DIR="${WORK_DIR}/logs"
mkdir -p "$LOG_DIR"

REF_FA="${DATA_DIR}/mm10.fa"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── 1. Wait for mm10.fa ───────────────────────────────────────────────────────
log "Waiting for mm10.fa download to complete..."
while [[ ! -f "$REF_FA" ]] || ls "${DATA_DIR}/Unconfirmed"*.crdownload 2>/dev/null | grep -q .; do
    sleep 30
    log "  Still waiting — $(ls -lh "$DATA_DIR" 2>/dev/null | grep -v '^total' | wc -l) files in data/"
done
log "All downloads complete."

# ── 2. Index mm10.fa ──────────────────────────────────────────────────────────
log "=== Indexing mm10.fa ==="

# samtools faidx
if [[ ! -f "${REF_FA}.fai" ]]; then
    log "  samtools faidx mm10.fa..."
    samtools faidx "$REF_FA" 2>"${LOG_DIR}/mm10_faidx.log"
    log "  faidx done."
else
    log "  mm10.fa.fai already exists — skipping."
fi

# GATK sequence dictionary
DICT="${DATA_DIR}/mm10.dict"
if [[ ! -f "$DICT" ]]; then
    log "  Creating GATK sequence dictionary..."
    gatk CreateSequenceDictionary -R "$REF_FA" -O "$DICT" \
        2>"${LOG_DIR}/mm10_dict.log"
    log "  Dict done."
else
    log "  mm10.dict already exists — skipping."
fi

# BWA index (takes ~90 min for 2.6 GB genome)
if [[ ! -f "${REF_FA}.bwt" ]]; then
    log "  BWA index mm10.fa (this takes ~60-90 min)..."
    bwa index "$REF_FA" 2>"${LOG_DIR}/mm10_bwa_index.log"
    log "  BWA index done."
else
    log "  BWA index already exists — skipping."
fi

# ── 3. Wait for reference_db files ───────────────────────────────────────────
log "=== Checking reference_db files ==="
for f in "${REF_DB}/merged_exons.bed" "${REF_DB}/proteins.fa"; do
    while [[ ! -f "$f" ]]; do
        log "  Waiting for $(basename $f)..."
        sleep 60
    done
    log "  ✔  $(basename $f)"
done

# ── 4. Check SnpEff GRCm38.86 database ───────────────────────────────────────
log "=== Checking SnpEff GRCm38.86 ==="
if snpEff databases 2>/dev/null | grep -q "GRCm38.86"; then
    log "  ✔  SnpEff GRCm38.86 database present."
else
    log "  Downloading SnpEff GRCm38.86 database..."
    snpEff download GRCm38.86 2>"${LOG_DIR}/snpeff_download.log"
    log "  SnpEff database downloaded."
fi

# ── 5. Launch pipeline ────────────────────────────────────────────────────────
log "=== All prerequisites met — launching neoantigen_pipeline.sh ==="
bash "${WORK_DIR}/neoantigen_pipeline.sh" 2>&1 | tee "${LOG_DIR}/pipeline_main.log"
