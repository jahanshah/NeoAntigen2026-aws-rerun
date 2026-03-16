#!/bin/bash
# Step 11 — Sequenza Tumor Purity & Copy Number Estimation
#
# Estimates tumor purity and ploidy per WES sample using matched-normal Sequenza.
# Inputs:  pre-processed WES BAMs (tumor + normal 423_D0_old)
# Outputs: purity/ploidy per sample → used to rerun PyClone-VI

set -euo pipefail
export PATH="/home/ec2-user/miniforge3/bin:$PATH"

BASE_DIR="/home/ec2-user"
REF="${BASE_DIR}/ref/mm10/mm10.fa"
GC_WIG="${BASE_DIR}/ref/mm10/mm10.gc50.wig.gz"
OUT_DIR="${BASE_DIR}/results/res_20260311_225555/sequenza"
TMP_DIR="${BASE_DIR}/tmp/neoantig_pipeline/sequenza_bam"
S3_WES="s3://neoantigen2026-rerun/pre-processed/wes"
S3_RESULTS="s3://neoantigen2026-rerun/results/res_20260311_225555/sequenza"
THREADS=8
LOG="${BASE_DIR}/results/res_20260311_225555/logs/sequenza_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "${OUT_DIR}" "${TMP_DIR}"
exec > >(tee -a "${LOG}") 2>&1

echo "[$(date)] === Step 11: Sequenza Purity Estimation ==="
echo "[$(date)] Reference: ${REF}"
echo "[$(date)] GC wiggle: ${GC_WIG}"

# Wait for GC wiggle if still generating
while [ ! -f "${GC_WIG}" ]; do
    echo "[$(date)] Waiting for GC wiggle file..."
    sleep 60
done

# ── Normal BAM ─────────────────────────────────────────────────────────────
NORMAL_BAM="${TMP_DIR}/423_D0_old.preproc.bam"
NORMAL_BAI="${TMP_DIR}/423_D0_old.preproc.bai"

if [ ! -f "${NORMAL_BAM}" ]; then
    echo "[$(date)] Downloading normal BAM (423_D0_old)..."
    aws s3 cp "${S3_WES}/423_D0_old.preproc.bam" "${NORMAL_BAM}"
    aws s3 cp "${S3_WES}/423_D0_old.preproc.bai" "${NORMAL_BAI}"
    # samtools needs filename.bam.bai or filename.bai — .bai works
fi
echo "[$(date)] Normal BAM: ${NORMAL_BAM}"

# ── Tumor samples ───────────────────────────────────────────────────────────
declare -A TUMOR_S3=(
    ["36_D99_new"]="${S3_WES}/36_D99_new.preproc.bam"
    ["38_D99_new"]="${S3_WES}/38_D99_new.preproc.bam"
    ["42_D122_old"]="${S3_WES}/42_D122_old.preproc.bam"
    ["428_D20_new"]="${S3_WES}/428_D20_new.preproc.bam"
    ["443_D21_new"]="${S3_WES}/443_D21_new.preproc.bam"
)

# Local BAMs for samples already downloaded
declare -A TUMOR_LOCAL=(
    ["428_D20_new"]="${BASE_DIR}/tmp/neoantig_pipeline/mutect2/428_D20_new.preproc.bam"
    ["443_D21_new"]="${BASE_DIR}/tmp/neoantig_pipeline/mutect2/443_D21_new.preproc.bam"
)

run_sequenza() {
    local SAMPLE="$1"
    local TUMOR_BAM="$2"

    local SEQZ_RAW="${OUT_DIR}/${SAMPLE}.seqz.gz"
    local SEQZ_BIN="${OUT_DIR}/${SAMPLE}.binned.seqz.gz"
    local SAMPLE_OUT="${OUT_DIR}/${SAMPLE}"

    echo "[$(date)] --- ${SAMPLE} ---"

    # ── seqz generation ─────────────────────────────────────────────────
    if [ ! -f "${SEQZ_RAW}" ]; then
        echo "[$(date)]   bam2seqz: ${TUMOR_BAM}"
        sequenza-utils bam2seqz \
            --normal  "${NORMAL_BAM}" \
            --tumor   "${TUMOR_BAM}" \
            --fasta   "${REF}" \
            --gc_file "${GC_WIG}" \
            --output  "${SEQZ_RAW}" \
            --parallel "${THREADS}"
        echo "[$(date)]   seqz done: $(ls -lh ${SEQZ_RAW} | awk '{print $5}')"
    else
        echo "[$(date)]   seqz exists, skipping"
    fi

    # ── binning ─────────────────────────────────────────────────────────
    if [ ! -f "${SEQZ_BIN}" ]; then
        sequenza-utils seqz_binning \
            --seqz "${SEQZ_RAW}" \
            --window 300 \
            --output "${SEQZ_BIN}"
        echo "[$(date)]   binned: $(ls -lh ${SEQZ_BIN} | awk '{print $5}')"
    fi

    # ── R analysis ──────────────────────────────────────────────────────
    mkdir -p "${SAMPLE_OUT}"
    Rscript - "${SEQZ_BIN}" "${SAMPLE_OUT}" "${SAMPLE}" << 'REOF'
args    <- commandArgs(trailingOnly=TRUE)
seqz    <- args[1]
outdir  <- args[2]
sname   <- args[3]

library(sequenza)

# Load and extract CNV + BAF data
seqzdata <- sequenza.extract(seqz, verbose=FALSE)

# Fit purity/ploidy grid
CP <- sequenza.fit(seqzdata)

# Save all outputs (segments, copy number, purity)
sequenza.results(
    sequenza.extract = seqzdata,
    cp.table         = CP,
    sample.id        = sname,
    out.dir          = outdir
)

# Extract best purity/ploidy estimate
best <- CP$values.matrix[which.max(CP$log.posterior), ]
purity <- best["cellularity"]
ploidy <- best["ploidy"]
cat(sprintf("PURITY_RESULT\t%s\t%.4f\t%.4f\n", sname, purity, ploidy))
REOF

    # Parse purity from R output (captured in log)
    echo "[$(date)]   R analysis complete for ${SAMPLE}"

    # Upload results to S3
    aws s3 cp "${SAMPLE_OUT}/" "${S3_RESULTS}/${SAMPLE}/" --recursive --quiet
    echo "[$(date)]   Uploaded to S3: ${S3_RESULTS}/${SAMPLE}/"
}

# ── Run each sample ─────────────────────────────────────────────────────────
for SAMPLE in "${!TUMOR_S3[@]}"; do
    # Use local BAM if available, otherwise download
    if [ -n "${TUMOR_LOCAL[$SAMPLE]+set}" ] && [ -f "${TUMOR_LOCAL[$SAMPLE]}" ]; then
        TUMOR_BAM="${TUMOR_LOCAL[$SAMPLE]}"
        echo "[$(date)] Using local BAM for ${SAMPLE}"
    else
        TUMOR_BAM="${TMP_DIR}/${SAMPLE}.preproc.bam"
        TUMOR_BAI="${TMP_DIR}/${SAMPLE}.preproc.bai"
        if [ ! -f "${TUMOR_BAM}" ]; then
            echo "[$(date)] Downloading ${SAMPLE} BAM..."
            aws s3 cp "${TUMOR_S3[$SAMPLE]}" "${TUMOR_BAM}"
            aws s3 cp "${TUMOR_S3[$SAMPLE]%.bam}.bai" "${TUMOR_BAI}"
        fi
    fi

    run_sequenza "${SAMPLE}" "${TUMOR_BAM}"
done

echo "[$(date)] === All Sequenza runs complete ==="
echo "[$(date)] Purity results saved to: ${OUT_DIR}"
