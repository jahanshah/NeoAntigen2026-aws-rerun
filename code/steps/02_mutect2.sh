#!/usr/bin/env bash
# =============================================================================
# Step 2 — Somatic Variant Calling with GATK4 Mutect2
# Builds Panel of Normals from D0 samples, then calls tumour vs normal.
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

S3_PREPROC="${S3_RESULTS}/preprocessed_bams"
OUTDIR="${RESULTS_DIR}/mutect2"
S3_OUT="${S3_RESULTS}/mutect2"
mkdir -p "${OUTDIR}/pon" "${OUTDIR}/vcf_raw" "${OUTDIR}/stats" "${TMP_DIR}/mutect2"

# Build -L interval args (use exon BED if available for WES speed)
L_ARGS=""
if [[ -f "${EXON_BED}" ]]; then
    L_ARGS="-L ${EXON_BED} --interval-padding 100"
    log "Interval list: ${EXON_BED} (exons +100 bp)"
else
    log "[WARN] No exon BED found — running whole-genome (slow)"
fi

# Download BAM: prefer preprocessed, fall back to raw
download_bam() {
    local SAMPLE="$1" DEST="$2"
    if [[ -f "${DEST}" ]]; then return 0; fi
    # Try preprocessed first
    if aws s3 ls "${S3_PREPROC}/${SAMPLE}.preproc.bam" &>/dev/null; then
        log "  Downloading preprocessed BAM: ${SAMPLE}"
        aws s3 cp "${S3_PREPROC}/${SAMPLE}.preproc.bam"     "${DEST}"
        aws s3 cp "${S3_PREPROC}/${SAMPLE}.preproc.bam.bai" "${DEST}.bai" 2>/dev/null || \
        aws s3 cp "${S3_PREPROC}/${SAMPLE}.preproc.bai"     "${DEST%.bam}.bai" 2>/dev/null || true
    else
        log "  Downloading raw BAM: ${SAMPLE}"
        aws s3 cp "${S3_BAM}/${SAMPLE}.bam"     "${DEST}"
        aws s3 cp "${S3_BAM}/${SAMPLE}.bam.bai" "${DEST}.bai" 2>/dev/null || true
    fi
}
# Backward-compatible alias
download_preproc() { download_bam "$@"; }

# =============================================================================
# 2a — Panel of Normals
# =============================================================================
PON_VCF="${OUTDIR}/pon/pon.vcf.gz"

if aws s3 ls "${S3_OUT}/pon/pon.vcf.gz" &>/dev/null; then
    skip "PoN (S3)"
    aws s3 cp "${S3_OUT}/pon/pon.vcf.gz"     "${PON_VCF}"
    aws s3 cp "${S3_OUT}/pon/pon.vcf.gz.tbi" "${PON_VCF}.tbi"
elif [[ ! -f "${PON_VCF}" ]]; then
    log "========================================"
    log "Step 2a — Building Panel of Normals"
    log "========================================"
    PON_V_ARGS=""
    for NORM in "${PON_NORMALS[@]}"; do
        NORM_VCF="${OUTDIR}/pon/${NORM}.vcf.gz"
        LOCAL_NORM="${TMP_DIR}/mutect2/${NORM}.preproc.bam"
        download_preproc "${NORM}" "${LOCAL_NORM}"

        if [[ ! -f "${NORM_VCF}" ]]; then
            log "  Mutect2 (tumor-only) for PoN: ${NORM}"
            ${GATK} --java-options "${JAVA_OPTS}" Mutect2 \
                -R "${REF}" -I "${LOCAL_NORM}" \
                ${L_ARGS} \
                -O "${NORM_VCF}" --max-mnp-distance 0 \
                --native-pair-hmm-threads ${THREADS} 2>&1 | grep -E "ProgressMeter|ERROR|Exception" | tail -3
        fi
        rm -f "${LOCAL_NORM}" "${LOCAL_NORM}.bai" "${LOCAL_NORM%.bam}.bai"
        PON_V_ARGS="${PON_V_ARGS} -V ${NORM_VCF}"
    done

    log "  GenomicsDBImport..."
    PON_DB="${OUTDIR}/pon/pon_db"
    rm -rf "${PON_DB}"
    DB_L_ARGS=""
    if [[ -f "${EXON_BED}" ]]; then
        DB_L_ARGS="-L ${EXON_BED} --interval-padding 100"
    else
        DB_L_ARGS="$(for c in $(seq 1 19) X Y; do echo "-L chr${c}"; done)"
    fi
    ${GATK} --java-options "${JAVA_OPTS}" GenomicsDBImport \
        -R "${REF}" ${PON_V_ARGS} \
        --genomicsdb-workspace-path "${PON_DB}" \
        ${DB_L_ARGS} 2>&1 | tail -3

    log "  CreateSomaticPanelOfNormals..."
    ${GATK} --java-options "${JAVA_OPTS}" CreateSomaticPanelOfNormals \
        -R "${REF}" -V "gendb://${PON_DB}" -O "${PON_VCF}" 2>&1 | tail -3

    aws s3 cp "${PON_VCF}"     "${S3_OUT}/pon/pon.vcf.gz"
    aws s3 cp "${PON_VCF}.tbi" "${S3_OUT}/pon/pon.vcf.gz.tbi"
    log "  PoN complete: ${PON_VCF}"
fi

# =============================================================================
# 2b — Download matched normal (keep for all tumour runs)
# =============================================================================
LOCAL_NORMAL="${TMP_DIR}/mutect2/${NORMAL_SAMPLE}.preproc.bam"
download_preproc "${NORMAL_SAMPLE}" "${LOCAL_NORMAL}"

# =============================================================================
# 2c — Mutect2 per tumour
# =============================================================================
for TUMOR in "${TUMOR_SAMPLES[@]}"; do
    RAW_VCF="${OUTDIR}/vcf_raw/${TUMOR}.vcf.gz"
    STATS="${OUTDIR}/stats/${TUMOR}.stats"

    if aws s3 ls "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz" &>/dev/null; then
        skip "${TUMOR} raw VCF (S3)"
        aws s3 cp "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz"     "${RAW_VCF}"
        aws s3 cp "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz.tbi" "${RAW_VCF}.tbi" 2>/dev/null || true
        aws s3 cp "${S3_OUT}/stats/${TUMOR}.stats"         "${STATS}"       2>/dev/null || true
        continue
    fi

    log "========================================"
    log "Step 2c — Mutect2: ${TUMOR} vs ${NORMAL_SAMPLE}"
    log "========================================"

    LOCAL_TUMOR="${TMP_DIR}/mutect2/${TUMOR}.preproc.bam"
    download_preproc "${TUMOR}" "${LOCAL_TUMOR}"

    ${GATK} --java-options "${JAVA_OPTS}" Mutect2 \
        -R "${REF}" \
        -I "${LOCAL_TUMOR}" \
        -I "${LOCAL_NORMAL}" \
        -normal "${NORMAL_SM}" \
        --panel-of-normals "${PON_VCF}" \
        ${L_ARGS} \
        -O "${RAW_VCF}" --stats "${STATS}" \
        --native-pair-hmm-threads ${THREADS} 2>&1 | grep -E "ProgressMeter|ERROR|Exception" | tail -3

    aws s3 cp "${RAW_VCF}"     "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz"
    aws s3 cp "${RAW_VCF}.tbi" "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz.tbi" 2>/dev/null || true
    aws s3 cp "${STATS}"       "${S3_OUT}/stats/${TUMOR}.stats"

    rm -f "${LOCAL_TUMOR}" "${LOCAL_TUMOR%.bam}.bai" "${LOCAL_TUMOR}.bai"
    log "  Done: ${TUMOR}"
done

rm -f "${LOCAL_NORMAL}" "${LOCAL_NORMAL%.bam}.bai" "${LOCAL_NORMAL}.bai"
log "Step 2 complete."
