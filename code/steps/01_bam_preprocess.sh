#!/usr/bin/env bash
# =============================================================================
# Step 1 — BAM Preprocessing
# Sort → MarkDuplicates → BaseQualityScoreRecalibration → ValidateSamFile
# Processes all 7 samples (normals + tumours) one at a time.
# Input BAMs streamed from S3; preprocessed BAMs uploaded back to S3.
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

OUTDIR="${RESULTS_DIR}/preprocessed_bams"
S3_PREPROC="${S3_RESULTS}/preprocessed_bams"
mkdir -p "${OUTDIR}" "${TMP_DIR}/preproc"

ALL_SAMPLES=("${NORMAL_SAMPLE}" "${PON_NORMALS[1]}" "${TUMOR_SAMPLES[@]}")

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    FINAL_BAM="${OUTDIR}/${SAMPLE}.preproc.bam"
    S3_OUT="${S3_PREPROC}/${SAMPLE}.preproc.bam"

    # Skip if already uploaded to S3
    if aws s3 ls "${S3_OUT}" &>/dev/null; then
        skip "${SAMPLE}.preproc.bam (S3)"
        continue
    fi

    log "========================================"
    log "Preprocessing: ${SAMPLE}"
    log "========================================"

    LOCAL_RAW="${TMP_DIR}/preproc/${SAMPLE}.bam"
    LOCAL_SORTED="${TMP_DIR}/preproc/${SAMPLE}.sorted.bam"
    LOCAL_DEDUP="${TMP_DIR}/preproc/${SAMPLE}.dedup.bam"
    METRICS="${OUTDIR}/${SAMPLE}.dup_metrics.txt"
    RECAL_TABLE="${OUTDIR}/${SAMPLE}.recal.table"

    # --- 1a. Download raw BAM ------------------------------------------------
    log "  Downloading from S3..."
    aws s3 cp "${S3_BAM}/${SAMPLE}.bam"     "${LOCAL_RAW}"
    aws s3 cp "${S3_BAM}/${SAMPLE}.bam.bai" "${LOCAL_RAW}.bai"

    # --- 1b. Sort (coordinate) — skip if already sorted ---------------------
    SORT_ORDER=$(${SAMTOOLS} view -H "${LOCAL_RAW}" | awk '/^@HD/{for(i=1;i<=NF;i++) if($i~/^SO:/) {sub("SO:","",$i); print $i}}')
    if [[ "${SORT_ORDER}" == "coordinate" ]]; then
        log "  Already coordinate-sorted — skipping sort"
        cp "${LOCAL_RAW}" "${LOCAL_SORTED}"
        ${SAMTOOLS} index "${LOCAL_SORTED}"
    else
        log "  Sorting..."
        ${SAMTOOLS} sort -@ ${THREADS} -o "${LOCAL_SORTED}" "${LOCAL_RAW}"
        ${SAMTOOLS} index "${LOCAL_SORTED}"
    fi
    rm -f "${LOCAL_RAW}" "${LOCAL_RAW}.bai"

    # --- 1c. MarkDuplicates (Picard) ----------------------------------------
    # Check if already marked (look for PG ID:MarkDuplicates in header)
    if ${SAMTOOLS} view -H "${LOCAL_SORTED}" | grep -q "ID:MarkDuplicates"; then
        log "  Duplicates already marked — skipping MarkDuplicates"
        cp "${LOCAL_SORTED}" "${LOCAL_DEDUP}"
        ${SAMTOOLS} index "${LOCAL_DEDUP}"
    else
        log "  Marking duplicates..."
        ${GATK} ${JAVA_OPTS:+--java-options "${JAVA_OPTS}"} MarkDuplicates \
            -I "${LOCAL_SORTED}" \
            -O "${LOCAL_DEDUP}" \
            -M "${METRICS}" \
            --REMOVE_DUPLICATES false \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY LENIENT \
            --TMP_DIR "${TMP_DIR}/preproc" 2>&1 | grep -E "INFO|WARN|ERROR" | tail -5 || true
        log "  Dup metrics: ${METRICS}"
    fi
    rm -f "${LOCAL_SORTED}" "${LOCAL_SORTED}.bai"

    # --- 1d. BQSR — BaseRecalibrator + ApplyBQSR ----------------------------
    log "  Running BaseRecalibrator..."
    BQSR_ARGS="-R ${REF} -I ${LOCAL_DEDUP} -O ${RECAL_TABLE}"
    if [[ -f "${KNOWN_SNPS}" && ! -f "${KNOWN_SNPS}.missing" ]]; then
        BQSR_ARGS="${BQSR_ARGS} --known-sites ${KNOWN_SNPS}"
    else
        log "  [WARN] No known SNPs — BQSR running without dbSNP"
        BQSR_ARGS="${BQSR_ARGS} --run-without-dbsnp-potentially-ruining-quality"
    fi
    ${GATK} ${JAVA_OPTS:+--java-options "${JAVA_OPTS}"} BaseRecalibrator \
        ${BQSR_ARGS} 2>&1 | grep -E "^23|INFO  Prog" | tail -3 || true

    log "  Applying BQSR..."
    ${GATK} ${JAVA_OPTS:+--java-options "${JAVA_OPTS}"} ApplyBQSR \
        -R "${REF}" \
        -I "${LOCAL_DEDUP}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${FINAL_BAM}" 2>&1 | grep -E "INFO  Prog" | tail -3 || true
    rm -f "${LOCAL_DEDUP}" "${LOCAL_DEDUP}.bai" "${LOCAL_DEDUP%.bam}.bai"

    # --- 1e. ValidateSamFile ------------------------------------------------
    log "  Validating BAM..."
    VALIDATE_OUT="${OUTDIR}/${SAMPLE}.validate.txt"
    ${GATK} ${JAVA_OPTS:+--java-options "${JAVA_OPTS}"} ValidateSamFile \
        -I "${FINAL_BAM}" \
        -O "${VALIDATE_OUT}" \
        --MODE SUMMARY 2>&1 | tail -3 || true
    if grep -q "No errors found" "${VALIDATE_OUT}"; then
        log "  Validation: PASS"
    else
        log "  [WARN] Validation issues — see ${VALIDATE_OUT}"
    fi

    # --- 1f. Upload to S3, clean local --------------------------------------
    log "  Uploading preprocessed BAM to S3..."
    aws s3 cp "${FINAL_BAM}"       "${S3_OUT}"
    aws s3 cp "${FINAL_BAM%.bam}.bai" "${S3_OUT%.bam}.bai" 2>/dev/null || \
    aws s3 cp "${FINAL_BAM}.bai"   "${S3_OUT}.bai" 2>/dev/null || true
    aws s3 cp "${METRICS}"         "${S3_PREPROC}/$(basename ${METRICS})" 2>/dev/null || true

    rm -f "${FINAL_BAM}" "${FINAL_BAM%.bam}.bai" "${FINAL_BAM}.bai"
    log "  Done: ${SAMPLE}"
done

log "Step 1 complete — all preprocessed BAMs on S3: ${S3_PREPROC}/"
