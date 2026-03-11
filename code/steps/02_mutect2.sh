#!/usr/bin/env bash
# =============================================================================
# Step 2 — Somatic Variant Calling (GATK4 Mutect2)
#
# Workflow per tumour sample:
#   2a. Panel of Normals  — Mutect2 tumour-only on both D0 normals (parallel)
#                           → GenomicsDBImport → CreateSomaticPanelOfNormals
#   2b. Common SNP sites prep for contamination estimation
#   2c. Download matched normal (423_D0_old, shared across all tumours)
#   2d. Per-tumour (parallel, MUTECT2_PARALLEL at a time):
#         - Mutect2 tumour-vs-normal  (--f1r2-tar-gz for orientation model)
#         - LearnReadOrientationModel (OxoG / FFPE orientation bias)
#         - GetPileupSummaries        (tumour + normal, for contamination)
#         - CalculateContamination
#         - Upload all artefacts → S3 (for FilterMutectCalls in Step 3)
#
# Preprocessed BAMs are pulled from the permanent store:
#   s3://neoantigen2026-rerun/data/bam/wes/preprocessed/
# The script polls that location and processes each tumour as soon as it
# appears — run it in parallel with ongoing Step 1 preprocessing.
#
# Env overrides:
#   MUTECT2_PARALLEL=2   (tumours called concurrently; default 2)
#   POLL_INTERVAL=10     (minutes between S3 polls; default 10)
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

OUTDIR="${RESULTS_DIR}/mutect2"
S3_OUT="${S3_RESULTS}/mutect2"
TMP_M="${TMP_DIR}/mutect2"
mkdir -p "${OUTDIR}/pon" "${OUTDIR}/vcf_raw" "${OUTDIR}/stats" \
         "${OUTDIR}/pileups" "${OUTDIR}/contamination" \
         "${OUTDIR}/orientation" "${TMP_M}"

MUTECT2_PARALLEL="${MUTECT2_PARALLEL:-2}"
POLL_INTERVAL="${POLL_INTERVAL:-10}"   # minutes

# Restrict calling to exon targets (WES — much faster, fewer false positives)
if [[ -f "${EXON_BED}" ]]; then
    L_ARGS="-L ${EXON_BED} --interval-padding 100"
    log "Intervals: ${EXON_BED} (+100 bp padding)"
else
    log "[WARN] No exon BED — calling whole genome"
    L_ARGS=""
fi

# =============================================================================
# Helper: download a preprocessed BAM from permanent store or run results.
# Waits up to MAX_WAIT_MIN minutes polling S3 if not yet available.
# Returns 1 (and logs skip) if still unavailable after the wait.
# =============================================================================
download_preproc() {
    local SAMPLE="$1" DEST="$2" MAX_WAIT_MIN="${3:-0}"
    local WAITED=0

    [[ -f "${DEST}" ]] && return 0

    while true; do
        # Permanent store first (canonical location across runs)
        if aws s3 ls "${S3_BAM_PREPROC}/${SAMPLE}.preproc.bam" &>/dev/null; then
            log "  Downloading ${SAMPLE} from permanent store..."
            aws s3 cp "${S3_BAM_PREPROC}/${SAMPLE}.preproc.bam" "${DEST}"
            aws s3 cp "${S3_BAM_PREPROC}/${SAMPLE}.preproc.bai" \
                      "${DEST%.bam}.bai" 2>/dev/null || \
            aws s3 cp "${S3_BAM_PREPROC}/${SAMPLE}.preproc.bam.bai" \
                      "${DEST}.bai"      2>/dev/null || true
            return 0
        fi
        # Fallback: current run's preprocessed_bams prefix
        if aws s3 ls "${S3_RESULTS}/preprocessed_bams/${SAMPLE}.preproc.bam" &>/dev/null; then
            log "  Downloading ${SAMPLE} from run results..."
            aws s3 cp "${S3_RESULTS}/preprocessed_bams/${SAMPLE}.preproc.bam" "${DEST}"
            aws s3 cp "${S3_RESULTS}/preprocessed_bams/${SAMPLE}.preproc.bai" \
                      "${DEST%.bam}.bai" 2>/dev/null || \
            aws s3 cp "${S3_RESULTS}/preprocessed_bams/${SAMPLE}.preproc.bam.bai" \
                      "${DEST}.bai"      2>/dev/null || true
            return 0
        fi

        if [[ "${WAITED}" -ge "${MAX_WAIT_MIN}" ]]; then
            log "  [SKIP] ${SAMPLE} not preprocessed after ${WAITED} min"
            return 1
        fi

        log "  [WAIT] ${SAMPLE} not yet on S3 — polling again in ${POLL_INTERVAL} min..."
        sleep $(( POLL_INTERVAL * 60 ))
        WAITED=$(( WAITED + POLL_INTERVAL ))
    done
}

# =============================================================================
# 2a — Panel of Normals (both D0 normals in parallel)
# =============================================================================
PON_VCF="${OUTDIR}/pon/pon.vcf.gz"

if aws s3 ls "${S3_OUT}/pon/pon.vcf.gz" &>/dev/null; then
    log "PoN already on S3 — downloading..."
    aws s3 cp "${S3_OUT}/pon/pon.vcf.gz"     "${PON_VCF}"
    aws s3 cp "${S3_OUT}/pon/pon.vcf.gz.tbi" "${PON_VCF}.tbi"
elif [[ ! -f "${PON_VCF}" ]]; then
    log "========================================"
    log "Step 2a — Building Panel of Normals"
    log "========================================"

    PON_PIDS=()
    for NORM in "${PON_NORMALS[@]}"; do
        NORM_VCF="${OUTDIR}/pon/${NORM}.vcf.gz"
        if [[ -f "${NORM_VCF}" ]]; then
            log "  PoN VCF already exists: ${NORM}"
            continue
        fi
        (
            LOCAL_NORM="${TMP_M}/pon_${NORM}.preproc.bam"
            download_preproc "${NORM}" "${LOCAL_NORM}" 120 || exit 0
            log "  Mutect2 tumour-only (PoN): ${NORM}"
            ${GATK} --java-options "${JAVA_OPTS}" Mutect2 \
                -R "${REF}" \
                -I "${LOCAL_NORM}" \
                --max-mnp-distance 0 \
                --native-pair-hmm-threads "${THREADS}" \
                ${L_ARGS} \
                -O "${NORM_VCF}" 2>&1 \
                | grep -E "ProgressMeter|ERROR|Exception" | tail -5 || true
            rm -f "${LOCAL_NORM}" "${LOCAL_NORM%.bam}.bai" "${LOCAL_NORM}.bai"
            log "  PoN Mutect2 done: ${NORM}"
        ) &
        PON_PIDS+=($!)
    done
    for PID in "${PON_PIDS[@]}"; do wait "${PID}"; done

    # Collect VCFs
    PON_V_ARGS=""
    for NORM in "${PON_NORMALS[@]}"; do
        NORM_VCF="${OUTDIR}/pon/${NORM}.vcf.gz"
        [[ -f "${NORM_VCF}" ]] && PON_V_ARGS="${PON_V_ARGS} -V ${NORM_VCF}"
    done
    [[ -z "${PON_V_ARGS}" ]] && { log "[ERROR] No PoN VCFs — cannot build PoN"; exit 1; }

    log "  GenomicsDBImport..."
    PON_DB="${OUTDIR}/pon/pon_db"
    rm -rf "${PON_DB}"
    ${GATK} --java-options "${JAVA_OPTS}" GenomicsDBImport \
        -R "${REF}" \
        ${PON_V_ARGS} \
        --genomicsdb-workspace-path "${PON_DB}" \
        ${L_ARGS} 2>&1 | grep -E "ProgressMeter|INFO.*Done|ERROR" | tail -5

    log "  CreateSomaticPanelOfNormals..."
    ${GATK} --java-options "${JAVA_OPTS}" CreateSomaticPanelOfNormals \
        -R "${REF}" \
        -V "gendb://${PON_DB}" \
        -O "${PON_VCF}" 2>&1 | grep -E "ProgressMeter|INFO.*Done|ERROR" | tail -3

    aws s3 cp "${PON_VCF}"     "${S3_OUT}/pon/pon.vcf.gz"
    aws s3 cp "${PON_VCF}.tbi" "${S3_OUT}/pon/pon.vcf.gz.tbi"
    log "  PoN complete → ${PON_VCF}"
fi

# =============================================================================
# 2b — Common SNP sites for GetPileupSummaries (contamination estimation)
#      Annotates MGP SNPs with AF=0.5 (all are known polymorphic sites).
# =============================================================================
COMMON_SITES="${OUTDIR}/common_snps.vcf.gz"
if [[ ! -f "${COMMON_SITES}" ]] && [[ -f "${KNOWN_SNPS}" ]]; then
    log "Preparing common SNP sites for GetPileupSummaries..."
    ${BCFTOOLS} view "${KNOWN_SNPS}" \
        | awk 'BEGIN{OFS="\t"}
               /^#/ { print; next }
               { $8 = ($8 == "." ? "AF=0.5" : $8";AF=0.5"); print }' \
        | ${CONDA_BIN}/bgzip -c > "${COMMON_SITES}"
    ${CONDA_BIN}/tabix -p vcf "${COMMON_SITES}"
    log "  Common sites ready: ${COMMON_SITES}"
fi

# =============================================================================
# 2c — Download matched normal (kept in TMP for all tumour calls)
# =============================================================================
LOCAL_NORMAL="${TMP_M}/${NORMAL_SAMPLE}.preproc.bam"
log "Downloading matched normal: ${NORMAL_SAMPLE} (SM: ${NORMAL_SM})"
download_preproc "${NORMAL_SAMPLE}" "${LOCAL_NORMAL}" 240

# =============================================================================
# 2d — Per-tumour calling function
# =============================================================================
call_tumor() {
    local TUMOR="$1"
    local TUMOR_SM="${TUMOR_SM_TAGS[${TUMOR}]}"

    local RAW_VCF="${OUTDIR}/vcf_raw/${TUMOR}.vcf.gz"
    local STATS_F="${OUTDIR}/stats/${TUMOR}.stats"
    local F1R2_F="${TMP_M}/${TUMOR}.f1r2.tar.gz"
    local ORIENT_F="${OUTDIR}/orientation/${TUMOR}.read_orientation_model.tar.gz"
    local PILEUP_T="${OUTDIR}/pileups/${TUMOR}.pileups.table"
    local PILEUP_N="${OUTDIR}/pileups/${NORMAL_SAMPLE}.pileups.table"
    local CONTAM_F="${OUTDIR}/contamination/${TUMOR}.contamination.table"
    local SEGMENT_F="${OUTDIR}/contamination/${TUMOR}.segments.table"

    # Skip if all outputs already on S3
    if aws s3 ls "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz"        &>/dev/null && \
       aws s3 ls "${S3_OUT}/contamination/${TUMOR}.contamination.table" &>/dev/null; then
        log "[SKIP] ${TUMOR} — all outputs already on S3"
        aws s3 cp "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz"         "${RAW_VCF}"   2>/dev/null || true
        aws s3 cp "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz.tbi"     "${RAW_VCF}.tbi" 2>/dev/null || true
        aws s3 cp "${S3_OUT}/stats/${TUMOR}.stats"             "${STATS_F}"   2>/dev/null || true
        aws s3 cp "${S3_OUT}/orientation/${TUMOR}.read_orientation_model.tar.gz" \
                  "${ORIENT_F}" 2>/dev/null || true
        aws s3 cp "${S3_OUT}/contamination/${TUMOR}.contamination.table" "${CONTAM_F}"  2>/dev/null || true
        aws s3 cp "${S3_OUT}/contamination/${TUMOR}.segments.table"      "${SEGMENT_F}" 2>/dev/null || true
        return 0
    fi

    log "========================================"
    log "Step 2d — Mutect2: ${TUMOR} (SM: ${TUMOR_SM}) vs ${NORMAL_SAMPLE} (SM: ${NORMAL_SM})"
    log "========================================"

    # Wait up to 4 h for preprocessing to finish
    local LOCAL_TUMOR="${TMP_M}/${TUMOR}.preproc.bam"
    download_preproc "${TUMOR}" "${LOCAL_TUMOR}" 240 || {
        log "[SKIP] ${TUMOR} — preprocessed BAM never appeared"
        return 0
    }

    # --- Mutect2 --------------------------------------------------------------
    log "  [${TUMOR}] Mutect2..."
    ${GATK} --java-options "${JAVA_OPTS}" Mutect2 \
        -R "${REF}" \
        -I "${LOCAL_TUMOR}" \
        -I "${LOCAL_NORMAL}" \
        -normal "${NORMAL_SM}" \
        --panel-of-normals "${PON_VCF}" \
        --f1r2-tar-gz "${F1R2_F}" \
        --native-pair-hmm-threads "${THREADS}" \
        ${L_ARGS} \
        -O "${RAW_VCF}" \
        --stats "${STATS_F}" 2>&1 \
        | grep -E "ProgressMeter|ERROR|Exception" | tail -5

    [[ ! -s "${RAW_VCF}" ]] && {
        log "[ERROR] Mutect2 produced no VCF for ${TUMOR}"
        rm -f "${LOCAL_TUMOR}" "${LOCAL_TUMOR%.bam}.bai" "${LOCAL_TUMOR}.bai"
        return 1
    }

    # --- LearnReadOrientationModel -------------------------------------------
    log "  [${TUMOR}] LearnReadOrientationModel..."
    ${GATK} --java-options "${JAVA_OPTS}" LearnReadOrientationModel \
        -I "${F1R2_F}" \
        -O "${ORIENT_F}" 2>&1 | tail -2 || true
    rm -f "${F1R2_F}"

    # --- GetPileupSummaries + CalculateContamination -------------------------
    if [[ -f "${COMMON_SITES}" ]]; then
        log "  [${TUMOR}] GetPileupSummaries (tumour)..."
        ${GATK} --java-options "-Xmx4g" GetPileupSummaries \
            -R "${REF}" \
            -I "${LOCAL_TUMOR}" \
            -V "${COMMON_SITES}" \
            -L "${COMMON_SITES}" \
            -O "${PILEUP_T}" 2>&1 | tail -2 || true

        if [[ ! -f "${PILEUP_N}" ]]; then
            log "  [${TUMOR}] GetPileupSummaries (normal)..."
            ${GATK} --java-options "-Xmx4g" GetPileupSummaries \
                -R "${REF}" \
                -I "${LOCAL_NORMAL}" \
                -V "${COMMON_SITES}" \
                -L "${COMMON_SITES}" \
                -O "${PILEUP_N}" 2>&1 | tail -2 || true
        fi

        if [[ -f "${PILEUP_T}" && -s "${PILEUP_T}" ]]; then
            log "  [${TUMOR}] CalculateContamination..."
            local CONTAM_ARGS="-I ${PILEUP_T} -O ${CONTAM_F} --tumor-segmentation ${SEGMENT_F}"
            [[ -f "${PILEUP_N}" && -s "${PILEUP_N}" ]] && \
                CONTAM_ARGS="${CONTAM_ARGS} --matched-normal ${PILEUP_N}"
            ${GATK} --java-options "-Xmx4g" CalculateContamination \
                ${CONTAM_ARGS} 2>&1 | tail -2 || true
        else
            log "  [${TUMOR}] [WARN] Pileup empty — skipping CalculateContamination"
        fi
    else
        log "  [${TUMOR}] [WARN] No common sites VCF — skipping contamination estimation"
    fi

    # --- Upload to S3 --------------------------------------------------------
    log "  [${TUMOR}] Uploading artefacts to S3..."
    aws s3 cp "${RAW_VCF}"     "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz"
    aws s3 cp "${RAW_VCF}.tbi" "${S3_OUT}/vcf_raw/${TUMOR}.vcf.gz.tbi"  2>/dev/null || true
    aws s3 cp "${STATS_F}"     "${S3_OUT}/stats/${TUMOR}.stats"          2>/dev/null || true
    [[ -f "${ORIENT_F}" ]]  && aws s3 cp "${ORIENT_F}"  "${S3_OUT}/orientation/${TUMOR}.read_orientation_model.tar.gz"
    [[ -f "${CONTAM_F}" ]]  && aws s3 cp "${CONTAM_F}"  "${S3_OUT}/contamination/${TUMOR}.contamination.table"
    [[ -f "${SEGMENT_F}" ]] && aws s3 cp "${SEGMENT_F}" "${S3_OUT}/contamination/${TUMOR}.segments.table"
    [[ -f "${PILEUP_T}" ]]  && aws s3 cp "${PILEUP_T}"  "${S3_OUT}/pileups/${TUMOR}.pileups.table" 2>/dev/null || true

    rm -f "${LOCAL_TUMOR}" "${LOCAL_TUMOR%.bam}.bai" "${LOCAL_TUMOR}.bai"
    log "  Done: ${TUMOR}"
}

# =============================================================================
# 2e — Dispatch tumours (MUTECT2_PARALLEL concurrent)
# =============================================================================
log "========================================"
log "Calling ${#TUMOR_SAMPLES[@]} tumours (${MUTECT2_PARALLEL} parallel)"
log "========================================"

PIDS=()
FAILED=()

for TUMOR in "${TUMOR_SAMPLES[@]}"; do
    # Throttle: wait until a slot opens
    while [[ ${#PIDS[@]} -ge ${MUTECT2_PARALLEL} ]]; do
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                wait "${PIDS[$i]}" || FAILED+=("slot_${i}")
                unset 'PIDS[$i]'
                PIDS=("${PIDS[@]}")
            fi
        done
        sleep 10
    done

    call_tumor "${TUMOR}" &
    PIDS+=($!)
    log "  Dispatched: ${TUMOR} (PID ${PIDS[-1]})"
done

for PID in "${PIDS[@]}"; do
    wait "${PID}" || FAILED+=("${PID}")
done

# Cleanup shared normal BAM
rm -f "${LOCAL_NORMAL}" "${LOCAL_NORMAL%.bam}.bai" "${LOCAL_NORMAL}.bai"

if [[ ${#FAILED[@]} -gt 0 ]]; then
    log "[ERROR] Failed jobs: ${FAILED[*]}"
    exit 1
fi

log "Step 2 complete."
log "  Raw VCFs     : ${S3_OUT}/vcf_raw/"
log "  PoN          : ${S3_OUT}/pon/"
log "  Orientation  : ${S3_OUT}/orientation/"
log "  Contamination: ${S3_OUT}/contamination/"
