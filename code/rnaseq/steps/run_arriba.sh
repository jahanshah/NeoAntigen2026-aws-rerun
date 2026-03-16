#!/usr/bin/env bash
# =============================================================================
# run_arriba.sh — Re-align with STAR (chimeric mode) + Arriba fusion calling
#
# Re-runs STAR with chimeric read flags required by Arriba (the original R02
# alignment omitted these). Processes samples in parallel (2 at a time).
#
# Usage:
#   RUN_ID=res_20260313_071626 bash run_arriba.sh [SAMPLE1 SAMPLE2 ...]
#   If no samples given, runs all samples in the manifest.
#
# Output (per sample):
#   S3: <S3_RNASEQ_OUT>/fusion/Arriba_<sample>.txt
#       <S3_RNASEQ_OUT>/fusion/Arriba_<sample>.discarded.txt
# =============================================================================

set -euo pipefail

source /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
ARRIBA="/home/ec2-user/miniforge3/bin/arriba"
ARRIBA_DB="/home/ec2-user/miniforge3/var/lib/arriba"
BLACKLIST="${ARRIBA_DB}/blacklist_mm10_GRCm38_v2.5.1.tsv.gz"
KNOWN_FUSIONS="${ARRIBA_DB}/known_fusions_mm10_GRCm38_v2.5.1.tsv.gz"
GENOME="/home/ec2-user/ref/mm10/mm10.fa"
GTF="/home/ec2-user/ref/mm10/mm10_chr.gtf"
STAR_INDEX="/home/ec2-user/ref/mm10/star_index"
STAR="/home/ec2-user/miniforge3/bin/STAR-avx2"
THREADS=8
PARALLEL=2   # STAR uses 8 threads each → 16 threads in flight on a 16-vCPU box

S3_FASTQ_BASE="${S3_ROOT}/data/fastq/rnaseq/Raw-GemOVCA-RNAseq"
S3_FUSION_OUT="${S3_RESULTS}/rnaseq/fusion"
TMP_FUSION="${TMP_DIR}/arriba"
mkdir -p "${TMP_FUSION}"

log "============================================================"
log "Arriba Fusion Calling  RUN_ID=${RUN_ID}"
log "  Arriba  : $(${ARRIBA} --version 2>&1 | head -1)"
log "  S3 out  : ${S3_FUSION_OUT}"
log "  Parallel: ${PARALLEL} concurrent jobs"
log "============================================================"

# ---------------------------------------------------------------------------
# Determine samples to process
# ---------------------------------------------------------------------------
if [[ $# -gt 0 ]]; then
    SAMPLES=("$@")
else
    SAMPLES=("${RNASEQ_SAMPLES[@]}")
fi

# ---------------------------------------------------------------------------
# Per-sample function
# ---------------------------------------------------------------------------
run_sample() {
    local SAMPLE="$1"
    local PATTERN="${RNASEQ_PATTERNS[${SAMPLE}]}"
    local SAMPLE_TMP="${TMP_FUSION}/${SAMPLE}"
    local S3_OUT="${S3_FUSION_OUT}/Arriba_${SAMPLE}.txt"

    # Skip if already done
    if aws s3 ls "${S3_OUT}" &>/dev/null; then
        log "[SKIP] ${SAMPLE}: Arriba output already on S3"
        return 0
    fi

    log "[START] ${SAMPLE} (pattern='${PATTERN}')"
    mkdir -p "${SAMPLE_TMP}"

    # -----------------------------------------------------------------------
    # 1. Discover + download FASTQs
    # -----------------------------------------------------------------------
    local ALL_KEYS
    ALL_KEYS=$(aws s3 ls --recursive "${S3_FASTQ_BASE}/" \
                | awk '{print $4}' \
                | grep "_R[12]_.*\.fastq\.gz$" \
                | grep -E "${PATTERN}" || true)

    if [[ -z "${ALL_KEYS}" ]]; then
        log "[WARN] ${SAMPLE}: no FASTQs found matching '${PATTERN}' — skipping"
        return 0
    fi

    local R1_KEYS R2_KEYS
    R1_KEYS=$(echo "${ALL_KEYS}" | grep "_R1_" | sort)
    R2_KEYS=$(echo "${ALL_KEYS}" | grep "_R2_" | sort)

    log "${SAMPLE}: downloading FASTQs..."
    local R1_LOCAL_LIST=() R2_LOCAL_LIST=()

    while IFS= read -r key; do
        local fname; fname=$(basename "${key}")
        local lpath="${SAMPLE_TMP}/${fname}"
        aws s3 cp "s3://neoantigen2026-rerun/${key}" "${lpath}" --quiet
        R1_LOCAL_LIST+=("${lpath}")
    done <<< "${R1_KEYS}"

    while IFS= read -r key; do
        local fname; fname=$(basename "${key}")
        local lpath="${SAMPLE_TMP}/${fname}"
        aws s3 cp "s3://neoantigen2026-rerun/${key}" "${lpath}" --quiet
        R2_LOCAL_LIST+=("${lpath}")
    done <<< "${R2_KEYS}"

    local R1_ARG R2_ARG
    R1_ARG=$(IFS=,; echo "${R1_LOCAL_LIST[*]}")
    R2_ARG=$(IFS=,; echo "${R2_LOCAL_LIST[*]}")

    # -----------------------------------------------------------------------
    # 2. STAR alignment with chimeric read flags (required by Arriba)
    # -----------------------------------------------------------------------
    local BAM="${SAMPLE_TMP}/Aligned.sortedByCoord.out.bam"

    if [[ -f "${BAM}" && -s "${BAM}" ]]; then
        log "${SAMPLE}: STAR BAM already present ($(du -sh "${BAM}" | cut -f1)) — skipping alignment"
    else
        log "${SAMPLE}: running STAR (chimeric mode)..."
    fi

    local T0=$SECONDS
    if [[ ! -f "${BAM}" || ! -s "${BAM}" ]]; then
    "${STAR}" \
        --runThreadN "${THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${R1_ARG}" "${R2_ARG}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMultimapNmax 50 \
        --alignSJoverhangMin 8 \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --runMode alignReads \
        --twopassMode Basic \
        --outFileNamePrefix "${SAMPLE_TMP}/" \
        --outBAMsortingThreadN 1 \
        --peOverlapNbasesMin 10 \
        --peOverlapMMp 0.1 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM SoftClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreMin 1 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50

    log "${SAMPLE}: STAR done in $(( SECONDS - T0 ))s."
    fi  # end STAR block

    if [[ ! -f "${BAM}" || ! -s "${BAM}" ]]; then
        log "[ERROR] ${SAMPLE}: BAM not found or empty after STAR — aborting"
        rm -rf "${SAMPLE_TMP}"
        return 1
    fi

    # -----------------------------------------------------------------------
    # 3. Run Arriba
    # -----------------------------------------------------------------------
    log "${SAMPLE}: running Arriba..."
    local FUSIONS_OUT="${SAMPLE_TMP}/Arriba_${SAMPLE}.txt"
    local DISCARDED_OUT="${SAMPLE_TMP}/Arriba_${SAMPLE}.discarded.txt"

    "${ARRIBA}" \
        -x "${BAM}" \
        -a "${GENOME}" \
        -g "${GTF}" \
        -b "${BLACKLIST}" \
        -k "${KNOWN_FUSIONS}" \
        -o "${FUSIONS_OUT}" \
        -O "${DISCARDED_OUT}"

    log "${SAMPLE}: Arriba complete. Fusions: $(wc -l < "${FUSIONS_OUT}") lines"

    # -----------------------------------------------------------------------
    # 4. Upload to S3
    # -----------------------------------------------------------------------
    log "${SAMPLE}: uploading to S3..."
    aws s3 cp "${FUSIONS_OUT}"    "${S3_FUSION_OUT}/Arriba_${SAMPLE}.txt"
    aws s3 cp "${DISCARDED_OUT}"  "${S3_FUSION_OUT}/Arriba_${SAMPLE}.discarded.txt"
    log "${SAMPLE}: uploaded to ${S3_FUSION_OUT}/"

    # -----------------------------------------------------------------------
    # 5. Clean up
    # -----------------------------------------------------------------------
    rm -rf "${SAMPLE_TMP}"
    log "[DONE] ${SAMPLE}"
}

# ---------------------------------------------------------------------------
# Main — run with parallelism
# ---------------------------------------------------------------------------
PIDS=()

for SAMPLE in "${SAMPLES[@]}"; do
    run_sample "${SAMPLE}" &
    PIDS+=($!)
    log "Dispatched ${SAMPLE} (PID ${PIDS[-1]})"

    while [[ ${#PIDS[@]} -ge ${PARALLEL} ]]; do
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                if wait "${PIDS[$i]}"; then
                    log "PID ${PIDS[$i]} finished successfully."
                else
                    log "[ERROR] PID ${PIDS[$i]} failed — continuing with remaining samples"
                fi
                unset "PIDS[$i]"
                PIDS=("${PIDS[@]}")
                break
            fi
        done
        sleep 5
    done
done

for PID in "${PIDS[@]}"; do
    if wait "${PID}"; then
        log "PID ${PID} finished successfully."
    else
        log "[ERROR] PID ${PID} failed."
        exit 1
    fi
done

log "============================================================"
log "Arriba fusion calling complete."
log "  Results: ${S3_FUSION_OUT}/"
log "============================================================"
