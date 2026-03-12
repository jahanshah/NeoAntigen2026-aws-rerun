#!/usr/bin/env bash
# =============================================================================
# R02_align.sh — Download FASTQs from S3 and align with STAR (2-pass)
#
# For each sample in RNASEQ_SAMPLES:
#   - Skip if final BAM already on S3
#   - Discover R1/R2 FASTQ paths on S3 using the fastq_pattern from the manifest
#   - Download FASTQs locally
#   - Run STAR alignment (2-pass, GeneCounts enabled)
#   - Index BAM with samtools
#   - Upload BAM, BAI, Log.final.out, ReadsPerGene.out.tab to S3
#   - Clean up local FASTQs and BAM
#
# Parallelism: STAR_ALIGN_PARALLEL samples are processed concurrently.
# =============================================================================

set -euo pipefail

source /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh

# FASTQ base path on S3 (directory containing timepoint subdirs)
S3_FASTQ_BASE="${S3_FASTQ}/Raw-GemOVCA-RNAseq"

# =============================================================================
# Helper: align one sample
# =============================================================================
align_sample() {
    local SAMPLE="$1"
    local PATTERN="${RNASEQ_PATTERNS[${SAMPLE}]}"

    local S3_BAM_OUT="${S3_RNASEQ_OUT}/bam/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    local LOCAL_SAMPLE_DIR="${RNASEQ_TMP}/${SAMPLE}"

    # ------------------------------------------------------------------
    # Skip if BAM already on S3
    # ------------------------------------------------------------------
    if aws s3 ls "${S3_BAM_OUT}" &>/dev/null; then
        log "[SKIP] ${SAMPLE}: BAM already on S3 (${S3_BAM_OUT})"
        return 0
    fi

    log "[START] ${SAMPLE}: pattern='${PATTERN}'"
    mkdir -p "${LOCAL_SAMPLE_DIR}"

    # ------------------------------------------------------------------
    # Discover FASTQ paths on S3 matching the pattern
    # ------------------------------------------------------------------
    log "${SAMPLE}: searching for FASTQs under ${S3_FASTQ_BASE}/ ..."
    local ALL_KEYS
    ALL_KEYS=$(aws s3 ls --recursive "${S3_FASTQ_BASE}/" \
                | awk '{print $4}' \
                | grep "_R[12]_.*\.fastq\.gz$" \
                | grep -E "${PATTERN}" || true)

    if [[ -z "${ALL_KEYS}" ]]; then
        log "[WARN] ${SAMPLE}: No FASTQs found matching pattern '${PATTERN}' — skipping."
        return 0
    fi

    # Separate R1 and R2 (sort for deterministic order)
    local R1_KEYS R2_KEYS
    R1_KEYS=$(echo "${ALL_KEYS}" | grep "_R1_" | sort)
    R2_KEYS=$(echo "${ALL_KEYS}" | grep "_R2_" | sort)

    if [[ -z "${R1_KEYS}" || -z "${R2_KEYS}" ]]; then
        log "[WARN] ${SAMPLE}: Could not find both R1 and R2 FASTQs — skipping."
        log "  Found keys: ${ALL_KEYS}"
        return 0
    fi

    # Check R1 and R2 counts match
    local R1_COUNT R2_COUNT
    R1_COUNT=$(echo "${R1_KEYS}" | wc -l)
    R2_COUNT=$(echo "${R2_KEYS}" | wc -l)
    if [[ "${R1_COUNT}" -ne "${R2_COUNT}" ]]; then
        log "[WARN] ${SAMPLE}: R1 count (${R1_COUNT}) != R2 count (${R2_COUNT}) — skipping."
        return 0
    fi

    # ------------------------------------------------------------------
    # Download FASTQs
    # ------------------------------------------------------------------
    log "${SAMPLE}: downloading ${R1_COUNT} R1+R2 FASTQ pairs..."
    local R1_LOCAL_LIST=() R2_LOCAL_LIST=()

    while IFS= read -r key; do
        local fname
        fname=$(basename "${key}")
        local local_path="${LOCAL_SAMPLE_DIR}/${fname}"
        aws s3 cp "s3://$(echo "${S3_FASTQ_BASE}" | sed 's|s3://||' | cut -d/ -f1)/$(echo "${key}")" \
            "${local_path}" --quiet || \
        aws s3 cp "${S3_ROOT%/data*}/$(echo "${key}")" "${local_path}" --quiet || \
        aws s3 cp "s3://neoantigen2026-rerun/${key}" "${local_path}" --quiet
        R1_LOCAL_LIST+=("${local_path}")
    done <<< "${R1_KEYS}"

    while IFS= read -r key; do
        local fname
        fname=$(basename "${key}")
        local local_path="${LOCAL_SAMPLE_DIR}/${fname}"
        aws s3 cp "s3://neoantigen2026-rerun/${key}" "${local_path}" --quiet
        R2_LOCAL_LIST+=("${local_path}")
    done <<< "${R2_KEYS}"

    # Combine multiple lanes with commas (STAR accepts comma-separated lists)
    local R1_ARG R2_ARG
    R1_ARG=$(IFS=,; echo "${R1_LOCAL_LIST[*]}")
    R2_ARG=$(IFS=,; echo "${R2_LOCAL_LIST[*]}")

    # ------------------------------------------------------------------
    # STAR 2-pass alignment
    # ------------------------------------------------------------------
    log "${SAMPLE}: running STAR alignment..."
    local STAR_OUTDIR="${LOCAL_SAMPLE_DIR}/"
    local T0=$SECONDS

    STAR \
        --runThreadN "${RNASEQ_THREADS}" \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${R1_ARG}" "${R2_ARG}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs RemoveNoncanonical \
        --runMode alignReads \
        --twopassMode Basic \
        --quantMode GeneCounts \
        --outFileNamePrefix "${STAR_OUTDIR}" \
        --outBAMsortingThreadN "${RNASEQ_THREADS}"

    log "${SAMPLE}: STAR alignment complete in $(( SECONDS - T0 ))s."

    local BAM="${STAR_OUTDIR}Aligned.sortedByCoord.out.bam"
    if [[ ! -f "${BAM}" ]]; then
        log "[ERROR] ${SAMPLE}: Expected BAM not found: ${BAM}"
        return 1
    fi

    # ------------------------------------------------------------------
    # Index BAM
    # ------------------------------------------------------------------
    log "${SAMPLE}: indexing BAM..."
    samtools index -@ "${RNASEQ_THREADS}" "${BAM}"

    # ------------------------------------------------------------------
    # Upload to S3
    # ------------------------------------------------------------------
    log "${SAMPLE}: uploading to S3..."
    aws s3 cp "${BAM}"                                            "${S3_RNASEQ_OUT}/bam/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    aws s3 cp "${BAM}.bai"                                        "${S3_RNASEQ_OUT}/bam/${SAMPLE}.Aligned.sortedByCoord.out.bam.bai"
    aws s3 cp "${STAR_OUTDIR}Log.final.out"                       "${S3_RNASEQ_OUT}/logs/${SAMPLE}.Log.final.out"
    aws s3 cp "${STAR_OUTDIR}ReadsPerGene.out.tab"                "${S3_RNASEQ_OUT}/counts/${SAMPLE}.ReadsPerGene.out.tab"
    log "${SAMPLE}: upload complete."

    # ------------------------------------------------------------------
    # Clean up local files to free disk space
    # ------------------------------------------------------------------
    log "${SAMPLE}: cleaning up local files..."
    rm -rf "${LOCAL_SAMPLE_DIR}"
    log "${SAMPLE}: done."
}

# =============================================================================
# Main — process all samples with STAR_ALIGN_PARALLEL concurrency
# =============================================================================
log "============================================================"
log "R02_align.sh — aligning ${#RNASEQ_SAMPLES[@]} samples"
log "  Parallelism : ${STAR_ALIGN_PARALLEL} concurrent jobs"
log "  STAR index  : ${STAR_INDEX}"
log "  S3 output   : ${S3_RNASEQ_OUT}/bam/"
log "============================================================"

PIDS=()

for SAMPLE in "${RNASEQ_SAMPLES[@]}"; do
    # Launch alignment in background
    align_sample "${SAMPLE}" &
    PIDS+=($!)
    log "Launched alignment for ${SAMPLE} (PID ${PIDS[-1]})"

    # Throttle: wait for a slot if we've reached STAR_ALIGN_PARALLEL
    while [[ ${#PIDS[@]} -ge ${STAR_ALIGN_PARALLEL} ]]; do
        # Wait for any one job to finish
        for i in "${!PIDS[@]}"; do
            if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                # Job finished — check exit code
                if wait "${PIDS[$i]}"; then
                    log "Job PID ${PIDS[$i]} completed successfully."
                else
                    log "[ERROR] Job PID ${PIDS[$i]} failed."
                    # Kill remaining jobs before exiting
                    for pid in "${PIDS[@]}"; do kill "${pid}" 2>/dev/null || true; done
                    exit 1
                fi
                unset "PIDS[$i]"
                PIDS=("${PIDS[@]}")  # re-index array
                break
            fi
        done
        sleep 5
    done
done

# Wait for any remaining jobs
for PID in "${PIDS[@]}"; do
    if wait "${PID}"; then
        log "Job PID ${PID} completed successfully."
    else
        log "[ERROR] Job PID ${PID} failed."
        exit 1
    fi
done

log "============================================================"
log "R02_align.sh complete — all samples aligned."
log "============================================================"
