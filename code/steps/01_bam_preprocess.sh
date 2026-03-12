#!/usr/bin/env bash
# =============================================================================
# Step 1 — BAM Preprocessing
# Sort → MarkDuplicates → BaseQualityScoreRecalibration → ValidateSamFile
#
# Preprocessed BAMs are stored permanently at:
#   s3://neoantigen2026-rerun/data/bam/wes/preprocessed/
# and reused across runs without reprocessing.
#
# Each run also gets a reference copy at:
#   s3://neoantigen2026-rerun/results/<RUN_ID>/preprocessed_bams/
# (S3 copy — no extra storage cost for the run reference)
#
# Parallel processing:
#   Set PARALLEL_JOBS=3 to process 3 samples concurrently.
#   Requires /scratch EBS volume (≥300GB). See setup_scratch_volume.sh.
#   Default: PARALLEL_JOBS=1 (sequential, safe on /tmp tmpfs).
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

OUTDIR="${RESULTS_DIR}/preprocessed_bams"
S3_PERMANENT="${S3_BAM_PREPROC}"          # permanent store — reused across runs
S3_RUN_PREPROC="${S3_RESULTS}/preprocessed_bams"   # per-run reference copy
mkdir -p "${OUTDIR}" "${TMP_DIR}/preproc"

# Parallel jobs: 1=sequential (safe with /tmp); 3=parallel (needs /scratch)
PARALLEL_JOBS="${PARALLEL_JOBS:-1}"
if [[ "${PARALLEL_JOBS}" -gt 1 ]] && [[ "${TMP_DIR}" == /tmp/* ]]; then
    log "[WARN] PARALLEL_JOBS=${PARALLEL_JOBS} requested but TMP_DIR is on /tmp (tmpfs)."
    log "[WARN] This may cause disk-full errors. Mount /scratch first (see setup_scratch_volume.sh)."
    log "[WARN] Falling back to PARALLEL_JOBS=1."
    PARALLEL_JOBS=1
fi

# SAMPLE_FILTER: space-separated list to restrict which samples to process.
# Set by run_pipeline.sh for phased execution (normals first, then tumors).
# Unset or empty = process all samples (default interactive behaviour).
if [[ -n "${SAMPLE_FILTER:-}" ]]; then
    read -ra ALL_SAMPLES <<< "${SAMPLE_FILTER}"
    log "Sample filter active: ${ALL_SAMPLES[*]}"
else
    ALL_SAMPLES=("${NORMAL_SAMPLE}" "${PON_NORMALS[1]}" "${TUMOR_SAMPLES[@]}")
fi

# NON_INTERACTIVE: when set, auto-answer prompts without user input.
#   done_permanent → copy_only   (fast S3-side copy to run results)
#   done_run / done_local → no   (already accessible, skip)
#   new → yes                    (always process)
NON_INTERACTIVE="${NON_INTERACTIVE:-}"

# --- Pre-run status check ----------------------------------------------------
log "========================================"
log "Checking preprocessing status for all samples..."
log "========================================"
declare -A SAMPLE_STATUS   # "done_permanent" | "done_run" | "done_local" | "new"
declare -A SAMPLE_INFO

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    S3_PERM_OUT="${S3_PERMANENT}/${SAMPLE}.preproc.bam"
    S3_RUN_OUT="${S3_RUN_PREPROC}/${SAMPLE}.preproc.bam"
    LOCAL_BAM="${OUTDIR}/${SAMPLE}.preproc.bam"

    if S3_META=$(aws s3 ls "${S3_PERM_OUT}" 2>/dev/null); then
        S3_DATE=$(echo "${S3_META}" | awk '{print $1, $2}')
        SAMPLE_STATUS["${SAMPLE}"]="done_permanent"
        SAMPLE_INFO["${SAMPLE}"]="preprocessed/indexed in permanent S3 store (${S3_DATE})"
    elif S3_META=$(aws s3 ls "${S3_RUN_OUT}" 2>/dev/null); then
        S3_DATE=$(echo "${S3_META}" | awk '{print $1, $2}')
        SAMPLE_STATUS["${SAMPLE}"]="done_run"
        SAMPLE_INFO["${SAMPLE}"]="preprocessed in run results (${S3_DATE})"
    elif [[ -f "${LOCAL_BAM}" ]]; then
        LOCAL_DATE=$(date -r "${LOCAL_BAM}" '+%Y-%m-%d %H:%M')
        SAMPLE_STATUS["${SAMPLE}"]="done_local"
        SAMPLE_INFO["${SAMPLE}"]="preprocessed locally (${LOCAL_DATE}), not yet on S3"
    else
        SAMPLE_STATUS["${SAMPLE}"]="new"
        SAMPLE_INFO["${SAMPLE}"]="not yet processed"
    fi
    log "  ${SAMPLE}: ${SAMPLE_INFO[${SAMPLE}]}"
done
log "========================================"

# --- Per-sample interactive confirmation -------------------------------------
declare -A RUN_SAMPLE   # "yes" | "no" | "copy_only"

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    if [[ "${SAMPLE_STATUS[${SAMPLE}]}" == "new" ]]; then
        RUN_SAMPLE["${SAMPLE}"]="yes"
    elif [[ "${SAMPLE_STATUS[${SAMPLE}]}" == "done_permanent" ]]; then
        if [[ -n "${NON_INTERACTIVE}" ]]; then
            RUN_SAMPLE["${SAMPLE}"]="copy_only"
            log "  [auto] ${SAMPLE} → copy to run results"
        else
            echo ""
            echo "  Sample '${SAMPLE}' was previously ${SAMPLE_INFO[${SAMPLE}]}."
            while true; do
                read -r -p "  Skip / copy to run results / reprocess? [s=skip / c=copy / r=reprocess]: " CHOICE
                case "${CHOICE}" in
                    s|S) RUN_SAMPLE["${SAMPLE}"]="no";        log "  -> Skipping ${SAMPLE}";              break ;;
                    c|C) RUN_SAMPLE["${SAMPLE}"]="copy_only"; log "  -> Will copy ${SAMPLE} to run results"; break ;;
                    r|R) RUN_SAMPLE["${SAMPLE}"]="yes";       log "  -> Will reprocess ${SAMPLE}";         break ;;
                    *) echo "  Please enter 's', 'c', or 'r'." ;;
                esac
            done
        fi
    else
        if [[ -n "${NON_INTERACTIVE}" ]]; then
            RUN_SAMPLE["${SAMPLE}"]="no"
            log "  [auto] ${SAMPLE} → skip (already exists)"
        else
            echo ""
            echo "  Sample '${SAMPLE}' was previously ${SAMPLE_INFO[${SAMPLE}]}."
            while true; do
                read -r -p "  Skip or reprocess? [s=skip / r=reprocess]: " CHOICE
                case "${CHOICE}" in
                    s|S) RUN_SAMPLE["${SAMPLE}"]="no";  log "  -> Skipping ${SAMPLE}";    break ;;
                    r|R) RUN_SAMPLE["${SAMPLE}"]="yes"; log "  -> Will reprocess ${SAMPLE}"; break ;;
                    *) echo "  Please enter 's' or 'r'." ;;
                esac
            done
        fi
    fi
done
echo ""

# =============================================================================
# preprocess_one_sample: run all preprocessing steps for a single sample.
# Called sequentially or in background depending on PARALLEL_JOBS.
# =============================================================================
preprocess_one_sample() {
    local SAMPLE="$1"
    local FINAL_BAM="${OUTDIR}/${SAMPLE}.preproc.bam"
    local S3_PERM_OUT="${S3_PERMANENT}/${SAMPLE}.preproc.bam"
    local S3_RUN_OUT="${S3_RUN_PREPROC}/${SAMPLE}.preproc.bam"

    local LOCAL_RAW="${TMP_DIR}/preproc/${SAMPLE}.bam"
    local LOCAL_SORTED="${TMP_DIR}/preproc/${SAMPLE}.sorted.bam"
    local LOCAL_DEDUP="${TMP_DIR}/preproc/${SAMPLE}.dedup.bam"
    local METRICS="${OUTDIR}/${SAMPLE}.dup_metrics.txt"
    local RECAL_TABLE="${OUTDIR}/${SAMPLE}.recal.table"

    log "========================================"
    log "Preprocessing: ${SAMPLE}"
    log "========================================"

    # --- Disk space check ----------------------------------------------------
    local AVAIL_GB
    AVAIL_GB=$(df -BG "${TMP_DIR}" | awk 'NR==2{gsub("G","",$4); print $4}')
    if [[ "${AVAIL_GB}" -lt 18 ]]; then
        log "[ERROR] Only ${AVAIL_GB}GB free in ${TMP_DIR} — need ≥18GB for ${SAMPLE}. Aborting."
        return 1
    fi
    log "  Disk: ${AVAIL_GB}GB free in ${TMP_DIR}"

    # --- 1a. Download raw BAM ------------------------------------------------
    log "  Downloading from S3..."
    aws s3 cp "${S3_BAM}/${SAMPLE}.bam"     "${LOCAL_RAW}"
    aws s3 cp "${S3_BAM}/${SAMPLE}.bam.bai" "${LOCAL_RAW}.bai"

    # --- 1b. Sort (coordinate) -----------------------------------------------
    local SORT_ORDER
    SORT_ORDER=$(${SAMTOOLS} view -H "${LOCAL_RAW}" | awk '/^@HD/{
        for(i=1;i<=NF;i++) if($i~/^SO:/) {sub("SO:","",$i); print $i}}')
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

    # --- 1c. MarkDuplicates --------------------------------------------------
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
            --TMP_DIR "${TMP_DIR}/preproc"
        [[ ! -s "${LOCAL_DEDUP}" ]] && { log "[ERROR] MarkDuplicates produced no output for ${SAMPLE}"; return 1; }
        log "  Dup metrics: ${METRICS}"
    fi
    rm -f "${LOCAL_SORTED}" "${LOCAL_SORTED}.bai"

    # --- 1d. BQSR (scattered across chromosomes, then merged) ----------------
    local KNOWN_SITES_ARG=""
    if [[ -f "${KNOWN_SNPS}" ]]; then
        KNOWN_SITES_ARG="--known-sites ${KNOWN_SNPS}"
    else
        log "  [WARN] No known SNPs — BQSR running without dbSNP"
        KNOWN_SITES_ARG="--run-without-dbsnp-potentially-ruining-quality"
    fi

    # Main chromosomes only (skips unplaced contigs for speed; captures >99% of reads)
    local CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
                  chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19
                  chrX chrY chrM)
    local SCATTER_DIR="${TMP_DIR}/preproc/${SAMPLE}_bqsr_scatter"
    mkdir -p "${SCATTER_DIR}"

    log "  Running BaseRecalibrator (scattered across ${#CHROMS[@]} chromosomes)..."
    local SCATTER_PIDS=() SCATTER_TABLES=()
    local BQSR_JAVA="-Xmx1500m"   # small heap per scatter job; I/O-bound not memory-bound

    for CHR in "${CHROMS[@]}"; do
        local CHR_TABLE="${SCATTER_DIR}/${SAMPLE}.${CHR}.recal.table"
        SCATTER_TABLES+=("${CHR_TABLE}")
        ${GATK} --java-options "${BQSR_JAVA}" BaseRecalibrator \
            -R "${REF}" -I "${LOCAL_DEDUP}" \
            -L "${CHR}" \
            ${KNOWN_SITES_ARG} \
            -O "${CHR_TABLE}" > "${SCATTER_DIR}/${CHR}.log" 2>&1 &
        SCATTER_PIDS+=($!)
    done

    # Wait for all scatter jobs; fail fast on any error
    local SCATTER_FAILED=0
    for i in "${!SCATTER_PIDS[@]}"; do
        if ! wait "${SCATTER_PIDS[$i]}"; then
            log "  [ERROR] BaseRecalibrator failed for ${CHROMS[$i]} — see ${SCATTER_DIR}/${CHROMS[$i]}.log"
            SCATTER_FAILED=1
        fi
    done
    [[ "${SCATTER_FAILED}" -eq 1 ]] && return 1

    log "  Merging BQSR scatter tables..."
    local INPUT_ARGS=""
    for T in "${SCATTER_TABLES[@]}"; do INPUT_ARGS="${INPUT_ARGS} -I ${T}"; done
    ${GATK} --java-options "${JAVA_OPTS}" GatherBQSRReports \
        ${INPUT_ARGS} -O "${RECAL_TABLE}" 2>&1 | tail -2 || true
    rm -rf "${SCATTER_DIR}"

    log "  Applying BQSR..."
    ${GATK} ${JAVA_OPTS:+--java-options "${JAVA_OPTS}"} ApplyBQSR \
        -R "${REF}" \
        -I "${LOCAL_DEDUP}" \
        --bqsr-recal-file "${RECAL_TABLE}" \
        -O "${FINAL_BAM}"
    rm -f "${LOCAL_DEDUP}" "${LOCAL_DEDUP}.bai" "${LOCAL_DEDUP%.bam}.bai"

    # Verify output BAM was actually created
    if [[ ! -s "${FINAL_BAM}" ]]; then
        log "[ERROR] ApplyBQSR produced no output for ${SAMPLE} — aborting"
        return 1
    fi

    # --- 1e. ValidateSamFile -------------------------------------------------
    log "  Validating BAM..."
    local VALIDATE_OUT="${OUTDIR}/${SAMPLE}.validate.txt"
    ${GATK} ${JAVA_OPTS:+--java-options "${JAVA_OPTS}"} ValidateSamFile \
        -I "${FINAL_BAM}" \
        -O "${VALIDATE_OUT}" \
        --MODE SUMMARY 2>&1 | tail -3 || true
    if grep -q "No errors found" "${VALIDATE_OUT}"; then
        log "  Validation: PASS"
    else
        log "  [WARN] Validation issues — see ${VALIDATE_OUT}"
    fi

    # --- 1f. Upload to permanent S3 store + run results ----------------------
    log "  Uploading to permanent S3 store (${S3_PERMANENT})..."
    aws s3 cp "${FINAL_BAM}"          "${S3_PERM_OUT}"
    aws s3 cp "${FINAL_BAM%.bam}.bai" "${S3_PERM_OUT%.bam}.bai" 2>/dev/null || \
    aws s3 cp "${FINAL_BAM}.bai"      "${S3_PERM_OUT}.bai"      2>/dev/null || true
    aws s3 cp "${METRICS}"            "${S3_PERMANENT}/$(basename ${METRICS})" 2>/dev/null || true
    log "  Copying to run results (${S3_RUN_PREPROC})..."
    aws s3 cp "${S3_PERM_OUT}"              "${S3_RUN_OUT}"
    aws s3 cp "${S3_PERM_OUT%.bam}.bai"     "${S3_RUN_OUT%.bam}.bai" 2>/dev/null || \
    aws s3 cp "${S3_PERM_OUT}.bai"          "${S3_RUN_OUT}.bai"      2>/dev/null || true

    rm -f "${FINAL_BAM}" "${FINAL_BAM%.bam}.bai" "${FINAL_BAM}.bai"
    log "  Done: ${SAMPLE}"
}

# =============================================================================
# Main dispatch loop — sequential or parallel depending on PARALLEL_JOBS
# =============================================================================
PIDS=()
FAILED_SAMPLES=()

for SAMPLE in "${ALL_SAMPLES[@]}"; do
    FINAL_BAM="${OUTDIR}/${SAMPLE}.preproc.bam"
    S3_PERM_OUT="${S3_PERMANENT}/${SAMPLE}.preproc.bam"
    S3_RUN_OUT="${S3_RUN_PREPROC}/${SAMPLE}.preproc.bam"

    # --- Copy from permanent store to run results (no reprocessing) ----------
    if [[ "${RUN_SAMPLE[${SAMPLE}]}" == "copy_only" ]]; then
        log "Copying ${SAMPLE} from permanent store to run results..."
        aws s3 cp "${S3_PERM_OUT}"              "${S3_RUN_OUT}"
        aws s3 cp "${S3_PERM_OUT%.bam}.bai"     "${S3_RUN_OUT%.bam}.bai" 2>/dev/null || true
        aws s3 cp "${S3_PERM_OUT%.bam}.bam.bai" "${S3_RUN_OUT%.bam}.bam.bai" 2>/dev/null || true
        log "  Done: ${SAMPLE} (copied)"
        continue
    fi

    if [[ "${RUN_SAMPLE[${SAMPLE}]}" == "no" ]]; then
        skip "${SAMPLE} (user chose to skip)"
        continue
    fi

    if [[ "${PARALLEL_JOBS}" -gt 1 ]]; then
        # Parallel: wait if at job limit
        while [[ ${#PIDS[@]} -ge ${PARALLEL_JOBS} ]]; do
            for i in "${!PIDS[@]}"; do
                if ! kill -0 "${PIDS[$i]}" 2>/dev/null; then
                    wait "${PIDS[$i]}" || FAILED_SAMPLES+=("${PIDS[$i]}")
                    unset 'PIDS[$i]'
                    PIDS=("${PIDS[@]}")
                fi
            done
            sleep 5
        done
        preprocess_one_sample "${SAMPLE}" &
        PIDS+=($!)
    else
        # Sequential
        preprocess_one_sample "${SAMPLE}"
    fi
done

# Wait for any remaining background jobs
for PID in "${PIDS[@]}"; do
    wait "${PID}" || FAILED_SAMPLES+=("${PID}")
done

if [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
    log "[ERROR] The following samples failed: ${FAILED_SAMPLES[*]}"
    exit 1
fi

log "Step 1 complete."
log "  Permanent store : ${S3_PERMANENT}/"
log "  Run results     : ${S3_RUN_PREPROC}/"
