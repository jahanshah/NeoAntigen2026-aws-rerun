#!/usr/bin/env bash
# =============================================================================
# index_bams.sh
# Downloads each BAM from S3 one at a time, runs samtools index,
# uploads the .bai back to S3, then removes the local BAM to save disk.
# =============================================================================

set -euo pipefail

# --- Config ------------------------------------------------------------------
S3_BAM_PREFIX="s3://bam-wes/NeoAntigen-aws/data/bam"
TMPDIR="/tmp/bam_index"
SAMTOOLS="/home/ec2-user/miniforge3/bin/samtools"
RESULTS_TSV="/home/ec2-user/results/bam_validation/bam_index_results.tsv"
THREADS=4

mkdir -p "${TMPDIR}"

# --- Verify tools ------------------------------------------------------------
for tool in "${SAMTOOLS}" aws; do
    command -v "${tool}" &>/dev/null || {
        echo "[ERROR] Required tool not found: ${tool}" >&2; exit 1; }
done
echo "[INFO] samtools: $(${SAMTOOLS} --version | head -1)"

# --- Collect BAM list --------------------------------------------------------
mapfile -t BAM_FILES < <(aws s3 ls "${S3_BAM_PREFIX}/" \
    | awk '{print $NF}' | grep '\.bam$')

echo "[INFO] Found ${#BAM_FILES[@]} BAM file(s) to index"
echo "[INFO] Processing one at a time to manage disk space"

# --- TSV header --------------------------------------------------------------
printf "sample\ts3_bam\ts3_bai\tsize_gb\tindex_status\tindex_time_min\tnotes\n" \
    > "${RESULTS_TSV}"

# --- Index each BAM ----------------------------------------------------------
for BAM_FILE in "${BAM_FILES[@]}"; do
    S3_BAM="${S3_BAM_PREFIX}/${BAM_FILE}"
    S3_BAI="${S3_BAM}.bai"
    SAMPLE="${BAM_FILE%.bam}"
    LOCAL_BAM="${TMPDIR}/${BAM_FILE}"
    LOCAL_BAI="${LOCAL_BAM}.bai"

    echo ""
    echo "========================================"
    echo "[INFO] Processing: ${BAM_FILE}"
    echo "========================================"

    STATUS="PASS"
    NOTES="none"
    START_TIME=$(date +%s)

    # Check if .bai already exists in S3
    if aws s3 ls "${S3_BAI}" &>/dev/null; then
        echo "  [SKIP] Index already exists: ${S3_BAI}"
        SIZE_BYTES=$(aws s3 ls "${S3_BAM}" | awk '{print $3}')
        SIZE_GB=$(awk "BEGIN {printf \"%.2f\", ${SIZE_BYTES}/1073741824}")
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${SAMPLE}" "${S3_BAM}" "${S3_BAI}" "${SIZE_GB}" \
            "SKIPPED" "0" "already_indexed" >> "${RESULTS_TSV}"
        continue
    fi

    # Get file size
    SIZE_BYTES=$(aws s3 ls "${S3_BAM}" | awk '{print $3}')
    SIZE_GB=$(awk "BEGIN {printf \"%.2f\", ${SIZE_BYTES}/1073741824}")
    echo "  Size: ${SIZE_GB} GB"

    # Check disk space (need SIZE_GB + 0.5 GB buffer)
    AVAIL_GB=$(df /tmp | awk 'NR==2 {printf "%.1f", $4/1048576}')
    NEEDED_GB=$(awk "BEGIN {printf \"%.1f\", ${SIZE_GB} + 0.5}")
    echo "  Disk available: ${AVAIL_GB} GB  |  Needed: ${NEEDED_GB} GB"
    if (( $(echo "${AVAIL_GB} < ${NEEDED_GB}" | bc -l) )); then
        echo "  [ERROR] Insufficient disk space — skipping"
        STATUS="FAIL"; NOTES="insufficient_disk"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${SAMPLE}" "${S3_BAM}" "${S3_BAI}" "${SIZE_GB}" \
            "${STATUS}" "NA" "${NOTES}" >> "${RESULTS_TSV}"
        continue
    fi

    # Download BAM
    echo "  Downloading from S3..."
    if ! aws s3 cp "${S3_BAM}" "${LOCAL_BAM}" 2>&1; then
        echo "  [ERROR] Download failed"
        STATUS="FAIL"; NOTES="download_failed"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${SAMPLE}" "${S3_BAM}" "${S3_BAI}" "${SIZE_GB}" \
            "${STATUS}" "NA" "${NOTES}" >> "${RESULTS_TSV}"
        continue
    fi
    echo "  Download complete."

    # samtools index
    echo "  Running samtools index (threads: ${THREADS})..."
    if ! "${SAMTOOLS}" index -@ ${THREADS} "${LOCAL_BAM}" 2>&1; then
        echo "  [ERROR] samtools index failed"
        STATUS="FAIL"; NOTES="index_failed"
        rm -f "${LOCAL_BAM}" "${LOCAL_BAI}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${SAMPLE}" "${S3_BAM}" "${S3_BAI}" "${SIZE_GB}" \
            "${STATUS}" "NA" "${NOTES}" >> "${RESULTS_TSV}"
        continue
    fi
    BAI_SIZE=$(ls -lh "${LOCAL_BAI}" | awk '{print $5}')
    echo "  Index created: ${LOCAL_BAI} (${BAI_SIZE})"

    # Upload .bai to S3
    echo "  Uploading .bai to S3..."
    if ! aws s3 cp "${LOCAL_BAI}" "${S3_BAI}" 2>&1; then
        echo "  [ERROR] Upload failed"
        STATUS="FAIL"; NOTES="upload_failed"
        rm -f "${LOCAL_BAM}" "${LOCAL_BAI}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${SAMPLE}" "${S3_BAM}" "${S3_BAI}" "${SIZE_GB}" \
            "${STATUS}" "NA" "${NOTES}" >> "${RESULTS_TSV}"
        continue
    fi
    echo "  Uploaded: ${S3_BAI}"

    # Clean up local files
    rm -f "${LOCAL_BAM}" "${LOCAL_BAI}"
    echo "  Local files removed."

    END_TIME=$(date +%s)
    ELAPSED_MIN=$(awk "BEGIN {printf \"%.1f\", (${END_TIME}-${START_TIME})/60}")
    echo "  Done in ${ELAPSED_MIN} min  |  QC: ${STATUS}"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SAMPLE}" "${S3_BAM}" "${S3_BAI}" "${SIZE_GB}" \
        "${STATUS}" "${ELAPSED_MIN}" "${NOTES}" >> "${RESULTS_TSV}"
done

echo ""
echo "========================================"
echo "[INFO] Indexing complete."
echo "[INFO] Results: ${RESULTS_TSV}"
echo "========================================"

# Verify all .bai files are now present in S3
echo ""
echo "[INFO] Verifying .bai files in S3:"
for BAM_FILE in "${BAM_FILES[@]}"; do
    S3_BAI="${S3_BAM_PREFIX}/${BAM_FILE}.bai"
    if aws s3 ls "${S3_BAI}" &>/dev/null; then
        echo "  ✓ ${BAM_FILE}.bai"
    else
        echo "  ✗ MISSING: ${BAM_FILE}.bai"
    fi
done
