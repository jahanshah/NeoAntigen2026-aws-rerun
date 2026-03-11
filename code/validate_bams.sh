#!/usr/bin/env bash
# =============================================================================
# validate_bams.sh
# Validates BAM files in s3://bam-wes/NeoAntigen-aws/data/bam/
# Checks: integrity (BGZF EOF), sort order, index presence, flagstat stats
# Outputs: TSV + per-sample flagstat -> parsed by summarize_bam_qc.R
# Strategy: streams BAMs from S3 via "aws s3 cp ... -" piped to samtools
# =============================================================================

set -euo pipefail

# --- Config ------------------------------------------------------------------
S3_BAM_PREFIX="s3://bam-wes/NeoAntigen-aws/data/bam"
OUTDIR="/home/ec2-user/results/bam_validation"
FLAGSTAT_DIR="${OUTDIR}/flagstat"
RESULTS_TSV="${OUTDIR}/bam_validation_results.tsv"
SAMTOOLS="/home/ec2-user/miniforge3/bin/samtools"
HEADER_BYTES=10485760   # 10 MB — enough to cover any BAM header

mkdir -p "${OUTDIR}" "${FLAGSTAT_DIR}"

# --- Verify tools ------------------------------------------------------------
for tool in "${SAMTOOLS}" aws python3; do
    command -v "${tool}" &>/dev/null || {
        echo "[ERROR] Required tool not found: ${tool}" >&2; exit 1; }
done
echo "[INFO] samtools: $(${SAMTOOLS} --version | head -1)"
echo "[INFO] AWS CLI:  $(aws --version 2>&1 | head -1)"

# --- Collect BAM list from S3 ------------------------------------------------
echo "[INFO] Listing BAM files in ${S3_BAM_PREFIX}/"
mapfile -t BAM_FILES < <(aws s3 ls "${S3_BAM_PREFIX}/" \
    | awk '{print $NF}' | grep '\.bam$')

[[ ${#BAM_FILES[@]} -eq 0 ]] && {
    echo "[ERROR] No BAM files found at ${S3_BAM_PREFIX}/" >&2; exit 1; }
echo "[INFO] Found ${#BAM_FILES[@]} BAM file(s)"

# --- TSV header --------------------------------------------------------------
printf "sample\ts3_path\tsize_gb\tquickcheck\tsort_order\treference_build\t\
read_group\tindex_present\ttotal_reads\tmapped_reads\tmapped_pct\t\
duplication_pct\tpassed_qc\tnotes\n" > "${RESULTS_TSV}"

# --- Helper: check BGZF EOF block (last 28 bytes) ----------------------------
# BGZF EOF magic = 1f8b0804000000000000ff0600424302001b0003000000000000000000
BGZF_EOF_HEX="1f8b0804000000000000ff0600424302001b0003000000000000000000"

check_bgzf_eof() {
    local s3_path="$1"
    local file_size="$2"
    local start_byte=$(( file_size - 28 ))
    # Use s3api to fetch last 28 bytes
    local eof_bytes
    eof_bytes=$(aws s3api get-object \
        --bucket bam-wes \
        --key "${s3_path#s3://bam-wes/}" \
        --range "bytes=${start_byte}-${file_size}" \
        /tmp/_bam_eof_check.bin 2>/dev/null \
        && xxd -p /tmp/_bam_eof_check.bin | tr -d '\n' || echo "error")
    rm -f /tmp/_bam_eof_check.bin
    if [[ "${eof_bytes}" == *"${BGZF_EOF_HEX}"* ]]; then
        echo "PASS"
    else
        echo "FAIL"
    fi
}

# --- Helper: stream first N bytes from S3 to a temp file --------------------
fetch_header_bytes() {
    local s3_path="$1"
    local n_bytes="$2"
    local outfile="$3"
    local file_size="$4"
    local end=$(( n_bytes < file_size ? n_bytes - 1 : file_size - 1 ))
    aws s3api get-object \
        --bucket bam-wes \
        --key "${s3_path#s3://bam-wes/}" \
        --range "bytes=0-${end}" \
        "${outfile}" > /dev/null 2>&1
}

# --- Validate each BAM -------------------------------------------------------
for BAM_FILE in "${BAM_FILES[@]}"; do
    S3_PATH="${S3_BAM_PREFIX}/${BAM_FILE}"
    SAMPLE="${BAM_FILE%.bam}"
    INDEX_PATH="${S3_PATH}.bai"

    echo ""
    echo "========================================"
    echo "[INFO] Validating: ${BAM_FILE}"
    echo "========================================"

    NOTES=""
    PASSED_QC="PASS"

    # 1. File size from S3 metadata
    SIZE_BYTES=$(aws s3 ls "${S3_PATH}" | awk '{print $3}')
    SIZE_GB=$(awk "BEGIN {printf \"%.2f\", ${SIZE_BYTES}/1073741824}")
    echo "  Size: ${SIZE_GB} GB"

    # 2. BGZF EOF integrity check (last 28 bytes)
    echo "  Checking BGZF EOF integrity..."
    QUICKCHECK=$(check_bgzf_eof "${S3_PATH}" "${SIZE_BYTES}")
    if [[ "${QUICKCHECK}" == "FAIL" ]]; then
        PASSED_QC="FAIL"
        NOTES="${NOTES}bgzf_eof_missing;"
        echo "  [WARN] BGZF EOF check FAILED — file may be truncated"
    else
        echo "  BGZF EOF: PASS"
    fi

    # 3. Fetch first 10 MB to read BAM header
    echo "  Fetching BAM header (first ${HEADER_BYTES} bytes)..."
    HEADER_TMP="/tmp/${SAMPLE}_header.bam"
    fetch_header_bytes "${S3_PATH}" "${HEADER_BYTES}" "${HEADER_TMP}" "${SIZE_BYTES}"

    HEADER=$("${SAMTOOLS}" view -H "${HEADER_TMP}" 2>/dev/null) || {
        echo "  [ERROR] Cannot read header — skipping"
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${SAMPLE}" "${S3_PATH}" "${SIZE_GB}" "${QUICKCHECK}" \
            "unknown" "unknown" "unknown" "unknown" \
            "NA" "NA" "NA" "NA" "FAIL" "header_read_error" \
            >> "${RESULTS_TSV}"
        rm -f "${HEADER_TMP}"
        continue
    }
    rm -f "${HEADER_TMP}"

    # Sort order
    SORT_ORDER=$(echo "${HEADER}" | awk '/^@HD/{
        for(i=1;i<=NF;i++) if($i~/^SO:/) {sub("SO:","",$i); print $i}
    }')
    SORT_ORDER="${SORT_ORDER:-unknown}"
    if [[ "${SORT_ORDER}" != "coordinate" ]]; then
        PASSED_QC="FAIL"
        NOTES="${NOTES}not_coordinate_sorted(${SORT_ORDER});"
        echo "  [WARN] Sort order: ${SORT_ORDER} (expected: coordinate)"
    else
        echo "  Sort order: coordinate [OK]"
    fi

    # Reference build
    REF_BUILD="unknown"
    SQ_NAMES=$(echo "${HEADER}" | awk '/^@SQ/{for(i=1;i<=NF;i++) if($i~/^SN:/) print $i}' | head -5)
    if echo "${HEADER}" | grep -qiE "mm10|GRCm38|mus_musculus"; then
        REF_BUILD="mm10/GRCm38"
    elif echo "${HEADER}" | grep -qiE "hg38|GRCh38"; then
        REF_BUILD="hg38/GRCh38"
    elif echo "${HEADER}" | grep -qiE "hg19|GRCh37"; then
        REF_BUILD="hg19/GRCh37"
    elif echo "${SQ_NAMES}" | grep -qiE "chr[0-9XYM]"; then
        REF_BUILD="ucsc_style"
    fi
    echo "  Reference build: ${REF_BUILD}"

    # Read groups
    RG_PRESENT=$(echo "${HEADER}" | grep -c "^@RG" || true)
    READ_GROUP=$([ "${RG_PRESENT}" -gt 0 ] && echo "YES(${RG_PRESENT})" || echo "NO")
    if [[ "${RG_PRESENT}" -eq 0 ]]; then
        NOTES="${NOTES}no_read_group;"
        echo "  [WARN] No @RG read group tags"
    else
        echo "  Read groups: ${READ_GROUP}"
    fi

    # 4. BAM index check
    echo "  Checking for BAM index (.bai)..."
    INDEX_PRESENT="NO"
    if aws s3 ls "${INDEX_PATH}" &>/dev/null \
       || aws s3 ls "${S3_BAM_PREFIX}/${SAMPLE}.bai" &>/dev/null; then
        INDEX_PRESENT="YES"
        echo "  Index (.bai): present"
    else
        NOTES="${NOTES}no_bai_index;"
        echo "  [WARN] No .bai index found"
    fi

    # 5. samtools flagstat — stream full BAM from S3 via pipe
    echo "  Running flagstat (streaming full BAM from S3)..."
    FLAGSTAT_FILE="${FLAGSTAT_DIR}/${SAMPLE}.flagstat.txt"
    if aws s3 cp "${S3_PATH}" - 2>/dev/null \
        | "${SAMTOOLS}" flagstat - > "${FLAGSTAT_FILE}" 2>/dev/null; then
        TOTAL_READS=$(awk '/in total/{print $1}'       "${FLAGSTAT_FILE}")
        MAPPED_READS=$(awk '/mapped \(/{print $1}'     "${FLAGSTAT_FILE}" | head -1)
        DUP_READS=$(awk '/duplicates/{print $1}'       "${FLAGSTAT_FILE}" | head -1)
        if [[ -n "${TOTAL_READS}" && "${TOTAL_READS}" -gt 0 ]]; then
            MAPPED_PCT=$(awk "BEGIN {printf \"%.2f\", 100*${MAPPED_READS}/${TOTAL_READS}}")
            DUP_PCT=$(awk "BEGIN {printf \"%.2f\", 100*${DUP_READS:-0}/${TOTAL_READS}}")
        else
            MAPPED_PCT="NA"; DUP_PCT="NA"
        fi
        echo "  Total reads  : ${TOTAL_READS}"
        echo "  Mapped       : ${MAPPED_PCT}%"
        echo "  Duplicates   : ${DUP_PCT}%"

        if [[ "${MAPPED_PCT}" != "NA" ]] \
           && (( $(echo "${MAPPED_PCT} < 80" | bc -l 2>/dev/null || echo 0) )); then
            PASSED_QC="WARN"
            NOTES="${NOTES}low_mapping_rate(${MAPPED_PCT}%);"
        fi
    else
        TOTAL_READS="NA"; MAPPED_READS="NA"; MAPPED_PCT="NA"; DUP_PCT="NA"
        NOTES="${NOTES}flagstat_failed;"
        echo "  [WARN] flagstat failed"
    fi

    echo "  QC status: ${PASSED_QC}"

    # 6. Write TSV row
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SAMPLE}" "${S3_PATH}" "${SIZE_GB}" "${QUICKCHECK}" \
        "${SORT_ORDER}" "${REF_BUILD}" "${READ_GROUP}" "${INDEX_PRESENT}" \
        "${TOTAL_READS}" "${MAPPED_READS}" "${MAPPED_PCT}" "${DUP_PCT}" \
        "${PASSED_QC}" "${NOTES:-none}" \
        >> "${RESULTS_TSV}"
done

echo ""
echo "========================================"
echo "[INFO] Validation complete."
echo "[INFO] Results TSV : ${RESULTS_TSV}"
echo "[INFO] Flagstat dir: ${FLAGSTAT_DIR}/"
echo "========================================"

# --- Run R summary -----------------------------------------------------------
echo "[INFO] Generating QC report with R..."
if command -v Rscript &>/dev/null; then
    Rscript /home/ec2-user/code/summarize_bam_qc.R "${RESULTS_TSV}" "${OUTDIR}" 2>&1
else
    echo "[WARN] Rscript not found. Run summarize_bam_qc.R manually."
fi
