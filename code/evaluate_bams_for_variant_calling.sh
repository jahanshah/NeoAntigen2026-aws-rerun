#!/usr/bin/env bash
# =============================================================================
# evaluate_bams_for_variant_calling.sh
#
# Evaluates BAMs in s3://bam-wes/NeoAntigen-aws/data/bam/ for suitability
# with GATK Mutect2 and VarScan somatic variant calling.
#
# Checks beyond basic validation:
#   - @RG completeness (SM, ID, LB, PL) — required by Mutect2
#   - Contig naming convention (chr prefix) — must match reference
#   - Estimated coverage (WES target)
#   - Properly paired % (>= 90%)
#   - Duplicate marking status
#   - BAM index present (.bai)
#   - 443_D21_new missing index warning
#   - Per-tool readiness: Mutect2 | VarScan
#
# Output: TSV + per-sample report printed to stdout
# Reuses existing flagstat files from results/bam_validation/flagstat/
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

# --- Parameters --------------------------------------------------------------
S3_BAM_PREFIX="${S3_BAM}"
OUTDIR="${RESULTS_DIR}/bam_evaluation"
FLAGSTAT_DIR="${RESULTS_DIR}/bam_validation/flagstat"   # reuse existing flagstats
RESULTS_TSV="${OUTDIR}/variant_calling_readiness.tsv"

# WES coverage estimation parameters (mm10)
EXOME_SIZE_BP=45000000      # ~45 Mb mouse WES capture
READ_LENGTH=150             # typical paired-end read length
MIN_COVERAGE_TUMOR=30       # Mutect2 recommended minimum
MIN_COVERAGE_NORMAL=20      # Mutect2 recommended minimum
MIN_COVERAGE_VARSCAN=8      # VarScan minimum usable depth

NORMAL_SAMPLES=("423_D0_old" "424_D0__old")
TUMOR_SAMPLES_LIST=("428_D20_new" "34_D52_old" "36_D99_new" "38_D99_new" "42_D122_old")

mkdir -p "${OUTDIR}"

# --- TSV header --------------------------------------------------------------
printf "sample\trole\tsize_gb\tindex_present\tsort_order\tref_build\t\
contig_naming\trg_present\trg_sm\trg_id\trg_lb\trg_pl\tdup_marked\t\
total_reads\tmapped_pct\tproperly_paired_pct\test_coverage_x\t\
mutect2_ready\tvarscan_ready\tnotes\n" > "${RESULTS_TSV}"

# --- Helper: stream first 10 MB of BAM from S3 -------------------------------
fetch_header() {
    local s3_path="$1"
    local size_bytes="$2"
    local tmp="$3"
    local end=$(( 10485760 < size_bytes ? 10485759 : size_bytes - 1 ))
    aws s3api get-object \
        --bucket bam-wes \
        --key "${s3_path#s3://bam-wes/}" \
        --range "bytes=0-${end}" \
        "${tmp}" > /dev/null 2>&1
}

# --- Helper: sample role -----------------------------------------------------
get_role() {
    local s="$1"
    for n in "${NORMAL_SAMPLES[@]}"; do [[ "$s" == "$n" ]] && echo "normal" && return; done
    for t in "${TUMOR_SAMPLES_LIST[@]}"; do [[ "$s" == "$t" ]] && echo "tumor" && return; done
    echo "unknown"
}

# --- Collect BAM list from S3 ------------------------------------------------
log "Listing BAMs in ${S3_BAM_PREFIX}/"
mapfile -t BAM_FILES < <(aws s3 ls "${S3_BAM_PREFIX}/" \
    | awk '{print $NF}' | grep '\.bam$' | sort)

[[ ${#BAM_FILES[@]} -eq 0 ]] && { log "[ERROR] No BAMs found"; exit 1; }
log "Found ${#BAM_FILES[@]} BAM(s)"

# =============================================================================
for BAM_FILE in "${BAM_FILES[@]}"; do
    S3_PATH="${S3_BAM_PREFIX}/${BAM_FILE}"
    SAMPLE="${BAM_FILE%.bam}"
    ROLE=$(get_role "${SAMPLE}")

    echo ""
    echo "========================================"
    log "Evaluating: ${SAMPLE}  [${ROLE}]"
    echo "========================================"

    NOTES=""
    MUTECT2_READY="PASS"
    VARSCAN_READY="PASS"

    # -------------------------------------------------------------------------
    # 1. Size + index
    # -------------------------------------------------------------------------
    SIZE_BYTES=$(aws s3 ls "${S3_PATH}" | awk '{print $3}')
    SIZE_GB=$(awk "BEGIN {printf \"%.2f\", ${SIZE_BYTES}/1073741824}")
    log "  Size: ${SIZE_GB} GB"

    INDEX_PRESENT="NO"
    if aws s3 ls "${S3_PATH}.bai" &>/dev/null \
       || aws s3 ls "${S3_BAM_PREFIX}/${SAMPLE}.bai" &>/dev/null; then
        INDEX_PRESENT="YES"
        log "  Index (.bai): present"
    else
        INDEX_PRESENT="NO"
        MUTECT2_READY="FAIL"
        VARSCAN_READY="FAIL"
        NOTES="${NOTES}NO_BAI_INDEX;"
        log "  [FAIL] No .bai index — required by both Mutect2 and VarScan"
    fi

    # -------------------------------------------------------------------------
    # 2. BAM header: sort order, reference, @RG tags, contigs, dup marking
    # -------------------------------------------------------------------------
    HEADER_TMP="/tmp/${SAMPLE}_hdr.bam"
    fetch_header "${S3_PATH}" "${SIZE_BYTES}" "${HEADER_TMP}"
    HEADER=$("${SAMTOOLS}" view -H "${HEADER_TMP}" 2>/dev/null)
    rm -f "${HEADER_TMP}"

    # Sort order
    SORT_ORDER=$(echo "${HEADER}" | awk '/^@HD/{
        for(i=1;i<=NF;i++) if($i~/^SO:/) {sub("SO:","",$i); print $i}
    }')
    SORT_ORDER="${SORT_ORDER:-unknown}"
    if [[ "${SORT_ORDER}" != "coordinate" ]]; then
        MUTECT2_READY="FAIL"; VARSCAN_READY="FAIL"
        NOTES="${NOTES}NOT_COORDINATE_SORTED;"
        log "  [FAIL] Sort order: ${SORT_ORDER} — must be coordinate"
    else
        log "  Sort order: coordinate [OK]"
    fi

    # Reference build
    REF_BUILD="unknown"
    if echo "${HEADER}" | grep -qiE "mm10|GRCm38|mus_musculus"; then
        REF_BUILD="mm10/GRCm38"
    elif echo "${HEADER}" | grep -qiE "hg38|GRCh38"; then REF_BUILD="hg38/GRCh38"
    elif echo "${HEADER}" | grep -qiE "hg19|GRCh37"; then REF_BUILD="hg19/GRCh37"
    fi
    log "  Reference build: ${REF_BUILD}"
    if [[ "${REF_BUILD}" != "mm10/GRCm38" ]]; then
        NOTES="${NOTES}REF_MISMATCH(${REF_BUILD});"
        log "  [WARN] Expected mm10/GRCm38 — pipeline reference is mm10"
        [[ "${MUTECT2_READY}" == "PASS" ]] && MUTECT2_READY="WARN"
        [[ "${VARSCAN_READY}" == "PASS" ]] && VARSCAN_READY="WARN"
    fi

    # Contig naming: chr1 vs 1
    FIRST_CONTIG=$(echo "${HEADER}" | awk '/^@SQ/{
        for(i=1;i<=NF;i++) if($i~/^SN:/) {sub("SN:","",$i); print $i; exit}
    }')
    if [[ "${FIRST_CONTIG}" == chr* ]]; then
        CONTIG_STYLE="UCSC(chr)"
        log "  Contig naming: UCSC (chr1, chr2...) [OK for mm10]"
    else
        CONTIG_STYLE="Ensembl(no-chr)"
        NOTES="${NOTES}CONTIG_NO_CHR_PREFIX;"
        log "  [WARN] Contig naming: Ensembl style (1, 2...) — reference uses chr prefix"
        [[ "${MUTECT2_READY}" == "PASS" ]] && MUTECT2_READY="WARN"
        [[ "${VARSCAN_READY}" == "PASS" ]] && VARSCAN_READY="WARN"
    fi

    # @RG completeness — Mutect2 requires SM; LB, PL, ID strongly recommended
    RG_COUNT=$(echo "${HEADER}" | grep -c "^@RG" || true)
    RG_SM=""; RG_ID=""; RG_LB=""; RG_PL=""

    if [[ "${RG_COUNT}" -gt 0 ]]; then
        RG_SM=$(echo "${HEADER}"  | grep "^@RG" | grep -oP "SM:\K\S+" | head -1 || true)
        RG_ID=$(echo "${HEADER}"  | grep "^@RG" | grep -oP "ID:\K\S+" | head -1 || true)
        RG_LB=$(echo "${HEADER}"  | grep "^@RG" | grep -oP "LB:\K\S+" | head -1 || true)
        RG_PL=$(echo "${HEADER}"  | grep "^@RG" | grep -oP "PL:\K\S+" | head -1 || true)
        log "  @RG present (${RG_COUNT}): SM=${RG_SM:-MISSING} ID=${RG_ID:-MISSING} LB=${RG_LB:-MISSING} PL=${RG_PL:-MISSING}"
        if [[ -z "${RG_SM}" ]]; then
            MUTECT2_READY="FAIL"
            NOTES="${NOTES}RG_SM_MISSING;"
            log "  [FAIL] @RG SM tag missing — required by Mutect2"
        fi
        if [[ -z "${RG_LB}" ]]; then
            NOTES="${NOTES}RG_LB_MISSING;"
            log "  [WARN] @RG LB (library) tag missing — recommended for MarkDuplicates"
        fi
        if [[ -z "${RG_PL}" ]]; then
            NOTES="${NOTES}RG_PL_MISSING;"
            log "  [WARN] @RG PL (platform) tag missing"
        fi
    else
        RG_SM="NO"; RG_ID="NO"; RG_LB="NO"; RG_PL="NO"
        MUTECT2_READY="FAIL"; VARSCAN_READY="WARN"
        NOTES="${NOTES}NO_READ_GROUPS;"
        log "  [FAIL] No @RG read group — required by Mutect2"
    fi

    # Duplicate marking status
    DUP_MARKED="NO"
    if echo "${HEADER}" | grep -q "ID:MarkDuplicates"; then
        DUP_MARKED="YES"
        log "  Duplicate marking: YES (MarkDuplicates in header)"
    else
        DUP_MARKED="NO"
        NOTES="${NOTES}DUPS_NOT_MARKED(raw_bam);"
        log "  Duplicate marking: NO — raw BAM, requires preprocessing before variant calling"
        [[ "${MUTECT2_READY}" == "PASS" ]] && MUTECT2_READY="WARN"
    fi

    # -------------------------------------------------------------------------
    # 3. Flagstat — reuse existing file or stream from S3
    # -------------------------------------------------------------------------
    FLAGSTAT_FILE="${FLAGSTAT_DIR}/${SAMPLE}.flagstat.txt"
    TOTAL_READS="NA"; MAPPED_PCT="NA"; PROPERLY_PAIRED_PCT="NA"; EST_COVERAGE="NA"

    if [[ -f "${FLAGSTAT_FILE}" ]]; then
        log "  Flagstat: using existing ${FLAGSTAT_FILE}"
    else
        log "  Flagstat: streaming from S3 (this may take a while)..."
        mkdir -p "${FLAGSTAT_DIR}"
        aws s3 cp "${S3_PATH}" - 2>/dev/null \
            | "${SAMTOOLS}" flagstat - > "${FLAGSTAT_FILE}" 2>/dev/null \
            || { log "  [WARN] flagstat failed"; NOTES="${NOTES}FLAGSTAT_FAILED;"; }
    fi

    if [[ -f "${FLAGSTAT_FILE}" ]]; then
        TOTAL_READS=$(awk '/in total/{print $1; exit}'           "${FLAGSTAT_FILE}")
        MAPPED_READS=$(awk '/primary mapped/{print $1; exit}'    "${FLAGSTAT_FILE}")
        PROPERLY_PAIRED=$(awk '/properly paired/{print $1; exit}' "${FLAGSTAT_FILE}")

        if [[ -n "${TOTAL_READS}" && "${TOTAL_READS}" -gt 0 ]]; then
            MAPPED_PCT=$(awk "BEGIN {printf \"%.1f\", 100*${MAPPED_READS}/${TOTAL_READS}}")
            PROPERLY_PAIRED_PCT=$(awk "BEGIN {printf \"%.1f\", 100*${PROPERLY_PAIRED}/${TOTAL_READS}}")
            EST_COVERAGE=$(awk "BEGIN {printf \"%.0f\", (${MAPPED_READS} * ${READ_LENGTH}) / ${EXOME_SIZE_BP}}")
            log "  Total reads       : ${TOTAL_READS}"
            log "  Mapped            : ${MAPPED_PCT}%"
            log "  Properly paired   : ${PROPERLY_PAIRED_PCT}%"
            log "  Est. coverage     : ~${EST_COVERAGE}x  (WES, ${EXOME_SIZE_BP}bp target, ${READ_LENGTH}bp reads)"

            # Properly paired threshold
            PP=$(echo "${PROPERLY_PAIRED_PCT}" | cut -d. -f1)
            if [[ "${PP}" -lt 90 ]]; then
                NOTES="${NOTES}LOW_PROPER_PAIRS(${PROPERLY_PAIRED_PCT}%);"
                log "  [WARN] Properly paired ${PROPERLY_PAIRED_PCT}% < 90%"
                [[ "${MUTECT2_READY}" == "PASS" ]] && MUTECT2_READY="WARN"
                [[ "${VARSCAN_READY}" == "PASS" ]] && VARSCAN_READY="WARN"
            fi

            # Coverage thresholds
            MIN_COV=$( [[ "${ROLE}" == "tumor" ]] \
                && echo "${MIN_COVERAGE_TUMOR}" \
                || echo "${MIN_COVERAGE_NORMAL}" )
            if [[ "${EST_COVERAGE}" -lt "${MIN_COVERAGE_VARSCAN}" ]]; then
                VARSCAN_READY="FAIL"
                NOTES="${NOTES}LOW_COVERAGE_VARSCAN(${EST_COVERAGE}x<${MIN_COVERAGE_VARSCAN}x);"
                log "  [FAIL] Coverage ~${EST_COVERAGE}x below VarScan minimum (${MIN_COVERAGE_VARSCAN}x)"
            elif [[ "${EST_COVERAGE}" -lt "${MIN_COV}" ]]; then
                [[ "${MUTECT2_READY}" == "PASS" ]] && MUTECT2_READY="WARN"
                NOTES="${NOTES}LOW_COVERAGE_MUTECT2(${EST_COVERAGE}x<${MIN_COV}x);"
                log "  [WARN] Coverage ~${EST_COVERAGE}x below recommended Mutect2 ${MIN_COV}x for ${ROLE}"
            else
                log "  Coverage check    : OK (${EST_COVERAGE}x >= ${MIN_COV}x)"
            fi
        fi
    fi

    # -------------------------------------------------------------------------
    # 4. Summary
    # -------------------------------------------------------------------------
    echo "  ----------------------------------------"
    log "  Mutect2 readiness : ${MUTECT2_READY}"
    log "  VarScan readiness : ${VARSCAN_READY}"
    echo "  ----------------------------------------"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SAMPLE}" "${ROLE}" "${SIZE_GB}" "${INDEX_PRESENT}" \
        "${SORT_ORDER}" "${REF_BUILD}" "${CONTIG_STYLE}" \
        "${RG_COUNT}" "${RG_SM:-NA}" "${RG_ID:-NA}" "${RG_LB:-NA}" "${RG_PL:-NA}" \
        "${DUP_MARKED}" \
        "${TOTAL_READS}" "${MAPPED_PCT}" "${PROPERLY_PAIRED_PCT}" "${EST_COVERAGE}" \
        "${MUTECT2_READY}" "${VARSCAN_READY}" "${NOTES:-none}" \
        >> "${RESULTS_TSV}"
done

# =============================================================================
# Summary table
# =============================================================================
echo ""
echo "========================================"
log "EVALUATION SUMMARY"
echo "========================================"
printf "%-18s %-8s %-6s %-10s %-12s %-10s\n" \
    "SAMPLE" "ROLE" "COV_X" "MUTECT2" "VARSCAN" "NOTES"
echo "----------------------------------------------------------------------"
tail -n +2 "${RESULTS_TSV}" | while IFS=$'\t' read -r \
    sample role size_gb index sort ref contig rg_n sm id lb pl dup \
    tot_reads mapped_pct pp_pct cov mutect2 varscan notes; do
    printf "%-18s %-8s %-6s %-10s %-12s %-10s\n" \
        "${sample}" "${role}" "${cov}x" "${mutect2}" "${varscan}" \
        "$(echo "${notes}" | cut -c1-30)"
done

echo ""
log "Results TSV: ${RESULTS_TSV}"
log "Done."
