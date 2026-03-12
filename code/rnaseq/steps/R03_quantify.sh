#!/usr/bin/env bash
# =============================================================================
# R03_quantify.sh — Run featureCounts on all aligned BAMs
#
# Steps:
#   1. Download all BAMs from S3 (if not already local)
#   2. Run featureCounts on all BAMs together (multi-sample count matrix)
#   3. Fallback: if featureCounts fails, merge STAR GeneCounts tabs instead
#   4. Upload counts output to S3
# =============================================================================

set -euo pipefail

source /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh

COUNTS_DIR="${RNASEQ_RESULTS}/counts"
BAM_DIR="${RNASEQ_TMP}/bams"
COUNTS_OUT="${COUNTS_DIR}/raw_counts.txt"
S3_COUNTS_OUT="${S3_RNASEQ_OUT}/counts/raw_counts.txt"

mkdir -p "${COUNTS_DIR}" "${BAM_DIR}"

# =============================================================================
# Skip if counts already on S3
# =============================================================================
if aws s3 ls "${S3_COUNTS_OUT}" &>/dev/null; then
    log "[SKIP] raw_counts.txt already on S3: ${S3_COUNTS_OUT}"
    # Still download for R04 to use
    if [[ ! -f "${COUNTS_OUT}" ]]; then
        log "Downloading existing counts from S3..."
        aws s3 cp "${S3_COUNTS_OUT}" "${COUNTS_OUT}"
    fi
    exit 0
fi

# =============================================================================
# 1. Download BAMs from S3
# =============================================================================
log "--- Downloading BAMs from S3 ---"
ALL_BAMS=()

for SAMPLE in "${RNASEQ_SAMPLES[@]}"; do
    S3_BAM="${S3_RNASEQ_OUT}/bam/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    LOCAL_BAM="${BAM_DIR}/${SAMPLE}.Aligned.sortedByCoord.out.bam"
    LOCAL_BAI="${LOCAL_BAM}.bai"

    if [[ -f "${LOCAL_BAM}" && -f "${LOCAL_BAI}" ]]; then
        log "[OK] ${SAMPLE}: BAM already local"
        ALL_BAMS+=("${LOCAL_BAM}")
        continue
    fi

    if aws s3 ls "${S3_BAM}" &>/dev/null; then
        log "${SAMPLE}: downloading BAM from S3..."
        aws s3 cp "${S3_BAM}"       "${LOCAL_BAM}"
        aws s3 cp "${S3_BAM}.bai"   "${LOCAL_BAI}" || \
            samtools index -@ "${RNASEQ_THREADS}" "${LOCAL_BAM}"
        ALL_BAMS+=("${LOCAL_BAM}")
    else
        log "[WARN] ${SAMPLE}: BAM not found on S3 (${S3_BAM}) — skipping in quantification"
    fi
done

if [[ ${#ALL_BAMS[@]} -eq 0 ]]; then
    log "[ERROR] No BAMs available for quantification. Run R02_align.sh first."
    exit 1
fi

log "Quantifying ${#ALL_BAMS[@]} BAMs."

# =============================================================================
# 2. Run featureCounts
# =============================================================================
log "--- Running featureCounts ---"
T0=$SECONDS

featureCounts_success=0

featureCounts \
    -T "${RNASEQ_THREADS}" \
    -p --countReadPairs \
    -s 2 \
    -a "${GTF}" \
    -o "${COUNTS_OUT}" \
    "${ALL_BAMS[@]}" \
    && featureCounts_success=1 \
    || log "[WARN] featureCounts failed — will fall back to STAR GeneCounts"

if [[ ${featureCounts_success} -eq 1 ]]; then
    log "featureCounts complete in $(( SECONDS - T0 ))s."
    log "Output: ${COUNTS_OUT}"
else
    # =========================================================================
    # 3. Fallback: merge STAR ReadsPerGene.out.tab files
    # =========================================================================
    log "--- Fallback: merging STAR GeneCounts tabs ---"
    STAR_COUNTS_DIR="${RNASEQ_TMP}/star_counts"
    mkdir -p "${STAR_COUNTS_DIR}"

    # Download ReadsPerGene tabs from S3
    STAR_TAB_FILES=()
    for SAMPLE in "${RNASEQ_SAMPLES[@]}"; do
        S3_TAB="${S3_RNASEQ_OUT}/counts/${SAMPLE}.ReadsPerGene.out.tab"
        LOCAL_TAB="${STAR_COUNTS_DIR}/${SAMPLE}.ReadsPerGene.out.tab"
        if aws s3 ls "${S3_TAB}" &>/dev/null; then
            aws s3 cp "${S3_TAB}" "${LOCAL_TAB}" --quiet
            STAR_TAB_FILES+=("${LOCAL_TAB}")
        else
            log "[WARN] ${SAMPLE}: STAR GeneCounts tab not on S3 — skipping"
        fi
    done

    if [[ ${#STAR_TAB_FILES[@]} -eq 0 ]]; then
        log "[ERROR] No STAR GeneCounts tabs available either. Cannot quantify."
        exit 1
    fi

    # Merge: column 2 = unstranded, column 3 = forward, column 4 = reverse
    # Use column 4 (reverse-stranded, matching featureCounts -s 2)
    # Header: Geneid + sample names
    log "Merging ${#STAR_TAB_FILES[@]} STAR GeneCounts tabs (column 4, reverse-stranded)..."

    # Build merged counts file using Python (more portable than awk for multiple files)
    python3 - "${STAR_TAB_FILES[@]}" <<'PYEOF'
import sys, os

files = sys.argv[1:]
counts = {}  # gene_id -> [counts per sample]
sample_names = []

for f in files:
    sample = os.path.basename(f).replace(".ReadsPerGene.out.tab", "")
    sample_names.append(sample)
    with open(f) as fh:
        for i, line in enumerate(fh):
            if i < 4:
                continue  # skip N_unmapped etc. summary lines
            parts = line.rstrip("\n").split("\t")
            gene_id = parts[0]
            count = parts[3]  # column 4 = reverse-stranded
            if gene_id not in counts:
                counts[gene_id] = []
            counts[gene_id].append(count)

# Write output in featureCounts-like format
out_file = sys.argv[0]  # placeholder — write to stdout and redirect
import sys
# Write to stdout; bash will redirect
header = "Geneid\t" + "\t".join(sample_names)
print(header)
for gene_id, cnt_list in counts.items():
    # Pad missing values with 0 for samples that joined late
    while len(cnt_list) < len(sample_names):
        cnt_list.append("0")
    print(gene_id + "\t" + "\t".join(cnt_list))
PYEOF
    # Note: The python block above prints to stdout; redirect to file:
    python3 - "${STAR_TAB_FILES[@]}" > "${COUNTS_OUT}" <<'PYEOF'
import sys, os

files = sys.argv[1:]
counts = {}
sample_names = []
gene_order = []

for f in files:
    sample = os.path.basename(f).replace(".ReadsPerGene.out.tab", "")
    sample_names.append(sample)
    with open(f) as fh:
        for i, line in enumerate(fh):
            if i < 4:
                continue
            parts = line.rstrip("\n").split("\t")
            gene_id = parts[0]
            count = int(parts[3])
            if gene_id not in counts:
                counts[gene_id] = {}
                gene_order.append(gene_id)
            counts[gene_id][sample] = count

print("Geneid\t" + "\t".join(sample_names))
for gene_id in gene_order:
    row = [str(counts[gene_id].get(s, 0)) for s in sample_names]
    print(gene_id + "\t" + "\t".join(row))
PYEOF

    log "STAR GeneCounts merge complete: ${COUNTS_OUT}"
fi

# =============================================================================
# 4. Upload counts to S3
# =============================================================================
log "--- Uploading counts to S3 ---"
aws s3 cp "${COUNTS_OUT}" "${S3_COUNTS_OUT}"
log "Uploaded: ${S3_COUNTS_OUT}"

# Also upload the featureCounts summary if present
if [[ -f "${COUNTS_OUT}.summary" ]]; then
    aws s3 cp "${COUNTS_OUT}.summary" "${S3_RNASEQ_OUT}/counts/raw_counts.txt.summary"
fi

# Clean up downloaded BAMs to free disk space
log "--- Cleaning up local BAMs ---"
rm -rf "${BAM_DIR}"
log "R03_quantify.sh complete."
