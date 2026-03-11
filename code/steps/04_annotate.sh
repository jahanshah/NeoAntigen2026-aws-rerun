#!/usr/bin/env bash
# =============================================================================
# Step 4 — SnpEff Annotation (GRCm38.86)
# Annotates PASS+exon VCFs with functional consequences.
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

S3_IN="${S3_RESULTS}/vcf_filtered"
INDIR="${RESULTS_DIR}/vcf_filtered"
OUTDIR="${RESULTS_DIR}/annotated"
S3_OUT="${S3_RESULTS}/annotated"
mkdir -p "${OUTDIR}" "${INDIR}"

for TUMOR in "${TUMOR_SAMPLES[@]}"; do
    ANN_VCF="${OUTDIR}/${TUMOR}.annotated.vcf.gz"
    ANN_STATS="${OUTDIR}/${TUMOR}.snpeff_stats.html"

    if aws s3 ls "${S3_OUT}/${TUMOR}.annotated.vcf.gz" &>/dev/null; then
        skip "${TUMOR} annotated VCF (S3)"
        continue
    fi

    log "========================================"
    log "Step 4 — Annotating: ${TUMOR}"
    log "========================================"

    EXON_VCF="${INDIR}/${TUMOR}.filtered.exon.vcf.gz"
    [[ ! -f "${EXON_VCF}" ]] && aws s3 cp "${S3_IN}/${TUMOR}.filtered.exon.vcf.gz" "${EXON_VCF}"
    [[ ! -f "${EXON_VCF}.tbi" ]] && ${BCFTOOLS} index -t "${EXON_VCF}"

    log "  Running SnpEff (${SNPEFF_DB})..."
    ${SNPEFF} \
        -Xmx8g \
        -v "${SNPEFF_DB}" \
        -stats "${ANN_STATS}" \
        -cancer \
        -noLog \
        "${EXON_VCF}" \
        | ${BCFTOOLS} view -O z -o "${ANN_VCF}"
    ${BCFTOOLS} index -t "${ANN_VCF}"

    # Summary of effect consequences
    CODING=$(${BCFTOOLS} view "${ANN_VCF}" | grep -o "ANN=[^;]*" \
        | grep -oP "(?<=\|)[A-Z_]+(?=\|)" | sort | uniq -c | sort -rn | head -10 || true)
    log "  Top consequences:"
    echo "${CODING}" | while read -r line; do log "    ${line}"; done

    s3up "${ANN_VCF}"     "annotated/$(basename ${ANN_VCF})"
    s3up "${ANN_VCF}.tbi" "annotated/$(basename ${ANN_VCF}).tbi"
    s3up "${ANN_STATS}"   "annotated/$(basename ${ANN_STATS})"

    rm -f "${EXON_VCF}" "${EXON_VCF}.tbi"
    log "  Done: ${TUMOR}"
done

log "Step 4 complete."
