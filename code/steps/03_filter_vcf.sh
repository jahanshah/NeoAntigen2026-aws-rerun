#!/usr/bin/env bash
# =============================================================================
# Step 3 — FilterMutectCalls + SelectVariants (exons only)
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

S3_IN="${S3_RESULTS}/mutect2"
INDIR="${RESULTS_DIR}/mutect2"
OUTDIR="${RESULTS_DIR}/vcf_filtered"
S3_OUT="${S3_RESULTS}/vcf_filtered"
mkdir -p "${OUTDIR}" "${INDIR}/vcf_raw" "${INDIR}/stats"

for TUMOR in "${TUMOR_SAMPLES[@]}"; do
    FILT_VCF="${OUTDIR}/${TUMOR}.filtered.vcf.gz"
    EXON_VCF="${OUTDIR}/${TUMOR}.filtered.exon.vcf.gz"

    if aws s3 ls "${S3_OUT}/${TUMOR}.filtered.exon.vcf.gz" &>/dev/null; then
        skip "${TUMOR} filtered exon VCF (S3)"
        continue
    fi

    log "========================================"
    log "Step 3 — Filtering: ${TUMOR}"
    log "========================================"

    # Ensure raw VCF and stats are local
    RAW_VCF="${INDIR}/vcf_raw/${TUMOR}.vcf.gz"
    STATS="${INDIR}/stats/${TUMOR}.stats"
    [[ ! -f "${RAW_VCF}" ]] && aws s3 cp "${S3_IN}/vcf_raw/${TUMOR}.vcf.gz"     "${RAW_VCF}"
    [[ ! -f "${RAW_VCF}.tbi" ]] && aws s3 cp "${S3_IN}/vcf_raw/${TUMOR}.vcf.gz.tbi" "${RAW_VCF}.tbi" 2>/dev/null || \
        ${BCFTOOLS} index -t "${RAW_VCF}"
    [[ ! -f "${STATS}" ]] && aws s3 cp "${S3_IN}/stats/${TUMOR}.stats" "${STATS}" 2>/dev/null || true

    # --- 3a. FilterMutectCalls ----------------------------------------------
    log "  FilterMutectCalls..."
    FILTER_ARGS="-R ${REF} -V ${RAW_VCF} -O ${FILT_VCF}"
    [[ -f "${STATS}" ]] && FILTER_ARGS="${FILTER_ARGS} --stats ${STATS}"
    ${GATK} --java-options "${JAVA_OPTS}" FilterMutectCalls \
        ${FILTER_ARGS} 2>&1 | tail -3

    # --- 3b. Keep PASS only -------------------------------------------------
    log "  Selecting PASS variants..."
    PASS_VCF="${OUTDIR}/${TUMOR}.pass.vcf.gz"
    ${BCFTOOLS} view -f PASS -O z -o "${PASS_VCF}" "${FILT_VCF}"
    ${BCFTOOLS} index -t "${PASS_VCF}"

    # --- 3c. SelectVariants — restrict to exon regions ----------------------
    if [[ -f "${EXON_BED}" ]]; then
        log "  Restricting to exons..."
        ${GATK} --java-options "${JAVA_OPTS}" SelectVariants \
            -R "${REF}" \
            -V "${PASS_VCF}" \
            -L "${EXON_BED}" \
            --interval-padding 50 \
            -O "${EXON_VCF}" 2>&1 | tail -3
        rm -f "${PASS_VCF}" "${PASS_VCF}.tbi"
    else
        log "  [WARN] No exon BED — using all PASS variants"
        mv "${PASS_VCF}"     "${EXON_VCF}"
        mv "${PASS_VCF}.tbi" "${EXON_VCF}.tbi"
    fi

    # --- Stats --------------------------------------------------------------
    TOTAL=$(${BCFTOOLS} stats "${FILT_VCF}" | awk '/^SN.*number of records/{print $NF}')
    PASS=$(${BCFTOOLS} stats  "${EXON_VCF}" | awk '/^SN.*number of records/{print $NF}')
    SNV=$(${BCFTOOLS} stats   "${EXON_VCF}" | awk '/^SN.*number of SNPs/{print $NF}')
    IND=$(${BCFTOOLS} stats   "${EXON_VCF}" | awk '/^SN.*number of indels/{print $NF}')
    log "  Total: ${TOTAL}  PASS+exon: ${PASS}  SNV: ${SNV}  Indel: ${IND}"

    # --- Upload, clean local ------------------------------------------------
    for F in "${FILT_VCF}" "${FILT_VCF}.tbi" "${EXON_VCF}" "${EXON_VCF}.tbi"; do
        s3up "${F}" "vcf_filtered/$(basename ${F})"
    done
    rm -f "${FILT_VCF}" "${FILT_VCF}.tbi" "${RAW_VCF}" "${RAW_VCF}.tbi"
    log "  Done: ${TUMOR}"
done

log "Step 3 complete."
