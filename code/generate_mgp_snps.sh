#!/usr/bin/env bash
# =============================================================================
# generate_mgp_snps.sh — Download and prepare MGP v5 SNPs for GATK BQSR
#
# Source: EBI Mouse Genomes Project REL-1505 (GRCm38, dbSNP 142)
# Strategy: Stream VCF → filter to exon regions → convert to UCSC chr names
#           → bgzip + tabix index. Output ~400 MB vs 21 GB full file.
# Output: /home/ec2-user/ref/mm10/mgp_snps.vcf.gz (UCSC chr-prefixed)
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

MGP_URL="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
OUTVCF="${REF_DIR}/mgp_snps.vcf.gz"
EXON_ENSEMBL="/tmp/exons_ensembl.bed"

log "generate_mgp_snps.sh started: $(date)"
log "Output: ${OUTVCF}"

# --- 1. Create Ensembl-style exon BED (strip chr prefix; chrM → MT) ----------
log "Creating Ensembl-style exon BED from ${EXON_BED}..."
awk 'BEGIN{OFS="\t"} {
    chrom = $1
    sub(/^chr/, "", chrom)
    if (chrom == "M") chrom = "MT"
    print chrom, $2, $3
}' "${EXON_BED}" | sort -k1,1 -k2,2n > "${EXON_ENSEMBL}"
log "  Exon intervals: $(wc -l < ${EXON_ENSEMBL})"

# --- 2. Stream VCF, filter to exons (-T), convert to UCSC, bgzip ------------
log "Streaming MGP VCF from EBI and filtering to exon targets..."
log "  URL: ${MGP_URL}"

wget -qO- "${MGP_URL}" \
    | ${BCFTOOLS} view -T "${EXON_ENSEMBL}" \
    | awk 'BEGIN{OFS="\t"}
        /^##contig=<ID=/ {
            # rename contigs in header: 1→chr1, MT→chrM
            sub(/^##contig=<ID=MT/, "##contig=<ID=chrM")
            sub(/^##contig=<ID=([0-9XY])/, "##contig=<ID=chr\\1")
            print; next
        }
        /^#/ { print; next }
        {
            # rename CHROM field
            if ($1 == "MT") $1 = "chrM"
            else if ($1 ~ /^[0-9XY]/) $1 = "chr" $1
            print
        }' \
    | ${CONDA_BIN}/bgzip -c > "${OUTVCF}"

# --- 3. Tabix index ----------------------------------------------------------
log "Indexing with tabix..."
${CONDA_BIN}/tabix -p vcf "${OUTVCF}"

# --- 4. Report ---------------------------------------------------------------
NSNP=$(${BCFTOOLS} stats "${OUTVCF}" | grep "^SN.*number of SNPs" | awk '{print $NF}')
SIZE=$(du -sh "${OUTVCF}" | cut -f1)
log "Done."
log "  File:  ${OUTVCF}  (${SIZE})"
log "  SNPs:  ${NSNP}"
log "generate_mgp_snps.sh complete: $(date)"
