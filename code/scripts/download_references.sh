#!/bin/bash
# =============================================================================
# download_references.sh
# Download the two reference files needed by neoantigen_pipeline.sh
# that are not supplied in data/ by the client:
#
#   reference_db/merged_exons.bed   — exon intervals for mm10 (UCSC chr-prefix)
#   reference_db/proteins.fa        — Ensembl GRCm38 protein FASTA
#
# Run once before starting the pipeline:
#   bash scripts/download_references.sh
# =============================================================================

set -euo pipefail
WORK_DIR="$(cd "$(dirname "$0")/.." && pwd)"
REF_DIR="${WORK_DIR}/reference_db"
mkdir -p "${REF_DIR}"

# ── 1. Exon BED (mm10, UCSC chr-prefix) ──────────────────────────────────────
# Source: GENCODE vM25 comprehensive annotation (GRCm38 = mm10)
#         Filtered to exon features, sorted and merged
EXON_BED="${REF_DIR}/merged_exons.bed"

if [[ ! -f "$EXON_BED" ]]; then
    echo "[1/2] Downloading GENCODE vM25 GTF for mm10..."
    GTF_GZ="${REF_DIR}/gencode.vM25.annotation.gtf.gz"
    if [[ ! -f "$GTF_GZ" ]]; then
        wget -q --show-progress -O "$GTF_GZ" \
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"
    fi

    echo "  Extracting exon intervals and merging..."
    zcat "$GTF_GZ" \
        | awk '$3=="exon" {
            split($0,a,"\t");
            # GENCODE uses chr-prefix already
            print a[1]"\t"(a[4]-1)"\t"a[5]
          }' \
        | sort -k1,1 -k2,2n \
        | bedtools merge -i - \
        > "$EXON_BED"

    echo "  Done → ${EXON_BED}"
    echo "  Lines: $(wc -l < "$EXON_BED")"
    # Optionally remove GTF to save space
    # rm "$GTF_GZ"
else
    echo "[1/2] merged_exons.bed already exists — skipping."
fi

# ── 2. Ensembl GRCm38 protein FASTA ──────────────────────────────────────────
# Source: Ensembl release 102 (GRCm38.p6) — same build as GRCm38.86 SnpEff DB
PROT_FA="${REF_DIR}/proteins.fa"
PROT_GZ="${REF_DIR}/Mus_musculus.GRCm38.pep.all.fa.gz"

if [[ ! -f "$PROT_FA" ]]; then
    echo "[2/2] Downloading Ensembl GRCm38 protein FASTA..."
    if [[ ! -f "$PROT_GZ" ]]; then
        wget -q --show-progress -O "$PROT_GZ" \
            "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz"
    fi

    echo "  Decompressing..."
    gunzip -c "$PROT_GZ" > "$PROT_FA"
    echo "  Done → ${PROT_FA}"
    echo "  Proteins: $(grep -c '^>' "$PROT_FA")"
    # Optionally remove compressed copy
    # rm "$PROT_GZ"
else
    echo "[2/2] proteins.fa already exists — skipping."
fi

echo ""
echo "Reference files ready:"
echo "  ${EXON_BED}"
echo "  ${PROT_FA}"
