#!/usr/bin/env bash
# =============================================================================
# R01_setup_ref.sh — Install tools, download GTF, build/restore STAR index
#
# All reference files are cached in S3 under:
#   s3://neoantigen2026-rerun/data/reference/rnaseq/
#     mm10_chr.gtf.gz      — UCSC chr-prefixed GTF (GRCm38.86)
#     star_index.tar.gz    — pre-built STAR 2.7 index
#
# Restore order (fast path first):
#   1. Check S3 for GTF → download; else fetch from Ensembl + convert
#   2. Check S3 for STAR index → extract; else build + upload to S3
# =============================================================================

set -euo pipefail

source /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh

MAMBA="/home/ec2-user/miniforge3/bin/mamba"

# S3 reference cache
S3_REF_BASE="${S3_ROOT}/data/reference/rnaseq"
S3_GTF="${S3_REF_BASE}/mm10_chr.gtf.gz"
S3_STAR_INDEX="${S3_REF_BASE}/star_index.tar.gz"

# Local raw (Ensembl) GTF — only needed if rebuilding chr GTF
GTF_RAW="/home/ec2-user/ref/mm10/mm10.gtf"
GTF_URL="https://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.gtf.gz"
GTF_GZ="${GTF_RAW}.gz"

# =============================================================================
# 1. Install required tools via mamba
# =============================================================================
log "--- Checking/installing tools ---"

install_if_missing() {
    local tool="$1"
    local pkg="${2:-$1}"
    if command -v "${tool}" &>/dev/null; then
        log "[OK] ${tool} already installed: $(command -v "${tool}")"
    else
        log "[INSTALL] Installing ${pkg} via mamba..."
        "${MAMBA}" install -y -c conda-forge -c bioconda "${pkg}"
        log "[OK] ${pkg} installed."
    fi
}

install_if_missing "STAR"          "star"
install_if_missing "featureCounts" "subread"
install_if_missing "fastp"         "fastp"

# =============================================================================
# 2. GTF — restore from S3 or build from Ensembl
# =============================================================================
log "--- Checking chr-prefixed GTF (GRCm38.86) ---"
mkdir -p "$(dirname "${GTF}")"

if [[ -f "${GTF}" ]]; then
    log "[SKIP] GTF already present locally: ${GTF}"

elif aws s3 ls "${S3_GTF}" &>/dev/null; then
    log "Restoring GTF from S3: ${S3_GTF}"
    aws s3 cp "${S3_GTF}" "${GTF}.gz"
    gunzip -f "${GTF}.gz"
    log "[OK] GTF restored: ${GTF}"

else
    # Download raw Ensembl GTF and convert chr names
    if [[ ! -f "${GTF_RAW}" ]]; then
        log "Downloading GTF from Ensembl release 86..."
        wget -q --show-progress -O "${GTF_GZ}" "${GTF_URL}"
        gunzip -f "${GTF_GZ}"
        log "[OK] Raw GTF: ${GTF_RAW}"
    fi

    log "Converting Ensembl chr names → UCSC (1→chr1, MT→chrM)..."
    awk '
/^#/ { print; next }
{
    chr = $1
    if (chr == "MT") chr = "chrM"
    else if (chr ~ /^[0-9]+$/ || chr == "X" || chr == "Y") chr = "chr" chr
    else next
    $1 = chr; print
}
OFS="\t"
' "${GTF_RAW}" > "${GTF}"
    log "[OK] chr-prefixed GTF built: ${GTF}"

    # Upload to S3 for future runs
    log "Uploading GTF to S3: ${S3_GTF}"
    gzip -c "${GTF}" | aws s3 cp - "${S3_GTF}"
    log "[OK] GTF uploaded to S3."
fi

# =============================================================================
# 3. STAR index — restore from S3 or build from scratch
# =============================================================================
log "--- Checking STAR index ---"
mkdir -p "${STAR_INDEX}"

if [[ -f "${STAR_INDEX}/SAindex" ]]; then
    log "[SKIP] STAR index already present locally: ${STAR_INDEX}"

elif aws s3 ls "${S3_STAR_INDEX}" &>/dev/null; then
    log "Restoring STAR index from S3: ${S3_STAR_INDEX}"
    local_tar="${RNASEQ_TMP}/star_index.tar.gz"
    mkdir -p "${RNASEQ_TMP}"
    aws s3 cp "${S3_STAR_INDEX}" "${local_tar}"
    log "Extracting to $(dirname "${STAR_INDEX}")..."
    tar -xzf "${local_tar}" -C "$(dirname "${STAR_INDEX}")"
    rm -f "${local_tar}"
    log "[OK] STAR index restored from S3."

else
    log "Building STAR index from scratch (this takes ~30 min)..."
    log "  Genome : ${REF}"
    log "  GTF    : ${GTF}"
    log "  Threads: ${RNASEQ_THREADS}"

    T0=$SECONDS
    STAR \
        --runMode genomeGenerate \
        --genomeDir "${STAR_INDEX}" \
        --genomeFastaFiles "${REF}" \
        --sjdbGTFfile "${GTF}" \
        --sjdbOverhang 149 \
        --genomeSAindexNbases 14 \
        --limitGenomeGenerateRAM 26000000000 \
        --runThreadN "${RNASEQ_THREADS}"
    log "STAR index built in $(( SECONDS - T0 ))s."

    # Upload to S3 — all future runs restore from here
    log "Uploading STAR index to S3: ${S3_STAR_INDEX}"
    T0=$SECONDS
    tar -czf - -C "$(dirname "${STAR_INDEX}")" "$(basename "${STAR_INDEX}")" \
        | aws s3 cp - "${S3_STAR_INDEX}"
    log "STAR index uploaded in $(( SECONDS - T0 ))s."
fi

log "R01_setup_ref.sh complete."
