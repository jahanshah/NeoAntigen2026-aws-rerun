#!/usr/bin/env bash
# =============================================================================
# R01_setup_ref.sh — Install tools, download GTF, build/restore STAR index
#
# Steps:
#   1. Install STAR, subread (featureCounts), fastp via mamba if not already present
#   2. Download GRCm38.86 GTF from Ensembl if not present locally
#   3. Check S3 for a pre-built STAR index tarball; download & extract if found
#   4. Otherwise build the STAR index from scratch and upload it to S3
# =============================================================================

set -euo pipefail

source /home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/config_rnaseq.sh

MAMBA="/home/ec2-user/miniforge3/bin/mamba"
S3_STAR_INDEX="${S3_ROOT}/data/reference/rnaseq/star_index.tar.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.gtf.gz"
GTF_RAW="/home/ec2-user/ref/mm10/mm10.gtf"      # Ensembl-style (1, X, MT)
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
# 2. Download GRCm38.86 GTF if not present
# =============================================================================
log "--- Checking GTF (GRCm38.86) ---"
mkdir -p "$(dirname "${GTF}")"

if [[ -f "${GTF_RAW}" ]]; then
    log "[SKIP] GTF (raw) already exists: ${GTF_RAW}"
else
    log "Downloading GTF from Ensembl release 86..."
    wget -q --show-progress -O "${GTF_GZ}" "${GTF_URL}"
    log "Decompressing GTF..."
    gunzip -f "${GTF_GZ}"
    log "[OK] GTF (raw) ready: ${GTF_RAW}"
fi

# Convert Ensembl chr names (1,X,MT) → UCSC (chr1,chrX,chrM) to match mm10.fa
if [[ -f "${GTF}" ]]; then
    log "[SKIP] chr-prefixed GTF already exists: ${GTF}"
else
    log "Building chr-prefixed GTF (Ensembl→UCSC chr names)..."
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
    log "[OK] chr-prefixed GTF: ${GTF}"
fi

# =============================================================================
# 3. Check S3 for pre-built STAR index
# =============================================================================
log "--- Checking for pre-built STAR index on S3 ---"
mkdir -p "${STAR_INDEX}"

if aws s3 ls "${S3_STAR_INDEX}" &>/dev/null; then
    if [[ -f "${STAR_INDEX}/SAindex" ]]; then
        log "[SKIP] STAR index already present locally: ${STAR_INDEX}"
    else
        log "Found STAR index on S3 — downloading and extracting..."
        local_tar="${RNASEQ_TMP}/star_index.tar.gz"
        mkdir -p "${RNASEQ_TMP}"
        aws s3 cp "${S3_STAR_INDEX}" "${local_tar}"
        log "Extracting to $(dirname "${STAR_INDEX}")..."
        tar -xzf "${local_tar}" -C "$(dirname "${STAR_INDEX}")"
        rm -f "${local_tar}"
        log "[OK] STAR index restored from S3."
    fi
else
    # ==========================================================================
    # 4. Build STAR index from scratch
    # ==========================================================================
    if [[ -f "${STAR_INDEX}/SAindex" ]]; then
        log "[SKIP] STAR index already built locally: ${STAR_INDEX}"
    else
        log "No pre-built index found — building STAR index (this takes ~30 min)..."
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
    fi

    # Upload the built index to S3 for future runs
    log "Uploading STAR index to S3: ${S3_STAR_INDEX}"
    T0=$SECONDS
    tar -czf - -C "$(dirname "${STAR_INDEX}")" "$(basename "${STAR_INDEX}")" \
        | aws s3 cp - "${S3_STAR_INDEX}"
    log "Upload complete in $(( SECONDS - T0 ))s."
fi

log "R01_setup_ref.sh complete."
