#!/usr/bin/env bash
# =============================================================================
# setup_variant_calling.sh
# Installs GATK4 and downloads mm10 reference genome for variant calling.
# Run once before variant_calling.sh
# =============================================================================

set -euo pipefail

CONDA="/home/ec2-user/miniforge3/bin/conda"
REF_DIR="/home/ec2-user/ref/mm10"
S3_RNASEQ="s3://bam-wes/NeoAntigen-aws/data/RNASeq"

mkdir -p "${REF_DIR}"

echo "========================================"
echo "[STEP 1] Install GATK4 via conda"
echo "========================================"
export GATK_LOCAL_JAR="/home/ec2-user/miniforge3/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"
if ! /home/ec2-user/miniforge3/bin/gatk --version &>/dev/null; then
    "${CONDA}" install -y -c bioconda -c conda-forge gatk4 2>&1 | tail -5
else
    echo "[SKIP] GATK4 already installed: $(/home/ec2-user/miniforge3/bin/gatk --version 2>&1 | head -1)"
fi

echo ""
echo "========================================"
echo "[STEP 2] Download mm10 reference FASTA from UCSC"
echo "========================================"
FA="${REF_DIR}/mm10.fa"
if [[ -f "${FA}" ]]; then
    echo "[SKIP] Reference already exists: ${FA}"
else
    FA_GZ="${REF_DIR}/mm10.fa.gz"
    echo "[INFO] Downloading mm10.fa.gz (~800 MB)..."
    curl -L --progress-bar \
        "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz" \
        -o "${FA_GZ}"
    echo "[INFO] Decompressing..."
    gunzip "${FA_GZ}"          # produces mm10.fa
    echo "[INFO] Reference FASTA: ${FA}"
fi

echo ""
echo "========================================"
echo "[STEP 3] Download dict and fai from S3"
echo "========================================"
for EXT in dict fa.fai; do
    DEST="${REF_DIR}/mm10.${EXT}"
    if [[ -f "${DEST}" ]]; then
        echo "[SKIP] ${DEST} already exists"
    else
        aws s3 cp "${S3_RNASEQ}/mm10.${EXT}" "${DEST}"
        echo "[INFO] Downloaded: ${DEST}"
    fi
done

echo ""
echo "========================================"
echo "[STEP 4] Verify reference files"
echo "========================================"
ls -lh "${REF_DIR}/"
echo ""
echo "[INFO] Setup complete. Ready to run variant_calling.sh"
