#!/usr/bin/env bash
# =============================================================================
# setup_pipeline.sh
# Install all pipeline dependencies and download reference files.
# Run once before run_pipeline.sh
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

MAMBA="${CONDA_BIN}/mamba"

# =============================================================================
# SECTION 1 — Bioinformatics tools via mamba
# =============================================================================
log "Installing bioinformatics tools..."
for pkg in samtools bcftools bwa gatk4 picard snpeff netmhcpan; do
    if ! "${CONDA_BIN}/${pkg}" --version &>/dev/null && \
       ! "${CONDA_BIN}/$(echo $pkg | sed 's/gatk4/gatk/;s/snpeff/snpEff/;s/netmhcpan/netMHCpan/')" --version &>/dev/null; then
        log "  Installing: ${pkg}"
        ${MAMBA} install -y -c bioconda -c conda-forge "${pkg}" 2>&1 | tail -3
    else
        log "  [SKIP] ${pkg} already installed"
    fi
done

# Fix GATK python shebang if needed
if head -1 "${GATK}" | grep -q "#!/usr/bin/env python$"; then
    sed -i 's|#!/usr/bin/env python$|#!/usr/bin/env python3|' "${GATK}"
    log "Fixed GATK python shebang"
fi
log "GATK: $(${GATK} --version 2>&1 | grep 'The Genome' || echo 'check GATK_LOCAL_JAR')"

# =============================================================================
# SECTION 2 — Python packages
# =============================================================================
log "Installing Python packages..."
${PYTHON} -m pip install --quiet --upgrade pip
${PYTHON} -m pip install --quiet pandas pysam pyensembl pyvcf3 pyyaml

# Install PyClone (original — beta-binomial MCMC)
if ! ${PYTHON} -c "import pyclone" &>/dev/null; then
    log "Installing PyClone..."
    ${PYTHON} -m pip install --quiet PyClone
fi
log "PyClone: $(${PYTHON} -c 'import pyclone; print(pyclone.__version__)' 2>/dev/null || echo 'installed')"

# =============================================================================
# SECTION 3 — Reference genome (mm10)
# =============================================================================
mkdir -p "${REF_DIR}"

# 3a. FASTA
if [[ -f "${REF}" ]]; then
    log "[SKIP] mm10.fa exists ($(du -sh ${REF} | cut -f1))"
else
    log "Downloading mm10.fa.gz from UCSC (~800 MB)..."
    curl -L --progress-bar \
        "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz" \
        -o "${REF}.gz"
    gunzip "${REF}.gz"
    log "Reference decompressed: ${REF}"
fi

# 3b. dict + fai from S3
for EXT in dict fa.fai; do
    DEST="${REF_DIR}/mm10.${EXT}"
    [[ -f "${DEST}" ]] && { log "[SKIP] mm10.${EXT}"; continue; }
    aws s3 cp "${S3_ROOT}/data/RNASeq/mm10.${EXT}" "${DEST}"
    log "Downloaded: ${DEST}"
done

# =============================================================================
# SECTION 4 — Known SNPs for BQSR (Mouse Genome Project)
# =============================================================================
if [[ -f "${KNOWN_SNPS}" ]]; then
    log "[SKIP] MGP SNPs exist"
else
    log "Downloading MGP v5 SNPs for BQSR..."
    MGP_URL="ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
    curl -L --progress-bar "${MGP_URL}" -o "${KNOWN_SNPS}" || {
        log "[WARN] MGP download failed — BQSR will use --run-without-dbsnp-potentially-ruining-quality"
        touch "${KNOWN_SNPS}.missing"
    }
    if [[ -f "${KNOWN_SNPS}" ]]; then
        ${BCFTOOLS} index -t "${KNOWN_SNPS}"
        # Ensure contig names match (UCSC chr-style)
        FIRST_CHR=$(${BCFTOOLS} view -h "${KNOWN_SNPS}" | grep "^##contig" | head -1 | grep -o 'ID=[^,>]*' | cut -d= -f2)
        if [[ "${FIRST_CHR}" != chr* ]]; then
            log "Adding chr prefix to MGP VCF contigs..."
            ${BCFTOOLS} annotate --rename-chrs <(for i in $(seq 1 19) X Y MT; do echo -e "${i}\tchr${i}"; done) \
                -o "${KNOWN_SNPS%.gz}.chrfix.vcf.gz" -O z "${KNOWN_SNPS}"
            mv "${KNOWN_SNPS%.gz}.chrfix.vcf.gz" "${KNOWN_SNPS}"
            ${BCFTOOLS} index -t "${KNOWN_SNPS}"
        fi
    fi
fi

# =============================================================================
# SECTION 5 — mm10 exon BED (for SelectVariants WES filtering)
# =============================================================================
if [[ -f "${EXON_BED}" ]]; then
    log "[SKIP] Exon BED exists"
else
    log "Building mm10 exon BED from Ensembl GTF..."
    GTF_GZ="${REF_DIR}/Mus_musculus.GRCm38.86.gtf.gz"
    if [[ ! -f "${GTF_GZ}" ]]; then
        curl -L --progress-bar \
            "ftp://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.gtf.gz" \
            -o "${GTF_GZ}"
    fi
    # Extract exon intervals, add chr prefix, sort, merge
    zcat "${GTF_GZ}" | awk '$3=="exon"' \
        | awk '{OFS="\t"; print $1,$4-1,$5}' \
        | awk '!/^#/' \
        | sed 's/^/chr/' \
        | sort -k1,1 -k2,2n \
        | ${CONDA_BIN}/bedtools merge -i - \
        > "${EXON_BED}"
    log "Exon BED: $(wc -l < ${EXON_BED}) regions"
fi

# =============================================================================
# SECTION 6 — SnpEff GRCm38.86 database
# =============================================================================
log "Checking SnpEff database..."
if ! ${SNPEFF} databases 2>/dev/null | grep -q "${SNPEFF_DB}"; then
    log "Downloading SnpEff ${SNPEFF_DB} database..."
    ${SNPEFF} download "${SNPEFF_DB}" 2>&1 | tail -3
else
    log "[SKIP] SnpEff ${SNPEFF_DB} database exists"
fi

# =============================================================================
# SECTION 7 — netMHCpan (verify installation)
# =============================================================================
log "Checking netMHCpan..."
if ! command -v netMHCpan &>/dev/null && [[ ! -f "${NETMHCPAN}" ]]; then
    log "[WARN] netMHCpan not found. Install manually from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
    log "       Then place the binary at ${NETMHCPAN}"
else
    log "netMHCpan: OK"
fi

# =============================================================================
# SECTION 8 — PyClone config directory
# =============================================================================
mkdir -p "${RESULTS_DIR}/pyclone/input" "${RESULTS_DIR}/pyclone/output"

log ""
log "========================================"
log "Setup complete. Run: bash /home/ec2-user/code/run_pipeline.sh"
log "========================================"
log "Reference: ${REF}"
log "Known SNPs: ${KNOWN_SNPS}"
log "Exon BED:  ${EXON_BED}"
