#!/usr/bin/env bash
# =============================================================================
# setup_aws_environment.sh
# =============================================================================
# One-time environment setup for the NeoAntigen pipeline on AWS Linux
# (Amazon Linux 2 / Ubuntu 22.04)
#
# Usage:
#   bash scripts/setup_aws_environment.sh [--refs-only] [--tools-only]
#
# What it installs:
#   - System packages (gcc, zlib, java, etc.)
#   - Miniconda3 (if not present)
#   - Conda environments: gatk-env, r-env, pyclone-env, snpeff-env
#   - Python packages: pysam, requests, pandas, biopython
#   - R packages: DESeq2, edgeR, limma, ggplot2, pheatmap, patchwork, ...
#   - Reference files: mm10.fa, gencode.vM25 GTF, proteins.fa, exon BED
#
# Estimated time: 45–90 min (mainly R Bioconductor + reference downloads)
# Disk required: ~50 GB (references + tools)
# =============================================================================

set -euo pipefail

WDIR="${WDIR:-/data/neoantigen}"
LOG="$WDIR/logs/setup.log"
REFS="$WDIR/reference_db"
DATA="$WDIR/data"

mkdir -p "$WDIR/logs" "$REFS" "$DATA/bam" "$DATA/rnaseq_bam" \
         "$WDIR/results" "$WDIR/figures" "$WDIR/reports" \
         "$WDIR/scripts" "$WDIR/peptides" "$WDIR/peptides_frameshift" \
         "$WDIR/peptides_fusion" "$WDIR/iedb_scores" \
         "$WDIR/pyclone_run/yaml" "$WDIR/pyclone_run/trace" \
         "$WDIR/pyclone_run/tables" "$WDIR/.done"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$LOG"; }

# ── Parse flags ───────────────────────────────────────────────────────────────
REFS_ONLY=false
TOOLS_ONLY=false
for arg in "$@"; do
  case $arg in
    --refs-only)  REFS_ONLY=true ;;
    --tools-only) TOOLS_ONLY=true ;;
  esac
done

# =============================================================================
# SECTION 1 — System packages
# =============================================================================
install_system_packages() {
    log "── System packages ──────────────────────────────"

    # Detect OS
    if grep -qi "ubuntu" /etc/os-release 2>/dev/null; then
        OS="ubuntu"
    elif grep -qi "amazon\|centos\|rhel" /etc/os-release 2>/dev/null; then
        OS="amazon"
    else
        OS="unknown"
    fi

    if [[ "$OS" == "ubuntu" ]]; then
        sudo apt-get update -qq
        sudo apt-get install -y -qq \
            wget curl git unzip bzip2 \
            gcc g++ make cmake \
            zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
            libncurses5-dev libreadline-dev \
            default-jre default-jdk \
            r-base \
            parallel pigz \
            awscli 2>>"$LOG"
    elif [[ "$OS" == "amazon" ]]; then
        sudo yum update -y -q
        sudo yum install -y -q \
            wget curl git unzip bzip2 \
            gcc gcc-c++ make cmake \
            zlib-devel bzip2-devel xz-devel libcurl-devel openssl-devel \
            ncurses-devel readline-devel \
            java-11-amazon-corretto \
            parallel pigz \
            awscli 2>>"$LOG"
    fi
    log "  System packages installed (OS: $OS)"
}

# =============================================================================
# SECTION 2 — Miniconda
# =============================================================================
install_miniconda() {
    log "── Miniconda ────────────────────────────────────"
    if command -v conda &>/dev/null; then
        log "  conda already installed: $(conda --version)"
        return
    fi

    MINICONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    wget -q "$MINICONDA_URL" -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
    rm -f /tmp/miniconda.sh

    # Initialize
    "$HOME/miniconda3/bin/conda" init bash
    source "$HOME/.bashrc" 2>/dev/null || true
    export PATH="$HOME/miniconda3/bin:$PATH"

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    conda update -y conda 2>>"$LOG"
    log "  Miniconda installed"
}

# =============================================================================
# SECTION 3 — Conda environments
# =============================================================================
install_conda_envs() {
    log "── Conda environments ───────────────────────────"

    # ── 3a. Main bioinformatics env ───────────────────────────────────────────
    log "  Creating gatk-env ..."
    conda create -y -n gatk-env -c bioconda -c conda-forge \
        gatk4=4.4.0.0 \
        picard=3.1.1 \
        samtools=1.19 \
        bcftools=1.19 \
        bedtools=2.31 \
        bwa=0.7.17 \
        2>>"$LOG" || log "  gatk-env: some packages may already exist"

    # ── 3b. SnpEff env ────────────────────────────────────────────────────────
    log "  Creating snpeff-env ..."
    conda create -y -n snpeff-env -c bioconda -c conda-forge \
        snpeff=5.1 \
        2>>"$LOG" || true

    # Download GRCm38.86 SnpEff database
    conda run -n snpeff-env snpEff download GRCm38.86 \
        2>>"$LOG" || log "  SnpEff DB: download manually if this fails"

    # ── 3c. PyClone env (Python 2.7) ──────────────────────────────────────────
    log "  Creating pyclone-env ..."
    conda create -y -n pyclone-env -c bioconda -c conda-forge \
        python=2.7 pyclone=0.13.1 \
        2>>"$LOG" || true

    # ── 3d. Python analysis env ───────────────────────────────────────────────
    log "  Creating python-env ..."
    conda create -y -n python-env -c bioconda -c conda-forge \
        python=3.11 \
        pysam=0.22 \
        pandas=2.1 \
        numpy=1.26 \
        scipy=1.11 \
        biopython=1.83 \
        requests=2.31 \
        openpyxl=3.1 \
        matplotlib=3.8 \
        seaborn=0.13 \
        2>>"$LOG" || true

    # ── 3e. R environment ─────────────────────────────────────────────────────
    log "  Creating r-env (this takes 15-30 min) ..."
    conda create -y -n r-env -c conda-forge -c bioconda \
        r-base=4.3.2 \
        r-ggplot2=3.4.4 \
        r-dplyr=1.1.4 \
        r-tidyr=1.3.0 \
        r-stringr=1.5.1 \
        r-readr=2.1.4 \
        r-ggrepel=0.9.4 \
        r-scales=1.3.0 \
        r-patchwork=1.1.3 \
        r-cowplot=1.1.3 \
        r-pheatmap=1.0.12 \
        r-rcolorbrewer=1.1_3 \
        r-ggpubr=0.6.0 \
        r-viridis=0.6.4 \
        r-gridextra=2.3 \
        r-knitr=1.45 \
        r-rmarkdown=2.25 \
        bioconductor-deseq2=1.42.0 \
        bioconductor-edger=4.0.2 \
        bioconductor-limma=3.58.1 \
        bioconductor-biocparallel=1.36.0 \
        bioconductor-complexheatmap=2.18.0 \
        bioconductor-apeglm=1.24.0 \
        2>>"$LOG" || log "  r-env: some packages may need manual install"

    log "  Conda environments complete"
}

# =============================================================================
# SECTION 4 — Python packages (pip into python-env)
# =============================================================================
install_python_packages() {
    log "── Python packages (pip) ────────────────────────"
    conda run -n python-env pip install -q \
        tqdm \
        pyaml \
        openpyxl \
        xlrd \
        2>>"$LOG"
    log "  Python pip packages installed"
}

# =============================================================================
# SECTION 5 — Reference files
# =============================================================================
download_references() {
    log "── Reference files ──────────────────────────────"
    cd "$REFS"

    # ── 5a. mm10 reference genome (UCSC) ──────────────────────────────────────
    if [[ ! -f "$DATA/mm10.fa" ]]; then
        log "  Downloading mm10.fa ..."
        wget -q "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz" \
            -O mm10.fa.gz
        gunzip mm10.fa.gz
        mv mm10.fa "$DATA/mm10.fa"

        log "  Building BWA index (takes ~60 min) ..."
        conda run -n gatk-env bwa index -a bwtsw "$DATA/mm10.fa" \
            2>>"$LOG" &
        BWA_PID=$!

        log "  Building samtools fai + dict ..."
        conda run -n gatk-env samtools faidx "$DATA/mm10.fa"
        conda run -n gatk-env picard CreateSequenceDictionary \
            R="$DATA/mm10.fa" O="$DATA/mm10.dict" 2>>"$LOG"

        wait $BWA_PID
        log "  mm10 indexed"
    else
        log "  mm10.fa already present"
    fi

    # ── 5b. GENCODE vM25 annotation (GRCm38) ──────────────────────────────────
    if [[ ! -f "$REFS/gencode.vM25.annotation.gtf.gz" ]]; then
        log "  Downloading GENCODE vM25 GTF ..."
        wget -q "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz" \
            -O "$REFS/gencode.vM25.annotation.gtf.gz"
    fi

    # ── 5c. Ensembl GRCm38 proteome ───────────────────────────────────────────
    if [[ ! -f "$REFS/Mus_musculus.GRCm38.pep.all.fa.gz" ]]; then
        log "  Downloading Ensembl GRCm38 proteome ..."
        wget -q "https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz" \
            -O "$REFS/Mus_musculus.GRCm38.pep.all.fa.gz"
    fi

    # Build proteins.fa (unmasked, primary transcripts)
    if [[ ! -f "$REFS/proteins.fa" ]]; then
        log "  Building proteins.fa ..."
        zcat "$REFS/Mus_musculus.GRCm38.pep.all.fa.gz" \
            | awk '/^>/{p=($0 !~ /PATCH|HSCHR/)} p' \
            > "$REFS/proteins.fa"
    fi

    # ── 5d. Exon BED (merged, protein-coding) ─────────────────────────────────
    if [[ ! -f "$REFS/merged_exons.bed" ]]; then
        log "  Building merged_exons.bed ..."
        zcat "$REFS/gencode.vM25.annotation.gtf.gz" \
            | awk '$3=="exon" && /gene_type "protein_coding"/ {
                  gsub("chr",""); print $1"\t"$4-1"\t"$5}' \
            | sort -k1,1 -k2,2n \
            | conda run -n gatk-env bedtools merge -i stdin \
            > "$REFS/merged_exons.bed"
    fi

    # ── 5e. Mouse gnomAD-equivalent (MGP SNP DB for Mutect2) ─────────────────
    if [[ ! -f "$REFS/mgp_REL2021_snps.vcf.gz" ]]; then
        log "  Downloading MGP SNP database ..."
        wget -q "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz" \
            -O "$REFS/mgp_REL2021_snps.vcf.gz" || \
        log "  MGP DB: download failed — set GNOMAD_MM10 path manually in pipeline config"

        [[ -f "$REFS/mgp_REL2021_snps.vcf.gz" ]] && \
            conda run -n gatk-env bcftools index "$REFS/mgp_REL2021_snps.vcf.gz"
    fi

    log "  References complete"
}

# =============================================================================
# SECTION 6 — Verify installation
# =============================================================================
verify_install() {
    log "── Verification ─────────────────────────────────"
    local ok=true

    check() {
        local name="$1" cmd="$2"
        if eval "$cmd" &>/dev/null; then
            log "  ✓ $name"
        else
            log "  ✗ $name — FAILED"
            ok=false
        fi
    }

    check "gatk"       "conda run -n gatk-env gatk --version"
    check "picard"     "conda run -n gatk-env picard --version"
    check "samtools"   "conda run -n gatk-env samtools --version"
    check "bcftools"   "conda run -n gatk-env bcftools --version"
    check "snpEff"     "conda run -n snpeff-env snpEff -version"
    check "PyClone"    "conda run -n pyclone-env PyClone --version"
    check "python"     "conda run -n python-env python --version"
    check "pysam"      "conda run -n python-env python -c 'import pysam'"
    check "pandas"     "conda run -n python-env python -c 'import pandas'"
    check "Rscript"    "conda run -n r-env Rscript --version"
    check "ggplot2"    "conda run -n r-env Rscript -e 'library(ggplot2)'"
    check "DESeq2"     "conda run -n r-env Rscript -e 'library(DESeq2)'"
    check "edgeR"      "conda run -n r-env Rscript -e 'library(edgeR)'"
    check "mm10.fa"    "test -f '$DATA/mm10.fa'"
    check "GTF"        "test -f '$REFS/gencode.vM25.annotation.gtf.gz'"
    check "proteins"   "test -f '$REFS/proteins.fa'"
    check "exon BED"   "test -f '$REFS/merged_exons.bed'"

    $ok && log "  All checks passed" || log "  Some checks failed — review log: $LOG"
}

# =============================================================================
# MAIN
# =============================================================================
log "══════════════════════════════════════════════════"
log " NeoAntigen AWS Environment Setup"
log " $(date)"
log " WDIR: $WDIR"
log "══════════════════════════════════════════════════"

if ! $REFS_ONLY; then
    install_system_packages
    install_miniconda
    install_conda_envs
    install_python_packages
fi

if ! $TOOLS_ONLY; then
    download_references
fi

verify_install

log ""
log "══════════════════════════════════════════════════"
log " Setup complete — $(date)"
log " Next: edit WDIR and sample paths in"
log "       scripts/aws_master_pipeline.sh"
log " Then: bash scripts/aws_master_pipeline.sh"
log "══════════════════════════════════════════════════"
