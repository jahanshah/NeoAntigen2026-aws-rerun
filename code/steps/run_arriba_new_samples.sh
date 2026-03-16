#!/bin/bash
# Run STAR (chimeric mode) + Arriba on D88_old, D99_old, D109_new
# These samples were previously aligned without --chimSegmentMin, so Arriba could not be run.
# This script downloads FASTQs, re-aligns with STAR chimeric flags, then runs Arriba.

set -euo pipefail
export PATH="/home/ec2-user/miniforge3/bin:$PATH"

BASE_DIR="/home/ec2-user"
REF_DIR="${BASE_DIR}/ref/mm10"
TMP_DIR="${BASE_DIR}/tmp/neoantig_pipeline"
FASTQ_DIR="${TMP_DIR}/rnaseq_fastq_new"
STAR_DIR="${TMP_DIR}/star_chimeric"
OUT_DIR="${TMP_DIR}/arriba_new"
S3_FASTQ="s3://neoantigen2026-rerun/data/fastq/rnaseq/Raw-GemOVCA-RNAseq"
S3_FUSION="s3://neoantigen2026-rerun/results/res_20260313_071626/rnaseq/fusion"

THREADS=8
STAR_INDEX="${REF_DIR}/star_index"
GENOME_FA="${REF_DIR}/mm10.fa"
GTF="${REF_DIR}/mm10_chr.gtf"
BLACKLIST="/home/ec2-user/miniforge3/var/lib/arriba/blacklist_mm10_GRCm38_v2.5.1.tsv.gz"
KNOWN_FUSIONS="/home/ec2-user/miniforge3/var/lib/arriba/known_fusions_mm10_GRCm38_v2.5.1.tsv.gz"

mkdir -p "${OUT_DIR}" "${FASTQ_DIR}" "${STAR_DIR}"

LOG="${BASE_DIR}/results/res_20260311_225555/logs/arriba_new_samples_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

echo "[$(date)] Starting STAR+Arriba for D88_old, D99_old, D109_new"

# ── Sample FASTQ locations ─────────────────────────────────────────────────
# D88_old: use the primary run (ODU060723-06 matches the existing BAM header)
process_sample() {
    local SAMPLE="$1"
    local S3_R1="$2"
    local S3_R2="$3"

    local ARRIBA_OUT="${OUT_DIR}/Arriba_${SAMPLE}.txt"
    local ARRIBA_DISC="${OUT_DIR}/Arriba_${SAMPLE}.discarded.txt"

    if [ -f "${ARRIBA_OUT}" ] && [ -s "${ARRIBA_OUT}" ]; then
        echo "[$(date)] [SKIP] ${SAMPLE} — Arriba output already exists"
        return 0
    fi

    echo "[$(date)] === Processing ${SAMPLE} ==="

    local SAMPLE_FQ="${FASTQ_DIR}/${SAMPLE}"
    local SAMPLE_STAR="${STAR_DIR}/${SAMPLE}"
    mkdir -p "${SAMPLE_FQ}" "${SAMPLE_STAR}"

    local LOCAL_R1="${SAMPLE_FQ}/$(basename ${S3_R1})"
    local LOCAL_R2="${SAMPLE_FQ}/$(basename ${S3_R2})"

    # ── Download FASTQs ──────────────────────────────────────────────────
    if [ ! -f "${LOCAL_R1}" ]; then
        echo "[$(date)]   Downloading R1: ${S3_R1}"
        aws s3 cp "${S3_R1}" "${LOCAL_R1}"
    fi
    if [ ! -f "${LOCAL_R2}" ]; then
        echo "[$(date)]   Downloading R2: ${S3_R2}"
        aws s3 cp "${S3_R2}" "${LOCAL_R2}"
    fi

    # ── STAR alignment with chimeric mode ────────────────────────────────
    local STAR_BAM="${SAMPLE_STAR}/Aligned.sortedByCoord.out.bam"

    if [ ! -f "${STAR_BAM}" ]; then
        echo "[$(date)]   Running STAR (chimeric mode)..."
        STAR \
            --runMode alignReads \
            --runThreadN ${THREADS} \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${LOCAL_R1}" "${LOCAL_R2}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${SAMPLE_STAR}/" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMstrandField intronMotif \
            --outSAMattributes NH HI AS NM MD \
            --outBAMsortingThreadN ${THREADS} \
            --outFilterMultimapNmax 20 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --alignSJoverhangMin 8 \
            --twopassMode Basic \
            --chimSegmentMin 10 \
            --chimOutType WithinBAM SoftClip \
            --chimJunctionOverhangMin 10 \
            --chimScoreMin 1 \
            --chimScoreSeparation 1 \
            --chimScoreJunctionNonGTAG 0 \
            --chimOutJunctionFormat 1

        if [ ! -f "${STAR_BAM}" ]; then
            echo "[$(date)]   ERROR: STAR failed for ${SAMPLE}"
            return 1
        fi
        samtools index -@ ${THREADS} "${STAR_BAM}"
        echo "[$(date)]   STAR done: ${STAR_BAM}"
    else
        echo "[$(date)]   STAR BAM exists: ${STAR_BAM}"
    fi

    # ── Arriba fusion calling ─────────────────────────────────────────────
    echo "[$(date)]   Running Arriba on ${SAMPLE}..."
    /home/ec2-user/miniforge3/bin/arriba \
        -x "${STAR_BAM}" \
        -a "${GENOME_FA}" \
        -g "${GTF}" \
        -b "${BLACKLIST}" \
        -k "${KNOWN_FUSIONS}" \
        -o "${ARRIBA_OUT}" \
        -O "${ARRIBA_DISC}"

    if [ -f "${ARRIBA_OUT}" ]; then
        local N
        N=$(tail -n +2 "${ARRIBA_OUT}" | wc -l)
        echo "[$(date)]   Arriba done: ${N} fusions → ${ARRIBA_OUT}"
        aws s3 cp "${ARRIBA_OUT}"  "${S3_FUSION}/Arriba_${SAMPLE}.txt" --quiet
        aws s3 cp "${ARRIBA_DISC}" "${S3_FUSION}/Arriba_${SAMPLE}.discarded.txt" --quiet
        echo "[$(date)]   Uploaded to S3"
    else
        echo "[$(date)]   ERROR: Arriba output not found for ${SAMPLE}"
        return 1
    fi

    # Clean up local FASTQs to save space
    rm -f "${LOCAL_R1}" "${LOCAL_R2}"
    echo "[$(date)]   Cleaned up local FASTQs"
}

# Run samples sequentially (heavy IO/compute)
process_sample "D88_old" \
    "${S3_FASTQ}/Day88/Sample_80_Day_88_Old/ODU060723-06_S34_L004_R1_001.fastq.gz" \
    "${S3_FASTQ}/Day88/Sample_80_Day_88_Old/ODU060723-06_S34_L004_R2_001.fastq.gz"

process_sample "D99_old" \
    "${S3_FASTQ}/D99/Sample_32_Day_99_Old-1/ODU060723-01_S29_L004_R1_001.fastq.gz" \
    "${S3_FASTQ}/D99/Sample_32_Day_99_Old-1/ODU060723-01_S29_L004_R2_001.fastq.gz"

# D109_new: find FASTQs dynamically
D109_R1=$(aws s3 ls "${S3_FASTQ}/Sample_2661_Day_109_New/" | grep _R1_ | awk '{print $4}' | head -1)
D109_R2=$(aws s3 ls "${S3_FASTQ}/Sample_2661_Day_109_New/" | grep _R2_ | awk '{print $4}' | head -1)
process_sample "D109_new" \
    "${S3_FASTQ}/Sample_2661_Day_109_New/${D109_R1}" \
    "${S3_FASTQ}/Sample_2661_Day_109_New/${D109_R2}"

echo "[$(date)] All samples complete. Arriba outputs in: ${OUT_DIR}"
