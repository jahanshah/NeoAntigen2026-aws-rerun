#!/bin/bash
# Sequenza pipeline: bam2seqz + seqz_binning + R analysis
# Runs for all 5 available tumor samples using matched normal 423_D0_old
set -euo pipefail
export PATH="/home/ec2-user/miniforge3/bin:$PATH"

RESULTS_DIR="/home/ec2-user/results/res_20260311_225555"
SEQ_DIR="${RESULTS_DIR}/sequenza"
TMP_DIR="/home/ec2-user/tmp/neoantig_pipeline/sequenza_bam"
S3_WES="s3://neoantigen2026-rerun/pre-processed/wes"
S3_OUT="s3://neoantigen2026-rerun/results/res_20260311_225555/sequenza"
REF="/home/ec2-user/ref/mm10/mm10.fa"
GC_WIG="/home/ec2-user/ref/mm10/mm10.gc50.wig.gz"
THREADS=8
LOG="${RESULTS_DIR}/logs/sequenza_$(date +%Y%m%d_%H%M%S).log"

mkdir -p "${SEQ_DIR}" "${TMP_DIR}"
exec > >(tee -a "${LOG}") 2>&1
echo "[$(date)] === Sequenza Purity Pipeline ==="

# ── Download normal BAM ─────────────────────────────────────────────────────
NORMAL_BAM="${TMP_DIR}/423_D0_old.preproc.bam"
NORMAL_BAI="${TMP_DIR}/423_D0_old.preproc.bai"
if [ ! -f "${NORMAL_BAM}" ]; then
    echo "[$(date)] Downloading normal BAM..."
    aws s3 cp "${S3_WES}/423_D0_old.preproc.bam" "${NORMAL_BAM}"
    aws s3 cp "${S3_WES}/423_D0_old.preproc.bai" "${NORMAL_BAI}"
fi

# ── Sample map: sample_id → use_local_bam? ─────────────────────────────────
declare -A LOCAL_BAMS=(
    ["428_D20_new"]="/home/ec2-user/tmp/neoantig_pipeline/mutect2/428_D20_new.preproc.bam"
    ["443_D21_new"]="/home/ec2-user/tmp/neoantig_pipeline/mutect2/443_D21_new.preproc.bam"
)
SAMPLES="36_D99_new 38_D99_new 42_D122_old 428_D20_new 443_D21_new"

for SAMPLE in ${SAMPLES}; do
    SEQZ_RAW="${SEQ_DIR}/${SAMPLE}.seqz.gz"
    SEQZ_BIN="${SEQ_DIR}/${SAMPLE}.binned.seqz.gz"
    SAMPLE_OUTDIR="${SEQ_DIR}/${SAMPLE}"

    echo "[$(date)] === ${SAMPLE} ==="

    # Get tumor BAM
    if [[ -n "${LOCAL_BAMS[$SAMPLE]+set}" ]] && [ -f "${LOCAL_BAMS[$SAMPLE]}" ]; then
        TUMOR_BAM="${LOCAL_BAMS[$SAMPLE]}"
        TUMOR_BAI="${TUMOR_BAM%.bam}.bai"
        # samtools needs .bam.bai or .bai — check
        [ ! -f "${TUMOR_BAI}" ] && TUMOR_BAI="${TUMOR_BAM%.preproc.bam}.preproc.bai"
    else
        TUMOR_BAM="${TMP_DIR}/${SAMPLE}.preproc.bam"
        TUMOR_BAI="${TMP_DIR}/${SAMPLE}.preproc.bai"
        if [ ! -f "${TUMOR_BAM}" ]; then
            echo "[$(date)]   Downloading ${SAMPLE}..."
            aws s3 cp "${S3_WES}/${SAMPLE}.preproc.bam" "${TUMOR_BAM}"
            aws s3 cp "${S3_WES}/${SAMPLE}.preproc.bai" "${TUMOR_BAI}"
        fi
    fi
    echo "[$(date)]   Tumor BAM: ${TUMOR_BAM}"

    # ── bam2seqz ──────────────────────────────────────────────────────────
    if [ ! -f "${SEQZ_RAW}" ]; then
        echo "[$(date)]   Running bam2seqz..."
        python3 -c "
import sys
sys.argv = ['bam2seqz',
    '-n', '${NORMAL_BAM}',
    '-t', '${TUMOR_BAM}',
    '-F', '${REF}',
    '--gc_file', '${GC_WIG}',
    '-o', '${SEQZ_RAW}',
    '--parallel', '${THREADS}']

import argparse
subparsers = argparse.ArgumentParser().add_subparsers()
from sequenza.programs.bam2seqz import bam2seqz_main, bam2seqz_args, add_parser
from sequenza.misc import SeqzLogger
log = SeqzLogger(level=30)
parser = add_parser(subparsers, 'bam2seqz')
args = bam2seqz_args(parser, sys.argv[1:])
bam2seqz_main(args, log)
"
        echo "[$(date)]   bam2seqz done: $(ls -lh ${SEQZ_RAW} | awk '{print \$5}')"
    else
        echo "[$(date)]   seqz exists, skipping bam2seqz"
    fi

    # ── seqz_binning ──────────────────────────────────────────────────────
    if [ ! -f "${SEQZ_BIN}" ]; then
        echo "[$(date)]   Running seqz_binning (300bp windows)..."
        python3 -c "
import sys
sys.argv = ['seqz_binning',
    '--seqz', '${SEQZ_RAW}',
    '-w', '300',
    '-o', '${SEQZ_BIN}']

import argparse
subparsers = argparse.ArgumentParser().add_subparsers()
from sequenza.programs.seqz_binning import seqz_binning, add_parser
from sequenza.misc import SeqzLogger
log = SeqzLogger(level=30)
parser = add_parser(subparsers, 'seqz_binning')

import argparse as ap2
p2 = ap2.ArgumentParser()
p2.add_argument('--seqz', required=True)
p2.add_argument('-w', '--window', type=int, default=300)
p2.add_argument('-o', '--output', required=True)
args = p2.parse_args(sys.argv[1:])
seqz_binning(None, 'seqz_binning', sys.argv[1:], log)
"
        echo "[$(date)]   Binned seqz: $(ls -lh ${SEQZ_BIN} | awk '{print \$5}')"
    fi

    # ── R sequenza analysis ───────────────────────────────────────────────
    mkdir -p "${SAMPLE_OUTDIR}"
    echo "[$(date)]   Running R sequenza analysis..."
    Rscript - "${SEQZ_BIN}" "${SAMPLE_OUTDIR}" "${SAMPLE}" << 'REOF'
args    <- commandArgs(trailingOnly=TRUE)
seqz_f  <- args[1]
outdir  <- args[2]
sname   <- args[3]

suppressPackageStartupMessages(library(sequenza))

cat(sprintf("[R] Loading seqz: %s\n", seqz_f))
seqzdata <- sequenza.extract(seqz_f, verbose=FALSE)

cat("[R] Fitting purity/ploidy grid...\n")
CP <- sequenza.fit(seqzdata)

cat("[R] Writing results...\n")
sequenza.results(
    sequenza.extract = seqzdata,
    cp.table         = CP,
    sample.id        = sname,
    out.dir          = outdir
)

# Extract and print best estimate
best_idx <- which.max(CP$log.posterior)
purity <- CP$values.matrix[best_idx, "cellularity"]
ploidy <- CP$values.matrix[best_idx, "ploidy"]
cat(sprintf("SEQUENZA_RESULT\t%s\tpurity=%.4f\tploidy=%.2f\n", sname, purity, ploidy))
REOF

    echo "[$(date)]   Uploading to S3..."
    aws s3 cp "${SAMPLE_OUTDIR}/" "${S3_OUT}/${SAMPLE}/" --recursive --quiet
    echo "[$(date)]   ${SAMPLE} complete"
done

echo "[$(date)] === All samples complete ==="
echo "[$(date)] Results: ${SEQ_DIR}"
