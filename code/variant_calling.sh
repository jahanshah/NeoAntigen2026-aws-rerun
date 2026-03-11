#!/usr/bin/env bash
# =============================================================================
# variant_calling.sh
# Somatic SNV/Indel calling with GATK4 Mutect2 (mm10/GRCm38, WES)
#
# Tumor-Normal pairing:
#   Normal : 423_D0_old (SM: D0_old)  -- matched normal baseline
#   PoN    : 423_D0_old + 424_D0__old -- Panel of Normals
#   Tumors : 34_D52_old, 36_D99_new, 38_D99_new, 428_D20_new, 42_D122_old
#
# Steps:
#   3a. Create Panel of Normals (PoN)
#   3b. Run Mutect2 per tumor (tumor vs matched normal + PoN)
#   3c. FilterMutectCalls
#   3d. Summarize with R
# =============================================================================

set -euo pipefail

# --- Config ------------------------------------------------------------------
S3_BAM="s3://bam-wes/NeoAntigen-aws/data/bam"
REF="/home/ec2-user/ref/mm10/mm10.fa"
OUTDIR="/home/ec2-user/results/variant_calling"
TMPDIR="/tmp/gatk_vc"
export GATK_LOCAL_JAR="/home/ec2-user/miniforge3/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar"
GATK="/home/ec2-user/miniforge3/bin/gatk"
SAMTOOLS="/home/ec2-user/miniforge3/bin/samtools"
THREADS=4

# Tumor-Normal manifest (tumor_bam  normal_sm  normal_bam)
declare -A TUMOR_SAMPLES=(
    ["34_D52_old"]="D0_old"
    ["36_D99_new"]="D0_old"
    ["38_D99_new"]="D0_old"
    ["428_D20_new"]="D0_old"
    ["42_D122_old"]="D0_old"
)
NORMAL_BAM="423_D0_old"
NORMAL_SM="D0_old"
PON_NORMALS=("423_D0_old" "424_D0__old")

mkdir -p "${OUTDIR}/pon" "${OUTDIR}/vcf_raw" "${OUTDIR}/vcf_filtered" \
         "${OUTDIR}/stats" "${TMPDIR}"

echo "[INFO] GATK: $(${GATK} --version 2>&1 | head -1)"
echo "[INFO] Reference: ${REF}"

# --- Verify reference --------------------------------------------------------
[[ -f "${REF}" ]]        || { echo "[ERROR] Reference not found: ${REF}"; exit 1; }
[[ -f "${REF}.fai" ]]    || { echo "[ERROR] Missing: ${REF}.fai"; exit 1; }
[[ -f "${REF%.fa}.dict" ]] || { echo "[ERROR] Missing: ${REF%.fa}.dict"; exit 1; }

# =============================================================================
# STEP 3a — Create Panel of Normals (PoN)
# =============================================================================
PON_DB="${OUTDIR}/pon/pon_db"
PON_VCF="${OUTDIR}/pon/pon.vcf.gz"

if [[ -f "${PON_VCF}" ]]; then
    echo "[SKIP] PoN already exists: ${PON_VCF}"
else
    echo ""
    echo "========================================"
    echo "[STEP 3a] Creating Panel of Normals"
    echo "========================================"

    PON_MUTECT2_ARGS=""
    for NORM_SAMPLE in "${PON_NORMALS[@]}"; do
        LOCAL_BAM="${TMPDIR}/${NORM_SAMPLE}.bam"
        LOCAL_BAI="${LOCAL_BAM}.bai"
        NORM_VCF="${OUTDIR}/pon/${NORM_SAMPLE}.vcf.gz"

        if [[ ! -f "${NORM_VCF}" ]]; then
            echo "  Downloading ${NORM_SAMPLE}.bam..."
            aws s3 cp "${S3_BAM}/${NORM_SAMPLE}.bam"     "${LOCAL_BAM}"
            aws s3 cp "${S3_BAM}/${NORM_SAMPLE}.bam.bai" "${LOCAL_BAI}"

            echo "  Running Mutect2 (tumor-only mode for PoN): ${NORM_SAMPLE}"
            ${GATK} Mutect2 \
                -R "${REF}" \
                -I "${LOCAL_BAM}" \
                -O "${NORM_VCF}" \
                --max-mnp-distance 0 \
                --native-pair-hmm-threads ${THREADS} 2>&1

            rm -f "${LOCAL_BAM}" "${LOCAL_BAI}"
        else
            echo "  [SKIP] PoN Mutect2 already done: ${NORM_SAMPLE}"
        fi
        PON_MUTECT2_ARGS="${PON_MUTECT2_ARGS} -V ${NORM_VCF}"
    done

    echo "  Running GenomicsDBImport for PoN..."
    rm -rf "${PON_DB}"
    ${GATK} GenomicsDBImport \
        -R "${REF}" \
        ${PON_MUTECT2_ARGS} \
        --genomicsdb-workspace-path "${PON_DB}" \
        -L chr1 -L chr2 -L chr3 -L chr4 -L chr5 \
        -L chr6 -L chr7 -L chr8 -L chr9 -L chr10 \
        -L chr11 -L chr12 -L chr13 -L chr14 -L chr15 \
        -L chr16 -L chr17 -L chr18 -L chr19 \
        -L chrX -L chrY 2>&1

    echo "  Creating PoN VCF..."
    ${GATK} CreateSomaticPanelOfNormals \
        -R "${REF}" \
        -V "gendb://${PON_DB}" \
        -O "${PON_VCF}" 2>&1

    echo "  PoN created: ${PON_VCF}"
fi

# =============================================================================
# STEP 3b — Download normal once, then run Mutect2 per tumor
# =============================================================================
echo ""
echo "========================================"
echo "[STEP 3b] Downloading matched normal: ${NORMAL_BAM}"
echo "========================================"
LOCAL_NORMAL="${TMPDIR}/${NORMAL_BAM}.bam"
if [[ ! -f "${LOCAL_NORMAL}" ]]; then
    aws s3 cp "${S3_BAM}/${NORMAL_BAM}.bam"     "${LOCAL_NORMAL}"
    aws s3 cp "${S3_BAM}/${NORMAL_BAM}.bam.bai" "${LOCAL_NORMAL}.bai"
    echo "  Normal downloaded: ${LOCAL_NORMAL}"
else
    echo "  [SKIP] Normal already local"
fi

S3_RESULTS="s3://bam-wes/NeoAntigen-aws/results/variant_calling"

RESULTS_TSV="${OUTDIR}/variant_calling_results.tsv"
printf "sample\tnormal\traw_vcf\tfiltered_vcf\ttotal_variants\tpass_variants\tsnv_pass\tindel_pass\tstatus\n" \
    > "${RESULTS_TSV}"

for TUMOR_SAMPLE in "${!TUMOR_SAMPLES[@]}"; do
    NORMAL_SM_TAG="${TUMOR_SAMPLES[$TUMOR_SAMPLE]}"
    RAW_VCF="${OUTDIR}/vcf_raw/${TUMOR_SAMPLE}.vcf.gz"
    FILT_VCF="${OUTDIR}/vcf_filtered/${TUMOR_SAMPLE}.filtered.vcf.gz"
    STATS="${OUTDIR}/stats/${TUMOR_SAMPLE}.mutect2.stats"

    echo ""
    echo "========================================"
    echo "[STEP 3b] Mutect2: ${TUMOR_SAMPLE} vs ${NORMAL_BAM}"
    echo "========================================"

    STATUS="PASS"
    LOCAL_TUMOR="${TMPDIR}/${TUMOR_SAMPLE}.bam"

    # Download tumor
    if [[ ! -f "${LOCAL_TUMOR}" ]]; then
        echo "  Downloading tumor BAM: ${TUMOR_SAMPLE}..."
        aws s3 cp "${S3_BAM}/${TUMOR_SAMPLE}.bam"     "${LOCAL_TUMOR}"
        aws s3 cp "${S3_BAM}/${TUMOR_SAMPLE}.bam.bai" "${LOCAL_TUMOR}.bai"
    fi

    # --- Mutect2 -------------------------------------------------------------
    if [[ -f "${RAW_VCF}" ]]; then
        echo "  [SKIP] Raw VCF exists: ${RAW_VCF}"
    else
        echo "  Running Mutect2..."
        ${GATK} Mutect2 \
            -R "${REF}" \
            -I "${LOCAL_TUMOR}" \
            -I "${LOCAL_NORMAL}" \
            -normal "${NORMAL_SM_TAG}" \
            --panel-of-normals "${PON_VCF}" \
            -O "${RAW_VCF}" \
            --stats "${STATS}" \
            --native-pair-hmm-threads ${THREADS} 2>&1 \
        || { STATUS="FAIL_MUTECT2"; echo "  [ERROR] Mutect2 failed"; }
    fi

    # --- FilterMutectCalls ---------------------------------------------------
    if [[ -f "${FILT_VCF}" ]]; then
        echo "  [SKIP] Filtered VCF exists: ${FILT_VCF}"
    elif [[ -f "${RAW_VCF}" ]]; then
        echo "  Running FilterMutectCalls..."
        ${GATK} FilterMutectCalls \
            -R "${REF}" \
            -V "${RAW_VCF}" \
            --stats "${STATS}" \
            -O "${FILT_VCF}" 2>&1 \
        || { STATUS="FAIL_FILTER"; echo "  [ERROR] FilterMutectCalls failed"; }
    fi

    # --- Count variants ------------------------------------------------------
    TOTAL_VARS="NA"; PASS_VARS="NA"; SNV_PASS="NA"; INDEL_PASS="NA"
    if [[ -f "${FILT_VCF}" ]]; then
        TOTAL_VARS=$(${GATK} CountVariants -V "${RAW_VCF}"  2>/dev/null | tail -1 || echo "NA")
        PASS_VARS=$( ${GATK} CountVariants -V "${FILT_VCF}" --variant-filter-name PASS 2>/dev/null | tail -1 || echo "NA")
        # Count SNVs vs indels in PASS calls
        SNV_PASS=$( /home/ec2-user/miniforge3/bin/bcftools view -f PASS "${FILT_VCF}" 2>/dev/null \
            | /home/ec2-user/miniforge3/bin/bcftools stats 2>/dev/null \
            | awk '/^SN.*number of SNPs/{print $NF}' || echo "NA")
        INDEL_PASS=$(/home/ec2-user/miniforge3/bin/bcftools view -f PASS "${FILT_VCF}" 2>/dev/null \
            | /home/ec2-user/miniforge3/bin/bcftools stats 2>/dev/null \
            | awk '/^SN.*number of indels/{print $NF}' || echo "NA")
        echo "  Total variants : ${TOTAL_VARS}"
        echo "  PASS variants  : ${PASS_VARS}"
        echo "  SNVs (PASS)    : ${SNV_PASS}"
        echo "  Indels (PASS)  : ${INDEL_PASS}"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${TUMOR_SAMPLE}" "${NORMAL_BAM}" "${RAW_VCF}" "${FILT_VCF}" \
        "${TOTAL_VARS}" "${PASS_VARS}" "${SNV_PASS}" "${INDEL_PASS}" "${STATUS}" \
        >> "${RESULTS_TSV}"

    # Upload VCFs and stats to S3, then remove local copies to free disk
    echo "  Uploading results to S3..."
    for F in "${RAW_VCF}" "${RAW_VCF}.tbi" \
              "${FILT_VCF}" "${FILT_VCF}.tbi" \
              "${STATS}"; do
        [[ -f "${F}" ]] && aws s3 cp "${F}" "${S3_RESULTS}/$(basename $(dirname ${F}))/$(basename ${F})" \
            && rm -f "${F}"
    done

    # Clean up tumor BAM (keep normal for next iteration)
    rm -f "${LOCAL_TUMOR}" "${LOCAL_TUMOR}.bai"
    echo "  Done: ${TUMOR_SAMPLE} [${STATUS}]"
done

# Clean up normal BAM
rm -f "${LOCAL_NORMAL}" "${LOCAL_NORMAL}.bai"

echo ""
echo "========================================"
echo "[INFO] Variant calling complete."
echo "[INFO] Results: ${RESULTS_TSV}"
echo "========================================"
echo ""
cat "${RESULTS_TSV}"

# --- Run R summary -----------------------------------------------------------
echo ""
echo "[INFO] Generating variant summary report with R..."
Rscript /home/ec2-user/code/summarize_variants.R \
    "${RESULTS_TSV}" \
    "${OUTDIR}" 2>&1
