#!/bin/bash
# run_443_pipeline.sh
# Complete the sample 443 (D21_new) processing pipeline:
#   AddReadGroups → MarkDuplicates → Mutect2 → Filter → AnnotateVCF → PyClone YAML
# Run AFTER 443_sorted.bam is fully written (PID 73899 exits).

set -euo pipefail

WDIR="/Users/jahanshah/Documents/Consultant-NGS/Ahmed/Projects2026/NeoAntigen"
PICARD="/opt/anaconda3/envs/target-env/bin/picard"
SNPEFF="/opt/anaconda3/envs/snpeff-env/bin/snpEff"
SNPEFF_DB="GRCm38.86"
EXON_BED="$WDIR/reference_db/merged_exons.bed"
REF="$WDIR/data/mm10.fa"
NORMAL_BAM="$WDIR/bam/424_marked_dup.bam"
SID="443"
MEM="8G"

BAM_DIR="$WDIR/bam"
LOG_DIR="$WDIR/logs"
MUTECT_DIR="$WDIR/results/02_mutect2"
FILT_DIR="$WDIR/results/03_filtered_vcf"
ANNO_DIR="$WDIR/annotated_vcf"
YAML_DIR="$WDIR/pyclone_run/yaml"
RESULTS_YAML="$WDIR/results/05_pyclone_yaml"
SCRIPTS="$WDIR/scripts"

mkdir -p "$LOG_DIR" "$MUTECT_DIR" "$FILT_DIR" "$ANNO_DIR" "$YAML_DIR" "$RESULTS_YAML"

echo "=== Starting 443 pipeline at $(date) ==="

# ── Step 0: Free source BAM (no longer needed once sorted) ─────────────────
SOURCE_BAM="$WDIR/data/${SID}_D21_new.bam"
if [[ -f "$SOURCE_BAM" ]]; then
    echo "[0] Removing source BAM $SOURCE_BAM ($(du -h "$SOURCE_BAM" | cut -f1))"
    rm -f "$SOURCE_BAM"
    echo "[0] Freed source BAM"
fi

# ── Step 1: AddOrReplaceReadGroups ──────────────────────────────────────────
SORTED_BAM="$BAM_DIR/${SID}_sorted.bam"
RG_BAM="$BAM_DIR/${SID}_rg.bam"

if [[ ! -f "$RG_BAM" ]]; then
    echo "[1] AddOrReplaceReadGroups: $SORTED_BAM → $RG_BAM"
    "$PICARD" "-Xmx${MEM}" AddOrReplaceReadGroups \
        I="$SORTED_BAM" O="$RG_BAM" \
        RGID="$SID" RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="$SID" \
        VALIDATION_STRINGENCY=LENIENT \
        2>"$LOG_DIR/${SID}_rg.log"
    echo "[1] Done: $RG_BAM ($(du -h "$RG_BAM" | cut -f1))"
else
    echo "[1] Skipping AddReadGroups ($RG_BAM exists)"
fi

# Delete sorted BAM to free space
if [[ -f "$SORTED_BAM" ]]; then
    echo "[1] Freeing sorted BAM..."
    rm -f "$SORTED_BAM"
fi

# ── Step 2: MarkDuplicates ───────────────────────────────────────────────────
DEDUP_BAM="$BAM_DIR/${SID}_marked_dup.bam"

if [[ ! -f "$DEDUP_BAM" ]]; then
    echo "[2] MarkDuplicates: $RG_BAM → $DEDUP_BAM"
    "$PICARD" "-Xmx${MEM}" MarkDuplicates \
        I="$RG_BAM" O="$DEDUP_BAM" \
        M="$BAM_DIR/${SID}_marked_dup_metrics.txt" \
        VALIDATION_STRINGENCY=LENIENT \
        USE_JDK_DEFLATER=true USE_JDK_INFLATER=true \
        2>"$LOG_DIR/${SID}_dedup.log"
    samtools index "$DEDUP_BAM"
    echo "[2] Done: $DEDUP_BAM ($(du -h "$DEDUP_BAM" | cut -f1))"
else
    echo "[2] Skipping MarkDuplicates ($DEDUP_BAM exists)"
fi

# Delete rg BAM to free space
if [[ -f "$RG_BAM" ]]; then
    echo "[2] Freeing rg BAM..."
    rm -f "$RG_BAM"
fi

# ── Step 3: Mutect2 ──────────────────────────────────────────────────────────
MUTECT_VCF="$MUTECT_DIR/${SID}_mutect2.vcf.gz"

if [[ ! -f "$MUTECT_VCF" ]]; then
    echo "[3] Mutect2: ${SID} vs normal 424"
    gatk Mutect2 \
        -R "$REF" \
        -I "$DEDUP_BAM" -tumor "$SID" \
        -I "$NORMAL_BAM" -normal 424 \
        -O "$MUTECT_VCF" \
        --f1r2-tar-gz "$MUTECT_DIR/${SID}_f1r2.tar.gz" \
        2>"$LOG_DIR/${SID}_mutect2.log"
    echo "[3] Done: $MUTECT_VCF"
else
    echo "[3] Skipping Mutect2 ($MUTECT_VCF exists)"
fi

# ── Step 4: FilterMutectCalls ────────────────────────────────────────────────
FILTERED_VCF="$MUTECT_DIR/${SID}_filtered.vcf.gz"

if [[ ! -f "$FILTERED_VCF" ]]; then
    echo "[4] FilterMutectCalls"
    gatk FilterMutectCalls \
        -R "$REF" \
        -V "$MUTECT_VCF" \
        -O "$FILTERED_VCF" \
        2>"$LOG_DIR/${SID}_filter.log"
    echo "[4] Done: $FILTERED_VCF"
else
    echo "[4] Skipping FilterMutectCalls ($FILTERED_VCF exists)"
fi

# Delete dedup BAM to free space
if [[ -f "$DEDUP_BAM" ]]; then
    echo "[4] Freeing marked_dup BAM..."
    rm -f "$DEDUP_BAM" "${DEDUP_BAM}.bai"
fi

# ── Step 5: Filter PASS + exons + AF/DP thresholds ──────────────────────────
EXON_VCF="$FILT_DIR/${SID}_filtered_exons.vcf.gz"

if [[ ! -f "$EXON_VCF" ]]; then
    echo "[5] Filtering to PASS/exon/AF>0.05/DP>10"
    bcftools view -f PASS -T "$EXON_BED" "$FILTERED_VCF" | \
        bcftools filter -i 'FORMAT/AF[0:0]>0.05 && FORMAT/DP[0]>10' | \
        bcftools view -O z -o "$EXON_VCF" \
        2>"$LOG_DIR/${SID}_exon_filter.log"
    bcftools index "$EXON_VCF"
    n_var=$(bcftools stats "$EXON_VCF" 2>/dev/null | grep "^SN.*number of records" | awk '{print $NF}')
    echo "[5] Done: $EXON_VCF ($n_var variants)"
else
    echo "[5] Skipping exon filter ($EXON_VCF exists)"
fi

# ── Step 6: SnpEff annotation ────────────────────────────────────────────────
ANNO_VCF="$ANNO_DIR/${SID}_annotated.vcf"

if [[ ! -f "$ANNO_VCF" ]]; then
    echo "[6] SnpEff annotation"
    "$SNPEFF" ann -Xmx4g -v "$SNPEFF_DB" \
        -csvStats "$ANNO_DIR/${SID}_stats.csv" \
        "$EXON_VCF" > "$ANNO_VCF" \
        2>"$LOG_DIR/${SID}_snpeff.log"
    echo "[6] Done: $ANNO_VCF ($(wc -l < "$ANNO_VCF") lines)"
else
    echo "[6] Skipping SnpEff ($ANNO_VCF exists)"
fi

# ── Step 7: PyClone YAML ──────────────────────────────────────────────────────
YAML_OUT="$YAML_DIR/${SID}_annotated.yaml"

if [[ ! -f "$YAML_OUT" ]]; then
    echo "[7] Preparing PyClone YAML"
    python3 "$SCRIPTS/prepare_pyclone_input.py" \
        --vcf       "$ANNO_VCF" \
        --sample_id "$SID" \
        --out       "$YAML_OUT" \
        2>"$LOG_DIR/${SID}_pyclone_prep.log"
    cp "$YAML_OUT" "$RESULTS_YAML/${SID}_annotated.yaml"
    echo "[7] Done: $YAML_OUT"
else
    echo "[7] Skipping PyClone YAML ($YAML_OUT exists)"
fi

echo ""
echo "=== 443 pipeline complete at $(date) ==="
echo "Disk usage:"
df -h "$WDIR"
