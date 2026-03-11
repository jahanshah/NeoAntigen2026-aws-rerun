#!/usr/bin/env bash
# =============================================================================
# aws_master_pipeline.sh
# =============================================================================
# Full neoantigen discovery pipeline — AWS Linux
#
# Blocks:
#   A. WES Preprocessing        (Picard, samtools)
#   B. Somatic Variant Calling  (GATK4 Mutect2, bcftools)
#   C. Variant Annotation       (SnpEff)
#   D. Clonal Evolution         (PyClone)
#   E. Neoantigen Extraction    (custom Python)
#   F. MHC Binding Prediction   (IEDB REST API)
#   G. RNA-seq Processing       (GATK SplitNCigarReads, ASEReadCounter)
#   H. Differential Expression  (DESeq2 via Rscript)
#   I. Fusion Analysis          (Arriba filter + scoring)
#   J. Integration & Report     (R + Python)
#
# Usage:
#   bash aws_master_pipeline.sh [BLOCK]
#   bash aws_master_pipeline.sh          # run all blocks
#   bash aws_master_pipeline.sh A        # run only block A
#   bash aws_master_pipeline.sh A B C    # run specific blocks
#
# AWS recommended instance:
#   Block A-D : m5.8xlarge  (32 vCPU, 128 GB RAM)
#   Block E-J : m5.4xlarge  (16 vCPU,  64 GB RAM)
#   Storage   : 2 TB gp3 EBS (or s3 mount via s3fs)
# =============================================================================

set -euo pipefail
IFS=$'\n\t'

# =============================================================================
# ── SECTION 0: CONFIGURATION — edit before running ───────────────────────────
# =============================================================================

# ── Project root (set to your AWS EBS mount or cloned repo) ──────────────────
WDIR="${WDIR:-/data/neoantigen}"
SCRIPTS="$WDIR/scripts"

# ── S3 bucket (optional — set empty string to skip S3 sync) ──────────────────
S3_BUCKET="${S3_BUCKET:-s3://your-bucket/neoantigen}"

# ── Reference files ───────────────────────────────────────────────────────────
GENOME="$WDIR/data/mm10.fa"
GENOME_DICT="$WDIR/data/mm10.dict"
GTF="$WDIR/reference_db/gencode.vM25.annotation.gtf.gz"
PROTEINS="$WDIR/reference_db/proteins.fa"
EXON_BED="$WDIR/reference_db/merged_exons.bed"
KNOWN_SITES="$WDIR/reference_db/mgp_REL2021_snps.vcf.gz"   # Mouse SNP DB for BQSR
GNOMAD_MM10="$WDIR/reference_db/gnomad.mm10.vcf.gz"        # for Mutect2 germline

# ── MHC allele ────────────────────────────────────────────────────────────────
ALLELE="H-2-Db"
PEP_LENGTHS="8,9,10,11"

# ── CPU / Memory ──────────────────────────────────────────────────────────────
THREADS=16
MEM="32g"

# ── Conda environments ────────────────────────────────────────────────────────
GATK_ENV="gatk-env"
PICARD_ENV="target-env"
SNPEFF_ENV="snpeff-env"
PYCLONE_ENV="pyclone-env"
R_ENV="r-env"

# ── Tool paths (resolved from conda envs) ────────────────────────────────────
PICARD="$(conda run -n $PICARD_ENV which picard 2>/dev/null || echo picard)"
SNPEFF="$(conda run -n $SNPEFF_ENV which snpEff 2>/dev/null || echo snpEff)"

# ── Normal sample (used as matched normal in Mutect2) ─────────────────────────
NORMAL_ID="424"
NORMAL_LABEL="D0_normal"
NORMAL_BAM="$WDIR/data/bam/424_marked_dup.bam"

# =============================================================================
# ── SAMPLE DEFINITION ─────────────────────────────────────────────────────────
# Format: "SAMPLE_ID|LABEL|WES_BAM|RNA_BAM|STAGE|BATCH"
# WES_BAM or RNA_BAM = "NA" if not available
# STAGE: Pre / Early / Trans / Peak / Late
# BATCH: old / new
# =============================================================================

declare -A WES_BAM RNA_BAM STAGE BATCH_GROUP

define_samples() {
    # WES + RNA
    WES_BAM[423]="$WDIR/data/bam/423_D0_old.bam";        RNA_BAM[423]="$WDIR/data/rnaseq_bam/D0_old.bam";       STAGE[423]="Pre";   BATCH_GROUP[423]="old"
    WES_BAM[428]="$WDIR/data/bam/428_D20_new.bam";       RNA_BAM[428]="$WDIR/data/rnaseq_bam/D20_new.bam";      STAGE[428]="Early"; BATCH_GROUP[428]="new"
    WES_BAM[443]="$WDIR/data/bam/443_D21_new.bam";       RNA_BAM[443]="$WDIR/data/rnaseq_bam/D21_new.bam";      STAGE[443]="Early"; BATCH_GROUP[443]="new"
    WES_BAM[34]="$WDIR/data/bam/34_D52_old.bam";         RNA_BAM[34]="$WDIR/data/rnaseq_bam/D52_old.bam";       STAGE[34]="Trans";  BATCH_GROUP[34]="old"
    WES_BAM[32]="$WDIR/data/bam/32_D99_old.bam";         RNA_BAM[32]="$WDIR/data/rnaseq_bam/D99_old.bam";       STAGE[32]="Peak";   BATCH_GROUP[32]="old"
    WES_BAM[36]="$WDIR/data/bam/36_D99_new.bam";         RNA_BAM[36]="$WDIR/data/rnaseq_bam/D99_36_new.bam";    STAGE[36]="Peak";   BATCH_GROUP[36]="new"
    WES_BAM[38]="$WDIR/data/bam/38_D99_new.bam";         RNA_BAM[38]="$WDIR/data/rnaseq_bam/D99_38_new.bam";    STAGE[38]="Peak";   BATCH_GROUP[38]="new"
    WES_BAM[42]="$WDIR/data/bam/42_D122_old.bam";        RNA_BAM[42]="$WDIR/data/rnaseq_bam/D122_old.bam";      STAGE[42]="Late";   BATCH_GROUP[42]="old"

    # RNA only (no WES BAM)
    WES_BAM[D88]="NA";                                   RNA_BAM[D88]="$WDIR/data/rnaseq_bam/D88_old.bam";      STAGE[D88]="Trans"; BATCH_GROUP[D88]="old"
    WES_BAM[D109]="NA";                                  RNA_BAM[D109]="$WDIR/data/rnaseq_bam/D109_new.bam";    STAGE[D109]="Late"; BATCH_GROUP[D109]="new"
}

# All WES tumor sample IDs (have BAM)
WES_SAMPLES=(423 428 443 34 32 36 38 42)
# All RNA samples
RNA_SAMPLES=(423 428 443 34 32 36 38 42 D88 D109)
# All samples with both WES + RNA
JOINT_SAMPLES=(423 428 443 34 32 36 38 42)

# =============================================================================
# ── HELPER FUNCTIONS ──────────────────────────────────────────────────────────
# =============================================================================

LOG_DIR="$WDIR/logs/pipeline"
DONE_DIR="$WDIR/.done"
mkdir -p "$LOG_DIR" "$DONE_DIR"

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_DIR/master.log"; }
done_flag() { echo "$DONE_DIR/$1"; }
is_done()   { [[ -f "$(done_flag "$1")" ]]; }
mark_done() { touch "$(done_flag "$1")"; log "DONE: $1"; }

skip_or_run() {
    # skip_or_run "flag_name" "description" cmd...
    local flag="$1" desc="$2"; shift 2
    if is_done "$flag"; then
        log "SKIP (already done): $desc"
    else
        log "START: $desc"
        "$@" && mark_done "$flag" || { log "FAILED: $desc"; exit 1; }
    fi
}

require_bam() {
    local bam="$1"
    if [[ "$bam" == "NA" || ! -f "$bam" ]]; then
        log "  SKIP — BAM not available: $bam"
        return 1
    fi
    return 0
}

s3_sync_out() {
    [[ -n "$S3_BUCKET" ]] && aws s3 sync "$1" "$S3_BUCKET/$2/" --quiet || true
}

# =============================================================================
# ── BLOCK A: WES PREPROCESSING ───────────────────────────────────────────────
# Per sample: AddReadGroups → Sort → MarkDuplicates → Index
# Tools: Picard, samtools
# AWS: m5.4xlarge (8 vCPU enough per sample)
# =============================================================================

block_A() {
    log "════════════════════════════════════════"
    log "BLOCK A — WES Preprocessing"
    log "════════════════════════════════════════"

    local out_dir="$WDIR/data/bam/processed"
    mkdir -p "$out_dir"

    for sid in "${WES_SAMPLES[@]}"; do
        local raw_bam="${WES_BAM[$sid]}"
        require_bam "$raw_bam" || continue

        local rg_bam="$out_dir/${sid}_rg.bam"
        local sorted_bam="$out_dir/${sid}_sorted.bam"
        local markdup_bam="$out_dir/${sid}_marked_dup.bam"
        local metrics="$out_dir/${sid}_dup_metrics.txt"
        local log_base="$LOG_DIR/A_${sid}"

        log "── Sample $sid ──"

        # A1. Add Read Groups
        skip_or_run "A1_rg_${sid}" "AddReadGroups $sid" \
            conda run -n $PICARD_ENV picard -Xmx${MEM} AddOrReplaceReadGroups \
                I="$raw_bam" O="$rg_bam" \
                RGID="$sid" RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM="$sid" \
                2>"${log_base}_rg.log"

        # A2. Sort
        skip_or_run "A2_sort_${sid}" "Sort $sid" \
            samtools sort -@ "$THREADS" -o "$sorted_bam" "$rg_bam" \
                2>"${log_base}_sort.log"
        rm -f "$rg_bam"

        # A3. Mark Duplicates
        skip_or_run "A3_markdup_${sid}" "MarkDuplicates $sid" \
            conda run -n $PICARD_ENV picard -Xmx${MEM} MarkDuplicates \
                I="$sorted_bam" O="$markdup_bam" M="$metrics" \
                USE_JDK_DEFLATER=true USE_JDK_INFLATER=true \
                2>"${log_base}_markdup.log"
        rm -f "$sorted_bam"

        # A4. Index
        skip_or_run "A4_index_${sid}" "Index $sid" \
            samtools index -@ "$THREADS" "$markdup_bam" \
                2>"${log_base}_index.log"

        # Update BAM path to processed
        WES_BAM[$sid]="$markdup_bam"
    done

    s3_sync_out "$out_dir" "bam/processed"
    log "Block A complete"
}

# =============================================================================
# ── BLOCK B: SOMATIC VARIANT CALLING ─────────────────────────────────────────
# Per sample: Mutect2 → FilterMutectCalls → bcftools exon filter
# Tools: GATK4, bcftools
# AWS: m5.8xlarge (Mutect2 is memory-heavy)
# =============================================================================

block_B() {
    log "════════════════════════════════════════"
    log "BLOCK B — Somatic Variant Calling"
    log "════════════════════════════════════════"

    local out_dir="$WDIR/results/vcf"
    mkdir -p "$out_dir"

    for sid in "${WES_SAMPLES[@]}"; do
        local tumor_bam="${WES_BAM[$sid]}"
        require_bam "$tumor_bam" || continue

        local log_base="$LOG_DIR/B_${sid}"
        local raw_vcf="$out_dir/${sid}_mutect2.vcf.gz"
        local filtered_vcf="$out_dir/${sid}_filtered.vcf.gz"
        local exon_vcf="$out_dir/${sid}_filtered_exons.vcf.gz"
        local f1r2="$out_dir/${sid}_f1r2.tar.gz"

        log "── Sample $sid ──"

        # B1. Mutect2 tumor-vs-normal
        skip_or_run "B1_mutect2_${sid}" "Mutect2 $sid" \
            conda run -n $GATK_ENV gatk Mutect2 \
                -R "$GENOME" \
                -I "$tumor_bam"  \
                -I "$NORMAL_BAM" \
                -normal "$NORMAL_ID" \
                --f1r2-tar-gz "$f1r2" \
                --germline-resource "$GNOMAD_MM10" \
                -O "$raw_vcf" \
                --native-pair-hmm-threads "$THREADS" \
                2>"${log_base}_mutect2.log"

        # B2. FilterMutectCalls
        skip_or_run "B2_filter_${sid}" "FilterMutectCalls $sid" \
            conda run -n $GATK_ENV gatk FilterMutectCalls \
                -R "$GENOME" \
                -V "$raw_vcf" \
                -O "$filtered_vcf" \
                2>"${log_base}_filter.log"

        # B3. bcftools: PASS + AF>5% + exon regions
        skip_or_run "B3_exon_${sid}" "Exon filter $sid" bash -c "
            bcftools view -f PASS '$filtered_vcf' \
                | bcftools filter -i 'FORMAT/AF[0:0]>0.05' \
                | bcftools view -R '$EXON_BED' -Oz -o '$exon_vcf' \
                2>'${log_base}_exon.log'
            bcftools index '$exon_vcf'
        "

        # B4. Variant stats
        skip_or_run "B4_stats_${sid}" "Variant stats $sid" bash -c "
            bcftools stats '$exon_vcf' > '${log_base}_vcf_stats.txt'
            echo 'Sample ${sid}: '$(bcftools view -H '$exon_vcf' | wc -l)' variants' >> '$LOG_DIR/variant_counts.txt'
        "
    done

    s3_sync_out "$out_dir" "results/vcf"
    log "Block B complete"
}

# =============================================================================
# ── BLOCK C: VARIANT ANNOTATION ──────────────────────────────────────────────
# Per sample: SnpEff → annotated VCF with ANN field
# Tools: SnpEff (GRCm38.86)
# =============================================================================

block_C() {
    log "════════════════════════════════════════"
    log "BLOCK C — Variant Annotation (SnpEff)"
    log "════════════════════════════════════════"

    local in_dir="$WDIR/results/vcf"
    local out_dir="$WDIR/results/annotated_vcf"
    mkdir -p "$out_dir"

    for sid in "${WES_SAMPLES[@]}"; do
        local exon_vcf="$in_dir/${sid}_filtered_exons.vcf.gz"
        [[ -f "$exon_vcf" ]] || { log "SKIP $sid — VCF not found"; continue; }

        local ann_vcf="$out_dir/${sid}_annotated.vcf"
        local stats_csv="$out_dir/${sid}_stats.csv"
        local log_file="$LOG_DIR/C_${sid}_snpeff.log"

        log "── Sample $sid ──"

        skip_or_run "C1_snpeff_${sid}" "SnpEff $sid" bash -c "
            bcftools view '$exon_vcf' \
                | conda run -n $SNPEFF_ENV snpEff -Xmx4g GRCm38.86 \
                    -csvStats '$stats_csv' \
                    -noLog \
                > '$ann_vcf' \
                2>'$log_file'
        "

        skip_or_run "C2_snpeff_count_${sid}" "SnpEff count $sid" bash -c "
            HIGH=\$(grep -c 'HIGH' '$ann_vcf' || true)
            MOD=\$(grep -c 'MODERATE' '$ann_vcf' || true)
            echo 'Sample $sid: HIGH='\$HIGH' MODERATE='\$MOD >> '$LOG_DIR/annotation_counts.txt'
        "
    done

    s3_sync_out "$out_dir" "results/annotated_vcf"
    log "Block C complete"
}

# =============================================================================
# ── BLOCK D: CLONAL EVOLUTION (PyClone) ──────────────────────────────────────
# Per sample: build YAML → Joint multi-sample PyClone run
# Tools: PyClone (conda env)
# Note: PyClone needs shared mutations across samples; run jointly
# =============================================================================

block_D() {
    log "════════════════════════════════════════"
    log "BLOCK D — Clonal Evolution (PyClone)"
    log "════════════════════════════════════════"

    local ann_dir="$WDIR/results/annotated_vcf"
    local vcf_dir="$WDIR/results/vcf"
    local yaml_dir="$WDIR/pyclone_run/yaml"
    local trace_dir="$WDIR/pyclone_run/trace"
    local table_dir="$WDIR/pyclone_run/tables"
    mkdir -p "$yaml_dir" "$trace_dir" "$table_dir"

    # D1. Build YAML per sample
    for sid in "${WES_SAMPLES[@]}"; do
        local ann_vcf="$ann_dir/${sid}_annotated.vcf"
        local exon_vcf="$vcf_dir/${sid}_filtered_exons.vcf.gz"
        local yaml_out="$yaml_dir/${sid}_annotated.yaml"
        [[ -f "$ann_vcf" ]] || continue

        skip_or_run "D1_yaml_${sid}" "PyClone YAML $sid" \
            python3 "$SCRIPTS/prepare_pyclone_input.py" \
                --vcf  "$exon_vcf" \
                --sample_id "$sid" \
                --out  "$yaml_out" \
                2>"$LOG_DIR/D_${sid}_yaml.log"
    done

    # D2. Build joint config
    skip_or_run "D2_pyclone_config" "Build PyClone config" bash -c "
        cat > '$WDIR/pyclone_run/config.yaml' << 'YAML'
num_iters: 10000
base_measure_params:
  alpha: 1.0
  beta: 1.0
concentration:
  value: 1.0
density: beta_binomial
beta_binomial_precision_params:
  value: 1000
  prior:
    rate: 0.001
    shape: 1.0
init_method: connected
working_dir: $WDIR/pyclone_run
trace_dir: $WDIR/pyclone_run/trace
samples:
$(for sid in \${WES_SAMPLES[@]}; do
    yaml='$yaml_dir'/${sid}_annotated.yaml
    [[ -f \$yaml ]] && printf '  %s:\n    mutations_file: %s\n    tumour_content: {value: 1.0}\n' \$sid \$yaml
done)
YAML
    "

    # D3. Run PyClone MCMC
    skip_or_run "D3_pyclone_run" "PyClone MCMC run" \
        conda run -n $PYCLONE_ENV PyClone run_analysis \
            --config_file "$WDIR/pyclone_run/config.yaml" \
            2>"$LOG_DIR/D_pyclone_run.log"

    # D4. Build output tables
    skip_or_run "D4_pyclone_tables" "PyClone build tables" bash -c "
        conda run -n $PYCLONE_ENV PyClone build_table \
            --config_file '$WDIR/pyclone_run/config.yaml' \
            --table_type cluster \
            --out_file '$table_dir/cluster.tsv' \
            2>>'$LOG_DIR/D_pyclone_run.log'

        conda run -n $PYCLONE_ENV PyClone build_table \
            --config_file '$WDIR/pyclone_run/config.yaml' \
            --table_type loci \
            --out_file '$table_dir/loci.tsv' \
            2>>'$LOG_DIR/D_pyclone_run.log'
    "

    s3_sync_out "$WDIR/pyclone_run" "pyclone_run"
    log "Block D complete"
}

# =============================================================================
# ── BLOCK E: NEOANTIGEN EXTRACTION ───────────────────────────────────────────
# Per sample:
#   E1. Missense/inframe/stop peptides  (extract_mutant_peptides.py)
#   E2. Frameshift junction peptides    (extract_frameshift_peptides.py)
#   E3. Fusion junction peptides        (from Arriba TSVs)
# =============================================================================

block_E() {
    log "════════════════════════════════════════"
    log "BLOCK E — Neoantigen Peptide Extraction"
    log "════════════════════════════════════════"

    local ann_dir="$WDIR/results/annotated_vcf"
    local pep_dir="$WDIR/peptides"
    local fs_dir="$WDIR/peptides_frameshift"
    local fus_dir="$WDIR/peptides_fusion"
    local arriba_dir="$WDIR/data/FusionCalling/arriba"
    mkdir -p "$pep_dir" "$fs_dir" "$fus_dir"

    for sid in "${WES_SAMPLES[@]}"; do
        local ann_vcf="$ann_dir/${sid}_annotated.vcf"
        [[ -f "$ann_vcf" ]] || continue
        log "── Sample $sid ──"

        # E1. Missense + inframe + stop peptides
        skip_or_run "E1_pep_${sid}" "Extract mutant peptides $sid" \
            python3 "$SCRIPTS/extract_mutant_peptides.py" \
                --vcf       "$ann_vcf" \
                --proteins  "$PROTEINS" \
                --lengths   "$PEP_LENGTHS" \
                --sample_id "$sid" \
                --out_fasta "$pep_dir/${sid}_peptides.fa" \
                --out_tsv   "$pep_dir/${sid}_peptides.tsv" \
                2>"$LOG_DIR/E_${sid}_peptides.log"

        # E2. Frameshift junction peptides
        skip_or_run "E2_fs_${sid}" "Extract frameshift peptides $sid" bash -c "
            HIGH=\$(grep -c 'frameshift_variant.*HIGH' '$ann_vcf' 2>/dev/null || echo 0)
            if [[ \$HIGH -gt 0 ]]; then
                python3 '$SCRIPTS/extract_frameshift_peptides.py' \
                    --vcf       '$ann_vcf' \
                    --genome    '$GENOME' \
                    --gtf       '$GTF' \
                    --proteins  '$PROTEINS' \
                    --lengths   '$PEP_LENGTHS' \
                    --sample_id '$sid' \
                    --out_fasta '$fs_dir/${sid}_fs_peptides.fa' \
                    --out_tsv   '$fs_dir/${sid}_fs_peptides.tsv' \
                    2>'$LOG_DIR/E_${sid}_frameshift.log'
            else
                echo 'No HIGH frameshift variants — skipping'
            fi
        "
    done

    # E3. Fusion junction peptides (from Arriba — all samples at once)
    skip_or_run "E3_fusion_peptides" "Extract fusion junction peptides" \
        python3 - <<'PYEOF'
import os, csv, re
from itertools import product

WDIR    = os.environ.get("WDIR", "/data/neoantigen")
in_dir  = os.path.join(WDIR, "data/FusionCalling/arriba")
out_tsv = os.path.join(WDIR, "peptides_fusion/fusion_peptides.tsv")
out_fa  = os.path.join(WDIR, "peptides_fusion/fusion_peptides.fa")

# Genes to exclude: Ig V(D)J recombination (normal B-cell biology)
IG_PAT = re.compile(r'^Ig[khl][vjc]', re.IGNORECASE)

LENGTHS = [8, 9, 10, 11]

def slide_windows(seq, lengths):
    for L in lengths:
        for i in range(len(seq) - L + 1):
            yield L, i, seq[i:i+L]

rows, seen = [], set()
with open(out_fa, "w") as fa_out, open(out_tsv, "w", newline="") as tsv_out:
    writer = csv.writer(tsv_out, delimiter="\t")
    writer.writerow(["Sample","Gene1","Gene2","Confidence","ReadFrame",
                     "Type","SplitReads","Breakpoint1","Breakpoint2",
                     "PepLength","PepSeq","JunctionPos"])

    for fname in sorted(os.listdir(in_dir)):
        if not fname.startswith("Arriba_") or not fname.endswith(".txt"):
            continue
        sample = fname.replace("Arriba_","").replace(".txt","")
        with open(os.path.join(in_dir, fname)) as f:
            reader = csv.DictReader(f, delimiter="\t")
            reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
            for row in reader:
                g1, g2   = row["gene1"], row["gene2"]
                conf     = row.get("confidence","")
                rf       = row.get("reading_frame","")
                pep_full = row.get("peptide_sequence","").strip()
                ftype    = row.get("type","")
                sr       = row.get("split_reads1","0")+"/"+row.get("split_reads2","0")
                bp1, bp2 = row.get("breakpoint1",""), row.get("breakpoint2","")

                # Filter: skip Ig V(D)J, require in-frame, require peptide
                if IG_PAT.match(g1) or IG_PAT.match(g2): continue
                if rf != "in-frame": continue
                if not pep_full or pep_full == ".": continue

                # Arriba peptide_sequence: upstream|downstream (| = junction)
                parts = pep_full.split("|")
                if len(parts) == 2:
                    junction_pos = len(parts[0])
                    full_seq = parts[0] + parts[1]
                else:
                    junction_pos = len(pep_full) // 2
                    full_seq = pep_full

                for L, start, pep in slide_windows(full_seq, LENGTHS):
                    # Only keep windows spanning the junction
                    end = start + L
                    if not (start < junction_pos <= end): continue
                    if pep in seen: continue
                    seen.add(pep)

                    hdr = f"FUSION|{sample}|{g1}--{g2}|{conf}|len{L}"
                    fa_out.write(f">{hdr}\n{pep}\n")
                    writer.writerow([sample,g1,g2,conf,rf,ftype,sr,
                                     bp1,bp2,L,pep,junction_pos])

print(f"Written {len(seen)} unique fusion junction peptides")
PYEOF

    s3_sync_out "$pep_dir" "peptides"
    s3_sync_out "$fs_dir"  "peptides_frameshift"
    s3_sync_out "$fus_dir" "peptides_fusion"
    log "Block E complete"
}

# =============================================================================
# ── BLOCK F: MHC BINDING PREDICTION (IEDB) ───────────────────────────────────
# Per sample per peptide type: FASTA → IEDB REST → scored TSV
# Tools: IEDB REST API (netMHCpan_EL)
# Note: API has rate limits — retry logic built in
# =============================================================================

block_F() {
    log "════════════════════════════════════════"
    log "BLOCK F — MHC Binding Prediction (IEDB)"
    log "════════════════════════════════════════"

    local pep_dir="$WDIR/peptides"
    local fs_dir="$WDIR/peptides_frameshift"
    local fus_dir="$WDIR/peptides_fusion"
    local score_dir="$WDIR/iedb_scores"
    mkdir -p "$score_dir"

    # F1. Missense/inframe/stop peptides — one TSV per sample
    for sid in "${WES_SAMPLES[@]}"; do
        local in_tsv="$pep_dir/${sid}_peptides.tsv"
        local out_tsv="$score_dir/${sid}_iedb_scores.tsv"
        [[ -f "$in_tsv" ]] || continue

        skip_or_run "F1_iedb_${sid}" "IEDB score missense/inframe $sid" \
            python3 "$SCRIPTS/score_indel_stop_iedb.py" \
                --in_tsv  "$in_tsv" \
                --allele  "$ALLELE" \
                --lengths "$PEP_LENGTHS" \
                --out_tsv "$out_tsv" \
                2>"$LOG_DIR/F_${sid}_iedb.log"
    done

    # F2. Frameshift peptides — one TSV per sample
    for sid in "${WES_SAMPLES[@]}"; do
        local in_tsv="$fs_dir/${sid}_fs_peptides.tsv"
        local out_tsv="$score_dir/${sid}_fs_iedb_scores.tsv"
        [[ -f "$in_tsv" ]] || continue

        skip_or_run "F2_fs_iedb_${sid}" "IEDB score frameshift $sid" \
            python3 "$SCRIPTS/score_frameshift_peptides.py" \
                --in_tsv  "$in_tsv" \
                --sample  "$sid" \
                --allele  "$ALLELE" \
                --lengths "$PEP_LENGTHS" \
                --out_tsv "$out_tsv" \
                2>"$LOG_DIR/F_${sid}_fs_iedb.log"
    done

    # F3. Fusion peptides — single combined file
    local fus_tsv="$fus_dir/fusion_peptides.tsv"
    local fus_out="$score_dir/fusion_iedb_scores.tsv"
    if [[ -f "$fus_tsv" ]]; then
        skip_or_run "F3_fusion_iedb" "IEDB score fusion peptides" \
            python3 "$SCRIPTS/score_indel_stop_iedb.py" \
                --in_tsv  "$fus_tsv" \
                --allele  "$ALLELE" \
                --lengths "$PEP_LENGTHS" \
                --out_tsv "$fus_out" \
                2>"$LOG_DIR/F_fusion_iedb.log"
    fi

    s3_sync_out "$score_dir" "iedb_scores"
    log "Block F complete"
}

# =============================================================================
# ── BLOCK G: RNA-seq PROCESSING ──────────────────────────────────────────────
# Per sample:
#   G1. SplitNCigarReads (STAR BAM → GATK-ready BAM)
#   G2. ASEReadCounter   (mutant allele expression at WES variant sites)
# Tools: GATK4
# =============================================================================

block_G() {
    log "════════════════════════════════════════"
    log "BLOCK G — RNA-seq Processing"
    log "════════════════════════════════════════"

    local vcf_dir="$WDIR/results/vcf"
    local rna_proc_dir="$WDIR/data/rnaseq_bam/processed"
    local ase_dir="$WDIR/results/ase"
    mkdir -p "$rna_proc_dir" "$ase_dir"

    for sid in "${JOINT_SAMPLES[@]}"; do
        local rna_bam="${RNA_BAM[$sid]}"
        require_bam "$rna_bam" || continue

        local split_bam="$rna_proc_dir/${sid}_splitN.bam"
        local exon_vcf="$vcf_dir/${sid}_filtered_exons.vcf.gz"
        local ase_out="$ase_dir/${sid}_ase.tsv"
        local log_base="$LOG_DIR/G_${sid}"

        log "── Sample $sid ──"

        # G1. SplitNCigarReads — split RNA splice junctions for GATK
        skip_or_run "G1_splitN_${sid}" "SplitNCigarReads $sid" \
            conda run -n $GATK_ENV gatk SplitNCigarReads \
                -R "$GENOME" \
                -I "$rna_bam" \
                -O "$split_bam" \
                2>"${log_base}_splitN.log"

        # G2. ASEReadCounter — count ref/alt reads at WES variant sites
        if [[ -f "$exon_vcf" ]]; then
            skip_or_run "G2_ase_${sid}" "ASEReadCounter $sid" \
                conda run -n $GATK_ENV gatk ASEReadCounter \
                    -R "$GENOME" \
                    -I "$split_bam" \
                    -V "$exon_vcf" \
                    -O "$ase_out" \
                    --min-mapping-quality 10 \
                    --min-base-quality 20 \
                    2>"${log_base}_ase.log"
        else
            log "  SKIP ASEReadCounter $sid — no WES VCF available"
        fi
    done

    # G3. RNA-only samples — expression check only (no WES VCF)
    for sid in D88 D109; do
        local rna_bam="${RNA_BAM[$sid]}"
        require_bam "$rna_bam" || continue
        local split_bam="$rna_proc_dir/${sid}_splitN.bam"

        skip_or_run "G1_splitN_${sid}" "SplitNCigarReads $sid" \
            conda run -n $GATK_ENV gatk SplitNCigarReads \
                -R "$GENOME" \
                -I "$rna_bam" \
                -O "$split_bam" \
                2>"$LOG_DIR/G_${sid}_splitN.log"
    done

    s3_sync_out "$ase_dir" "results/ase"
    log "Block G complete"
}

# =============================================================================
# ── BLOCK H: DIFFERENTIAL EXPRESSION & EXPRESSION VALIDATION ─────────────────
# H1. DESeq2 DE analysis (Pre+Early vs Peak, Peak vs Late)
# H2. Expression validation per neoantigen source gene
# H3. ASE integration — flag neoantigens with RNA mutant allele evidence
# Tools: DESeq2, R
# =============================================================================

block_H() {
    log "════════════════════════════════════════"
    log "BLOCK H — Differential Expression & RNA Validation"
    log "════════════════════════════════════════"

    local counts_dir="$WDIR/data/RNASeq/count_files"
    local ase_dir="$WDIR/results/ase"
    local neo_dir="$WDIR/neoantigen_results"
    local fig_dir="$WDIR/figures_R"
    mkdir -p "$neo_dir" "$fig_dir"

    # H1. QC + batch correction + DE (existing R script)
    skip_or_run "H1_rnaseq_qc" "RNA-seq QC and DE" \
        conda run -n $R_ENV Rscript "$SCRIPTS/rnaseq_qc_exploration.R" \
            2>"$LOG_DIR/H_rnaseq_qc.log"

    skip_or_run "H2_rnaseq_de" "RNA-seq WES integration + DE" \
        conda run -n $R_ENV Rscript "$SCRIPTS/rnaseq_wes_integration.R" \
            2>"$LOG_DIR/H_rnaseq_de.log"

    # H3. Expression validation per neoantigen
    skip_or_run "H3_expr_validation" "Expression validation of neoantigens" \
        python3 - <<'PYEOF'
import os, csv
import pandas as pd

WDIR = os.environ.get("WDIR", "/data/neoantigen")
counts_f = os.path.join(WDIR, "data/RNASeq/count_files",
               "Longitudinal GEM-Mouse Project (RNA)_normalizedCountsWithAnnotations.txt")
neo_f    = os.path.join(WDIR, "neoantigen_results/all_neoantigens_presence_absence.tsv")
out_f    = os.path.join(WDIR, "neoantigen_results/neoantigens_with_expression.tsv")

counts = pd.read_csv(counts_f, sep="\t", index_col=0)
neo    = pd.read_csv(neo_f, sep="\t")

# Map gene symbol to expression (max across peak samples)
PEAK_COLS = [c for c in counts.columns if "D99" in c or "Peak" in c.title()]

def get_expr(gene):
    match = counts[counts.index.str.contains(gene, na=False, regex=False)]
    if match.empty: return None
    return float(match[PEAK_COLS].max(axis=1).iloc[0]) if PEAK_COLS else None

neo["peak_expression"] = neo.get("gene_id_approx", neo.get("Protein_ID","")).apply(
    lambda g: get_expr(str(g)) if pd.notna(g) else None)
neo["expressed_at_peak"] = neo["peak_expression"].apply(
    lambda x: True if (x is not None and x > 10) else False)

neo.to_csv(out_f, sep="\t", index=False)
print(f"Written {len(neo)} rows with expression annotation → {out_f}")
PYEOF

    # H4. ASE integration — annotate neoantigens with RNA mutant allele evidence
    skip_or_run "H4_ase_integration" "ASE integration" \
        python3 - <<'PYEOF'
import os, glob, csv
import pandas as pd

WDIR   = os.environ.get("WDIR", "/data/neoantigen")
ase_dir = os.path.join(WDIR, "results/ase")
neo_f   = os.path.join(WDIR, "neoantigen_results/neoantigens_with_expression.tsv")
out_f   = os.path.join(WDIR, "neoantigen_results/neoantigens_final.tsv")

neo = pd.read_csv(neo_f, sep="\t")
ase_frames = []
for f in glob.glob(os.path.join(ase_dir, "*_ase.tsv")):
    sid = os.path.basename(f).replace("_ase.tsv","")
    df  = pd.read_csv(f, sep="\t", comment="#")
    df["sample"] = sid
    ase_frames.append(df)

if ase_frames:
    ase = pd.concat(ase_frames)
    # Flag positions with alt allele fraction > 0.05 in RNA
    ase["rna_mutant_expressed"] = ase["altCount"] / (ase["totalCount"] + 1) > 0.05
    ase_pos = set(ase[ase["rna_mutant_expressed"]][["position"]].squeeze())
    neo["rna_mutant_evidence"] = neo.get("POS", neo.get("Variant_Position","")).isin(ase_pos)
else:
    neo["rna_mutant_evidence"] = False

neo.to_csv(out_f, sep="\t", index=False)
print(f"Written {len(neo)} neoantigens with ASE annotation → {out_f}")
PYEOF

    s3_sync_out "$neo_dir" "neoantigen_results"
    log "Block H complete"
}

# =============================================================================
# ── BLOCK I: FUSION ANALYSIS ──────────────────────────────────────────────────
# I1. Filter Arriba output (remove Ig V(D)J, apply confidence filters)
# I2. Longitudinal fusion presence/absence matrix
# I3. Merge fusion neoantigens into final neoantigen table
# =============================================================================

block_I() {
    log "════════════════════════════════════════"
    log "BLOCK I — Fusion Analysis"
    log "════════════════════════════════════════"

    local arriba_dir="$WDIR/data/FusionCalling/arriba"
    local fus_dir="$WDIR/results/fusions"
    local neo_dir="$WDIR/neoantigen_results"
    mkdir -p "$fus_dir"

    # I1. Filter + classify all Arriba calls
    skip_or_run "I1_fusion_filter" "Filter Arriba fusions" \
        python3 - <<'PYEOF'
import os, csv, re
from collections import defaultdict
import pandas as pd

WDIR      = os.environ.get("WDIR", "/data/neoantigen")
in_dir    = os.path.join(WDIR, "data/FusionCalling/arriba")
out_f     = os.path.join(WDIR, "results/fusions/filtered_fusions.tsv")
recur_f   = os.path.join(WDIR, "results/fusions/recurrent_fusions.tsv")

IG_PAT    = re.compile(r'^Ig[khl][vjc]', re.IGNORECASE)
ONCO_GENES = {"Trp53","Myc","Ras","Nras","Kras","Hras","Braf","Akt1",
               "Nsd3","Kat6a","Nsd2","Meis1","Runx1","Etv6","Bcr","Abl1"}

all_rows = []
for fname in sorted(os.listdir(in_dir)):
    if not fname.startswith("Arriba_") or not fname.endswith(".txt"): continue
    sample = fname.replace("Arriba_","").replace(".txt","")
    with open(os.path.join(in_dir, fname)) as f:
        reader = csv.DictReader(f, delimiter="\t")
        reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
        for row in reader:
            g1, g2 = row["gene1"], row["gene2"]
            if IG_PAT.match(g1) or IG_PAT.match(g2): continue
            row["sample"]    = sample
            row["is_oncogenic"] = "yes" if (g1 in ONCO_GENES or g2 in ONCO_GENES) else "no"
            row["pair"] = f"{g1}--{g2}"
            all_rows.append(row)

df = pd.DataFrame(all_rows)
df.to_csv(out_f, sep="\t", index=False)

# Recurrent fusions (≥2 samples)
recur = df.groupby("pair").agg(
    n_samples=("sample","nunique"),
    samples=("sample", lambda x: "|".join(sorted(set(x)))),
    best_confidence=("confidence", lambda x: "high" if "high" in x.values else ("medium" if "medium" in x.values else "low")),
    any_inframe=("reading_frame", lambda x: "yes" if "in-frame" in x.values else "no"),
    is_oncogenic=("is_oncogenic", lambda x: "yes" if "yes" in x.values else "no"),
).reset_index()
recur = recur[recur["n_samples"] >= 2].sort_values("n_samples", ascending=False)
recur.to_csv(recur_f, sep="\t", index=False)

print(f"Filtered: {len(df)} non-Ig fusions | {len(recur)} recurrent pairs")
PYEOF

    # I2. Merge fusion neoantigens into final table
    skip_or_run "I2_fusion_merge" "Merge fusion neoantigens" \
        python3 "$SCRIPTS/merge_new_binders.py" \
            --scores  "$WDIR/iedb_scores/fusion_iedb_scores.tsv" \
            --meta    "$WDIR/peptides_fusion/fusion_peptides.tsv" \
            --current "$neo_dir/neoantigens_final.tsv" \
            --out     "$neo_dir/neoantigens_final.tsv" \
            2>"$LOG_DIR/I_fusion_merge.log" || \
        log "  Note: merge_new_binders.py may need --mutation_type fusion arg — check manually"

    s3_sync_out "$fus_dir"  "results/fusions"
    s3_sync_out "$neo_dir"  "neoantigen_results"
    log "Block I complete"
}

# =============================================================================
# ── BLOCK J: INTEGRATION & FINAL REPORT ──────────────────────────────────────
# J1. Build presence/absence matrix (filter_and_merge.py)
# J2. Priority scoring (R)
# J3. Final figures (R)
# J4. HTML report (build_html.py)
# J5. Sync all outputs to S3
# =============================================================================

block_J() {
    log "════════════════════════════════════════"
    log "BLOCK J — Integration & Final Report"
    log "════════════════════════════════════════"

    local score_dir="$WDIR/iedb_scores"
    local pep_dir="$WDIR/peptides"
    local pyclone_dir="$WDIR/pyclone_run/tables"
    local neo_dir="$WDIR/neoantigen_results"
    local fig_dir="$WDIR/figures_R"
    local report_dir="$WDIR/reports"
    mkdir -p "$neo_dir" "$fig_dir" "$report_dir"

    # J1. Merge all scored peptides into unified presence/absence matrix
    skip_or_run "J1_merge" "filter_and_merge — build presence/absence matrix" \
        python3 "$SCRIPTS/filter_and_merge.py" \
            --netmhc_dir  "$score_dir" \
            --peptide_dir "$pep_dir" \
            --loci_tsv    "$pyclone_dir/loci.tsv" \
            --cluster_tsv "$pyclone_dir/cluster.tsv" \
            --thr_sb 0.5 --thr_wb 2.0 --thr_pb 10.0 \
            --out_dir     "$neo_dir" \
            2>"$LOG_DIR/J_merge.log"

    # J2. Add mutation type annotation
    skip_or_run "J2_mut_type" "Add mutation type labels" \
        python3 "$SCRIPTS/add_mutation_type.py" \
            2>"$LOG_DIR/J_muttype.log" || true

    # J3. Neoantigen figures (per-type analysis)
    skip_or_run "J3_figures_per_type" "Plot per-type figures" \
        conda run -n $R_ENV Rscript "$SCRIPTS/plot_per_type_figures.R" \
            2>"$LOG_DIR/J_figures_per_type.log"

    skip_or_run "J3_tables_per_type" "Generate per-type tables" \
        conda run -n $R_ENV Rscript "$SCRIPTS/generate_per_type_tables.R" \
            2>"$LOG_DIR/J_tables_per_type.log"

    skip_or_run "J3_final_analysis" "Final integrated analysis" \
        conda run -n $R_ENV Rscript "$SCRIPTS/final_analysis.R" \
            2>"$LOG_DIR/J_final_analysis.log"

    # J4. Priority scoring
    skip_or_run "J3_priority" "Neoantigen priority scoring" \
        conda run -n $R_ENV Rscript "$SCRIPTS/rnaseq_wes_integration.R" \
            2>"$LOG_DIR/J_priority.log"

    # J5. HTML report
    skip_or_run "J5_report" "Build HTML report" bash -c "
        DEV=\$(ls '$WDIR'/dev*/scripts/build_html.py 2>/dev/null | sort -V | tail -1)
        if [[ -n \"\$DEV\" ]]; then
            python3 \"\$DEV\" && \
            cp \"\$(dirname \$DEV)/../reports/neoantigen_analysis.html\" \
               '$report_dir/neoantigen_analysis.html'
        else
            python3 '$SCRIPTS/build_final_report_html.py'
        fi
        2>'$LOG_DIR/J_report.log'
    "

    # J6. Final S3 sync — everything
    if [[ -n "$S3_BUCKET" ]]; then
        log "Syncing all results to $S3_BUCKET ..."
        aws s3 sync "$WDIR/neoantigen_results" "$S3_BUCKET/neoantigen_results/" --quiet
        aws s3 sync "$WDIR/figures_R"          "$S3_BUCKET/figures/"            --quiet
        aws s3 sync "$WDIR/reports"            "$S3_BUCKET/reports/"            --quiet
        aws s3 sync "$WDIR/results"            "$S3_BUCKET/results/"            --quiet
        aws s3 sync "$WDIR/logs"               "$S3_BUCKET/logs/"               --quiet
    fi

    log "Block J complete"
    log ""
    log "══════════════════════════════════════════════════"
    log " Pipeline complete — $(date)"
    log " Report: $report_dir/neoantigen_analysis.html"
    log "══════════════════════════════════════════════════"
}

# =============================================================================
# ── ENVIRONMENT SETUP (run once on a fresh AWS instance) ─────────────────────
# =============================================================================

setup_environment() {
    log "Setting up conda environments and tools ..."

    # Base tools
    conda install -y -c bioconda -c conda-forge \
        samtools bcftools gatk4 picard snpeff \
        star arriba pysam requests pandas 2>"$LOG_DIR/setup_conda.log"

    # R environment with Bioconductor
    conda create -y -n r-env -c conda-forge -c bioconda \
        r-base r-ggplot2 r-dplyr r-tidyr r-pheatmap \
        bioconductor-deseq2 bioconductor-limma \
        r-ggrepel r-patchwork r-cowplot \
        2>>"$LOG_DIR/setup_conda.log"

    # PyClone
    conda create -y -n pyclone-env -c bioconda -c conda-forge python=2.7 \
        pyclone 2>>"$LOG_DIR/setup_conda.log"

    # Download reference files
    bash "$SCRIPTS/download_references.sh" 2>"$LOG_DIR/setup_refs.log"

    log "Setup complete"
}

# =============================================================================
# ── MAIN ──────────────────────────────────────────────────────────────────────
# =============================================================================

main() {
    log "══════════════════════════════════════════════════"
    log " NeoAntigen Master Pipeline — AWS Linux"
    log " Start: $(date)"
    log " WDIR:  $WDIR"
    log "══════════════════════════════════════════════════"

    define_samples

    # Determine which blocks to run
    local BLOCKS=("${@:-A B C D E F G H I J}")
    if [[ $# -eq 0 ]]; then
        BLOCKS=(A B C D E F G H I J)
    else
        BLOCKS=("$@")
    fi

    for block in "${BLOCKS[@]}"; do
        case "$block" in
            setup) setup_environment ;;
            A)     block_A ;;
            B)     block_B ;;
            C)     block_C ;;
            D)     block_D ;;
            E)     block_E ;;
            F)     block_F ;;
            G)     block_G ;;
            H)     block_H ;;
            I)     block_I ;;
            J)     block_J ;;
            *)     log "Unknown block: $block"; exit 1 ;;
        esac
    done
}

main "$@"
