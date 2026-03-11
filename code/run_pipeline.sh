#!/usr/bin/env bash
# =============================================================================
# NeoAntigen Pipeline — Main Orchestrator
# Runs all 10 steps in order; each step is skipped if its S3 output exists.
#
# Usage:
#   ./run_pipeline.sh [START_STEP] [END_STEP]
#   ./run_pipeline.sh            # run all steps 1-10
#   ./run_pipeline.sh 5          # run steps 5-10
#   ./run_pipeline.sh 3 7        # run steps 3-7
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

STEPS_DIR="/home/ec2-user/code/steps"
START=${1:-1}
END=${2:-12}

log "============================================================"
log "NeoAntigen Pipeline"
log "Samples: ${TUMOR_SAMPLES[*]}"
log "Normal:  ${NORMAL_SAMPLE}"
log "S3:      ${S3_RESULTS}"
log "Steps:   ${START} → ${END}"
log "============================================================"

run_step() {
    local num="$1"
    local script="$2"
    local desc="$3"

    [[ $num -lt $START || $num -gt $END ]] && return 0

    log ""
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Step ${num}: ${desc}"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    local ext="${script##*.}"
    case "$ext" in
        sh)  bash   "${STEPS_DIR}/${script}" ;;
        py)  "${PYTHON}" "${STEPS_DIR}/${script}" ;;
        R)   Rscript "${STEPS_DIR}/${script}" ;;
        *)   log "[ERROR] Unknown extension: ${ext}"; exit 1 ;;
    esac

    log "Step ${num} done."
}

# ── Steps ─────────────────────────────────────────────────────────────────────
run_step  1  "01_bam_preprocess.sh"   "BAM sort, MarkDuplicates, BQSR, Validate"
run_step  2  "02_mutect2.sh"          "Somatic variant calling (Mutect2 + PoN)"
run_step  3  "03_filter_vcf.sh"       "FilterMutectCalls + SelectVariants (exons)"
run_step  4  "04_annotate.sh"         "SnpEff functional annotation (GRCm38.86)"
run_step  5  "05_pyclone_prep.py"     "Convert annotated VCFs → PyClone YAML input"
run_step  6  "06_pyclone_run.sh"      "PyClone MCMC (beta-binomial, multi-sample)"
run_step  7  "07_pyclone_tables.R"    "PyClone cluster/loci tables + figures"
run_step  8  "08_peptide_extract.py"   "Extract 8/9/10-mer missense/frameshift peptides"
run_step  9  "08b_fusion_peptides.py" "Extract fusion junction peptides (Arriba)"
run_step 10  "09_netmhcpan.sh"        "netMHCpan H-2Kb/H-2Db binding prediction"
run_step 11  "11_rnaseq_expression.R" "RNASeq batch correction + expression matrix"
run_step 12  "10_filter_results.R"    "Rank candidates (binding+clonality+expression)"

log ""
log "============================================================"
log "Pipeline complete. Results: ${RESULTS_DIR}"
log "S3: ${S3_RESULTS}"
log "============================================================"
