#!/usr/bin/env bash
# =============================================================================
# NeoAntigen Pipeline — Main Orchestrator
#
# Sample manifest : code/samples.tsv  (edit there to add/remove samples)
#
# Execution model — optimised for minimal wall-clock time:
#
#   Phase 1 │ Preprocess PoN normals (role=normal,pon) ← needed by PoN + Mutect2
#   ─────────┼──────────────────────────────────────────────────────────────
#   Phase 2a │ Build PoN (step 2, background)       ← starts immediately
#   Phase 2b │ Preprocess tumour samples (step 1)   ← runs in foreground
#            │   As each tumour BAM hits S3, step 2 picks it up and starts
#            │   Mutect2 for that sample without waiting for the others.
#   ─────────┼──────────────────────────────────────────────────────────────
#   Phase 3  │ Wait for all Mutect2 jobs to finish
#   Phase 4  │ Steps 3–12 sequentially
#
# Usage:
#   bash run_pipeline.sh             # full pipeline
#   bash run_pipeline.sh 2           # start from step 2 (normals already done)
#   bash run_pipeline.sh 3 12        # steps 3-12 only
#   bash run_pipeline.sh 1 1         # step 1 (preprocessing) only
# =============================================================================

set -euo pipefail

export RUN_ID="${RUN_ID:-res_$(date '+%Y%m%d_%H%M%S')}"
source /home/ec2-user/code/config.sh

STEPS_DIR="/home/ec2-user/code/steps"
LOG_DIR="${RESULTS_DIR}"
mkdir -p "${LOG_DIR}"

START=${1:-1}
END=${2:-12}

log "============================================================"
log "NeoAntigen Pipeline  RUN_ID=${RUN_ID}"
log "Steps:   ${START} → ${END}"
log "Manifest: ${SAMPLES_META}"
log "Normal:  ${NORMAL_SAMPLE} (SM: ${NORMAL_SM})"
log "PoN:     ${PON_NORMALS[*]}"
log "Tumours: ${TUMOR_SAMPLES[*]}"
log "S3:      ${S3_RESULTS}"
log "============================================================"

# =============================================================================
# Helper: run a step script (bash/python/R) with timing and S3 log upload
# =============================================================================
run_step() {
    local NUM="$1" SCRIPT="$2" DESC="$3"
    [[ ${NUM} -lt ${START} || ${NUM} -gt ${END} ]] && return 0

    log ""
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Step ${NUM}: ${DESC}"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    local T0=$SECONDS

    local EXT="${SCRIPT##*.}"
    case "${EXT}" in
        sh)  bash      "${STEPS_DIR}/${SCRIPT}" ;;
        py)  "${PYTHON}" "${STEPS_DIR}/${SCRIPT}" ;;
        R)   Rscript   "${STEPS_DIR}/${SCRIPT}" ;;
        *)   log "[ERROR] Unknown extension: ${EXT}"; exit 1 ;;
    esac

    log "Step ${NUM} done in $(( SECONDS - T0 ))s."
}

# =============================================================================
# Step 1 — BAM Preprocessing (phased)
# =============================================================================
if [[ ${START} -le 1 && ${END} -ge 1 ]]; then
    log ""
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Step 1 — BAM Preprocessing (phased: normals first, then tumours)"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # ── Phase 1a: preprocess PoN normals ─────────────────────────────────────
    # All samples with role=pon (from samples.tsv) are needed before PoN can
    # be built. Process them first, sequentially (disk constraint), before
    # launching any tumour work.
    log "[Phase 1a] Preprocessing PoN normals: ${PON_NORMALS[*]}"
    SAMPLE_FILTER="${PON_NORMALS[*]}" \
    NON_INTERACTIVE=1 \
        bash "${STEPS_DIR}/01_bam_preprocess.sh"

    log "[Phase 1a] Normals ready."
fi

# =============================================================================
# Step 2 — Somatic Variant Calling (launched in background after normals)
# =============================================================================
MUTECT2_PID=""
if [[ ${START} -le 2 && ${END} -ge 2 ]]; then
    MUTECT2_LOG="${LOG_DIR}/mutect2.log"
    log ""
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    log "Step 2 — Mutect2 launched in background (polls for each tumour BAM)"
    log "  Log: ${MUTECT2_LOG}"
    log "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # Run step 2 in background; it will build PoN immediately then call each
    # tumour as soon as its BAM appears in the permanent S3 store.
    bash "${STEPS_DIR}/02_mutect2.sh" >"${MUTECT2_LOG}" 2>&1 &
    MUTECT2_PID=$!
    log "  Step 2 running (PID ${MUTECT2_PID}) — tail -f ${MUTECT2_LOG}"
fi

# =============================================================================
# Step 1 (continued) — Preprocess tumour samples
# Step 2 picks up each tumour as it finishes and lands on S3.
# =============================================================================
if [[ ${START} -le 1 && ${END} -ge 1 ]]; then
    log ""
    log "[Phase 1b] Preprocessing tumour samples..."
    log "  Step 2 is running in background and will call each tumour as it appears."
    log ""

    SAMPLE_FILTER="${TUMOR_SAMPLES[*]}" \
    NON_INTERACTIVE=1 \
        bash "${STEPS_DIR}/01_bam_preprocess.sh"

    log "[Phase 1b] All tumour BAMs preprocessed."
fi

# =============================================================================
# Wait for Step 2 (Mutect2) to finish before proceeding to downstream steps
# =============================================================================
if [[ -n "${MUTECT2_PID}" ]]; then
    log ""
    log "Waiting for Step 2 (Mutect2, PID ${MUTECT2_PID}) to complete..."
    if wait "${MUTECT2_PID}"; then
        log "Step 2 complete."
    else
        log "[ERROR] Step 2 failed — check ${LOG_DIR}/mutect2.log"
        exit 1
    fi
fi

# =============================================================================
# Steps 3–12 — Sequential downstream analysis
# =============================================================================
run_step  3  "03_filter_vcf.sh"       "FilterMutectCalls + SelectVariants (exons)"
run_step  4  "04_annotate.sh"         "SnpEff functional annotation (GRCm38.86)"
run_step  5  "05_pyclone_prep.py"     "Convert annotated VCFs → PyClone YAML input"
run_step  6  "06_pyclone_run.sh"      "PyClone MCMC (beta-binomial, multi-sample)"
run_step  7  "07_pyclone_tables.R"    "PyClone cluster/loci tables + figures"
run_step  8  "08_peptide_extract.py"  "Extract 8/9/10-mer missense/frameshift peptides"
run_step  9  "08b_fusion_peptides.py" "Extract fusion junction peptides (Arriba)"
run_step 10  "09_netmhcpan.sh"        "netMHCpan H-2Kb/H-2Db binding prediction"
# Step 11 (11_rnaseq_expression.R) has been replaced by the parallel RNASeq
# pipeline (code/rnaseq/run_rnaseq_pipeline.sh), which runs STAR alignment,
# featureCounts quantification, DESeq2 normalization, and limma batch correction
# concurrently with the WES steps above.  run_rnaseq_pipeline.sh also polls for
# the netMHCpan output and triggers step 12 itself once both data sources are
# ready.  If you are running run_pipeline.sh standalone (without the RNASeq
# pipeline), step 12 below will still execute using whatever expression data is
# already present in ${RNASEQ_RESULTS}/expression_matrix.tsv; it handles missing
# expression data gracefully.
run_step 12  "10_filter_results.R"    "Rank candidates (binding + clonality + expression)"

log ""
log "============================================================"
log "Pipeline complete."
log "  Local  : ${RESULTS_DIR}"
log "  S3     : ${S3_RESULTS}"
log "============================================================"
