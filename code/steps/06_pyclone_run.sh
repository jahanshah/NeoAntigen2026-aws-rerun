#!/usr/bin/env bash
# =============================================================================
# Step 6 — PyClone MCMC (beta-binomial, multi-sample)
# Runs clonal population inference across all tumour timepoints.
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

INDIR="${RESULTS_DIR}/pyclone/input"
OUTDIR="${RESULTS_DIR}/pyclone/output"
PYCLONE_CFG="${OUTDIR}/config.yaml"
S3_OUT="${S3_RESULTS}/pyclone"
PYCLONE="${PYTHON} -m pyclone"

mkdir -p "${OUTDIR}"

# Collect YAML files
YAML_FILES=()
for TUMOR in "${TUMOR_SAMPLES[@]}"; do
    YAML="${INDIR}/${TUMOR}.yaml"
    if [[ ! -f "${YAML}" ]]; then
        aws s3 cp "${S3_OUT}/input/${TUMOR}.yaml" "${YAML}" 2>/dev/null || {
            log "[WARN] Missing YAML for ${TUMOR} — skipping"
            continue
        }
    fi
    YAML_FILES+=("${YAML}")
done

[[ ${#YAML_FILES[@]} -eq 0 ]] && { log "[ERROR] No YAML input files found"; exit 1; }
log "Found ${#YAML_FILES[@]} sample YAMLs"

# Check if already done
if aws s3 ls "${S3_OUT}/output/tables/cluster.tsv" &>/dev/null; then
    skip "PyClone output (S3)"
    exit 0
fi

# Build PyClone config
log "Building PyClone config: ${PYCLONE_CFG}"
cat > "${PYCLONE_CFG}" << YAML
working_dir: ${OUTDIR}
trace_dir:   ${OUTDIR}/trace
num_iters:   10000
burn_in:     1000
thin:        1
base_measure_params:
  alpha: 1.0
  beta:  1.0
beta_binomial_precision_params:
  value: 1000
  prior:
    shape:  1.0
    rate:   0.0001
density: beta_binomial
init_method: disconnected
samples:
$(for YAML in "${YAML_FILES[@]}"; do
    SNAME=$(basename "${YAML}" .yaml)
    printf "  %s:\n    tumour_content:\n      value: 1.0\n    error_rate: 0.001\n    mutations_file: %s\n" \
        "${SNAME}" "${YAML}"
done)
YAML

# Run MCMC
log "Running PyClone MCMC (multi-sample, beta-binomial)..."
log "  Samples: ${YAML_FILES[*]##*/}"
log "  Iterations: 10000  burn-in: 1000"

${PYCLONE} run_analysis_pipeline \
    --config_file "${PYCLONE_CFG}" \
    --num_processors ${THREADS} 2>&1

log "PyClone complete. Output: ${OUTDIR}"

# Upload outputs
aws s3 sync "${OUTDIR}/" "${S3_OUT}/output/" --exclude "*.bam" 2>&1 | tail -3

log "Step 6 complete."
