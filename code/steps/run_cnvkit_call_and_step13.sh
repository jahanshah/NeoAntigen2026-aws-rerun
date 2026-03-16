#!/bin/bash
# Waits for cnvkit batch to finish, then runs cnvkit call + step 13
set -euo pipefail
export PATH="/home/ec2-user/miniforge3/bin:$PATH"

CNVKIT_DIR="/home/ec2-user/results/res_20260311_225555/sequenza/cnvkit"
RESULTS_DIR="/home/ec2-user/results/res_20260311_225555"
LOG="${RESULTS_DIR}/logs/cnvkit_call_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "${RESULTS_DIR}/logs"

exec > >(tee -a "${LOG}") 2>&1

echo "[$(date)] Waiting for cnvkit batch (PID 296018) to finish..."
while kill -0 296018 2>/dev/null; do
    sleep 30
    echo "[$(date)] cnvkit batch still running..."
done

echo "[$(date)] cnvkit batch done. Listing outputs:"
ls -lh "${CNVKIT_DIR}/"*.cns 2>/dev/null || echo "No .cns files yet"

# Run cnvkit call on each .cns to get integer copy numbers
echo "[$(date)] Running cnvkit call..."
for cns in "${CNVKIT_DIR}"/*.cns; do
    [ -f "${cns}" ] || continue
    # Skip if already a .call.cns
    [[ "${cns}" == *.call.cns ]] && continue
    base="${cns%.cns}"
    echo "[$(date)]   cnvkit call: ${cns}"
    cnvkit.py call \
        "${cns}" \
        --output "${base}.call.cns" \
        --method clonal \
        2>/dev/null || echo "[$(date)]   WARNING: cnvkit call failed for ${cns}"
done

echo "[$(date)] cnvkit call complete. Running step 13..."
cd /home/ec2-user
python3 /home/ec2-user/NeoAntigen2026-aws-rerun/code/steps/13_cnvkit_purity_pyclone.py

echo "[$(date)] Step 13 complete."

# Also add 34_D52_old if its BAM downloaded and batch missed it
BAM_34="/home/ec2-user/tmp/neoantig_pipeline/sequenza_bam/34_D52_old.preproc.bam"
if [ -f "${BAM_34}" ] && [ ! -f "${CNVKIT_DIR}/34_D52_old.preproc.cns" ]; then
    echo "[$(date)] Running cnvkit for 34_D52_old separately..."
    REF="${CNVKIT_DIR}/cnvkit_reference.cnn"
    if [ -f "${REF}" ]; then
        cnvkit.py coverage "${BAM_34}" \
            "${RESULTS_DIR}/sequenza/targets_auto.bed" \
            -o "${CNVKIT_DIR}/34_D52_old.preproc.targetcoverage.cnn" -p 8
        cnvkit.py coverage "${BAM_34}" \
            "${RESULTS_DIR}/sequenza/antitargets_auto.bed" \
            -o "${CNVKIT_DIR}/34_D52_old.preproc.antitargetcoverage.cnn" -p 8
        cnvkit.py fix \
            "${CNVKIT_DIR}/34_D52_old.preproc.targetcoverage.cnn" \
            "${CNVKIT_DIR}/34_D52_old.preproc.antitargetcoverage.cnn" \
            "${REF}" \
            -o "${CNVKIT_DIR}/34_D52_old.preproc.cnr"
        cnvkit.py segment \
            "${CNVKIT_DIR}/34_D52_old.preproc.cnr" \
            -o "${CNVKIT_DIR}/34_D52_old.preproc.cns" -p 8
        cnvkit.py call \
            "${CNVKIT_DIR}/34_D52_old.preproc.cns" \
            -o "${CNVKIT_DIR}/34_D52_old.preproc.call.cns" --method clonal
        echo "[$(date)] 34_D52_old cnvkit complete. Rerunning step 13..."
        python3 /home/ec2-user/NeoAntigen2026-aws-rerun/code/steps/13_cnvkit_purity_pyclone.py
    else
        echo "[$(date)] Reference .cnn not found, skipping 34_D52_old"
    fi
fi

echo "[$(date)] All done."
