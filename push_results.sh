#!/usr/bin/env bash
# =============================================================================
# push_results.sh
#
# Syncs result files from a pipeline run into the git repo and pushes.
# Commits: .html, .pdf, .png, .svg, .tsv, .csv, .txt, .xlsx, .vcf, .md
# Skips:   .bam, .bai, .gz, .vcf.gz, large intermediates
#
# Usage:
#   bash push_results.sh                    # auto-detect latest run
#   bash push_results.sh res_20260311_184352
#   bash push_results.sh res_20260311_184352 --from-s3
# =============================================================================

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_BASE="/home/ec2-user/results"
S3_RESULTS="s3://neoantigen2026-rerun/results"

# --- Determine RUN_ID --------------------------------------------------------
RUN_ID="${1:-}"
FROM_S3=false
[[ "${2:-}" == "--from-s3" ]] && FROM_S3=true

if [[ -z "${RUN_ID}" ]]; then
    # Auto-detect latest local run
    RUN_ID=$(ls -1t "${RESULTS_BASE}" 2>/dev/null | grep '^res_' | head -1 || true)
    if [[ -z "${RUN_ID}" ]]; then
        # Fall back to latest on S3
        RUN_ID=$(aws s3 ls "${S3_RESULTS}/" \
            | awk '{print $NF}' | tr -d '/' | grep '^res_' \
            | sort | tail -1 || true)
        FROM_S3=true
    fi
fi

[[ -z "${RUN_ID}" ]] && { echo "[ERROR] No run ID found. Pass one explicitly."; exit 1; }

echo "[push_results] Run ID : ${RUN_ID}"
echo "[push_results] Source : $( [[ "${FROM_S3}" == true ]] && echo "S3 (${S3_RESULTS}/${RUN_ID}/)" || echo "local (${RESULTS_BASE}/${RUN_ID}/)" )"

REPO_RUN_DIR="${REPO_DIR}/results/${RUN_ID}"
mkdir -p "${REPO_RUN_DIR}"

# --- File types to include ---------------------------------------------------
INCLUDE=(
    --include="*.html"
    --include="*.pdf"
    --include="*.png"
    --include="*.svg"
    --include="*.jpg"
    --include="*.tsv"
    --include="*.csv"
    --include="*.txt"
    --include="*.xlsx"
    --include="*.xls"
    --include="*.md"
    --include="*.vcf"
    --include="*.annotated.vcf"
    --include="*.R"
    --include="*.py"
    --include="*.sh"
)

# --- Sync files into repo ----------------------------------------------------
if [[ "${FROM_S3}" == true ]]; then
    echo "[push_results] Syncing from S3..."
    aws s3 sync "${S3_RESULTS}/${RUN_ID}/" "${REPO_RUN_DIR}/" \
        --exclude="*" \
        "${INCLUDE[@]}" \
        --no-progress
else
    echo "[push_results] Syncing from local..."
    rsync -av --delete \
        --exclude="*.bam" \
        --exclude="*.bai" \
        --exclude="*.gz" \
        --exclude="*.zip" \
        --exclude="*.recal.table" \
        --exclude="*.dup_metrics.txt" \
        --exclude="*.tmp" \
        --include="*.html" \
        --include="*.pdf" \
        --include="*.png" \
        --include="*.svg" \
        --include="*.jpg" \
        --include="*.tsv" \
        --include="*.csv" \
        --include="*.txt" \
        --include="*.xlsx" \
        --include="*.xls" \
        --include="*.md" \
        --include="*.vcf" \
        --include="*.R" \
        --include="*.py" \
        --include="*.sh" \
        --include="*/" \
        --exclude="*" \
        "${RESULTS_BASE}/${RUN_ID}/" \
        "${REPO_RUN_DIR}/"
fi

# --- Check for anything to commit --------------------------------------------
cd "${REPO_DIR}"
git add "results/${RUN_ID}/" .gitignore

if git diff --cached --quiet; then
    echo "[push_results] Nothing new to commit for ${RUN_ID}."
    exit 0
fi

# --- Summarise what's being committed ----------------------------------------
FILE_COUNT=$(git diff --cached --name-only | wc -l)
echo "[push_results] Committing ${FILE_COUNT} file(s)..."
git diff --cached --name-only | sed 's/^/  + /'

# --- Commit and push ---------------------------------------------------------
git commit -m "$(cat <<EOF
Add results for ${RUN_ID}

Includes: reports, tables, figures, HTML outputs, annotated VCFs.
Excludes: BAMs, BAIs, compressed intermediates.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"

git push origin main

echo ""
echo "[push_results] Done. Results for ${RUN_ID} pushed to GitHub."
