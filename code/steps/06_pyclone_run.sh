#!/usr/bin/env bash
# =============================================================================
# Step 6 — PyClone-VI clonal inference (multi-sample)
# Uses pyclone-vi (variational inference, beta-binomial) instead of the legacy
# MCMC-based PyClone, which requires Python 2 and is no longer maintained.
#
# Input:  ${RESULTS_DIR}/pyclone/input/all_samples_mutations.tsv  (from step 5)
#         per-sample TSVs: ${RESULTS_DIR}/pyclone/input/<sample>.mutations.tsv
# Output: ${RESULTS_DIR}/pyclone/output/results.h5    — fitted model
#         ${RESULTS_DIR}/pyclone/output/tables/loci.tsv
#         ${RESULTS_DIR}/pyclone/output/tables/cluster.tsv
# =============================================================================

set -euo pipefail
source /home/ec2-user/code/config.sh

INDIR="${RESULTS_DIR}/pyclone/input"
OUTDIR="${RESULTS_DIR}/pyclone/output"
TABLES="${OUTDIR}/tables"
S3_OUT="${S3_RESULTS}/pyclone"
PYCLONE_VI="${CONDA_BIN}/pyclone-vi"

mkdir -p "${OUTDIR}" "${TABLES}"

# Check if already done
if aws s3 ls "${S3_OUT}/output/tables/loci.tsv" &>/dev/null; then
    skip "PyClone-VI output (S3)"
    exit 0
fi

# ── 1. Build pyclone-vi input TSV from step 5 mutation files ─────────────────
# pyclone-vi multi-sample requires every mutation to have an entry for every
# sample. Mutations not detected in a sample get alt_counts=0 and
# ref_counts=depth (using median depth from samples where it was detected).
PCVI_INPUT="${OUTDIR}/pyclone_vi_input.tsv"

log "Building pyclone-vi input TSV from per-sample mutation files..."

# Collect missing TSVs from S3
for TUMOR in "${TUMOR_SAMPLES[@]}"; do
    MUT_TSV="${INDIR}/${TUMOR}.mutations.tsv"
    if [[ ! -f "${MUT_TSV}" ]]; then
        aws s3 cp "${S3_OUT}/input/${TUMOR}.mutations.tsv" "${MUT_TSV}" 2>/dev/null || {
            log "[WARN] Missing mutations TSV for ${TUMOR} — skipping"
        }
    fi
done

# Build complete input via Python (handles cross-sample filling)
"${PYTHON}" - << INNEREOF
import os, sys
import pandas as pd

indir  = "${INDIR}"
outdir = "${OUTDIR}"
tumor_samples = "${TUMOR_SAMPLES[@]}".split()

all_per_sample = {}
for samp in tumor_samples:
    tsv = f"{indir}/{samp}.mutations.tsv"
    if not os.path.exists(tsv):
        print(f"  [SKIP] {samp}: TSV not found", flush=True)
        continue
    df = pd.read_csv(tsv, sep="\t")
    df = df[df["alt_count"] >= 3].copy()
    df["mutation_id"] = df["chrom"].astype(str) + ":" + df["pos"].astype(str) + ":" + df["ref"] + ":" + df["alt"]
    all_per_sample[samp] = df
    print(f"  {samp}: {len(df)} mutations", flush=True)

if not all_per_sample:
    print("ERROR: No sample data found", flush=True)
    sys.exit(1)

# All unique mutation IDs across all samples
all_mutations = pd.concat([
    df[["mutation_id","ref_count","alt_count","depth"]].assign(sample_id=s)
    for s, df in all_per_sample.items()
]).reset_index(drop=True)

# Median depth per mutation (for filling 0-count rows)
med_depth = all_mutations.groupby("mutation_id")["depth"].median().reset_index()
med_depth.columns = ["mutation_id","med_depth"]

all_mut_ids = all_mutations["mutation_id"].unique()
all_samps   = list(all_per_sample.keys())

rows = []
for mut_id in all_mut_ids:
    md = int(med_depth.loc[med_depth["mutation_id"]==mut_id, "med_depth"].values[0])
    for samp in all_samps:
        sub = all_mutations[(all_mutations["mutation_id"]==mut_id) & (all_mutations["sample_id"]==samp)]
        if len(sub) > 0:
            rc = int(sub.iloc[0]["ref_count"])
            ac = int(sub.iloc[0]["alt_count"])
        else:
            rc = md  # not detected: use median depth as ref_count, 0 alt
            ac = 0
        rows.append({"mutation_id": mut_id, "sample_id": samp,
                     "ref_counts": rc, "alt_counts": ac,
                     "normal_cn": 2, "major_cn": 2, "minor_cn": 0, "tumour_content": 1.0})

out = pd.DataFrame(rows)
out.to_csv(f"{outdir}/pyclone_vi_input.tsv", sep="\t", index=False)
print(f"pyclone-vi input: {len(out)} rows, {len(all_mut_ids)} mutations, {len(all_samps)} samples", flush=True)
INNEREOF

N_MUTS=$(awk 'NR>1' "${PCVI_INPUT}" | wc -l)
N_UNIQUE=$(awk 'NR>1{print $1}' "${PCVI_INPUT}" | sort -u | wc -l)
log "pyclone-vi input: ${N_MUTS} rows, ${N_UNIQUE} unique mutations"
[[ ${N_MUTS} -eq 0 ]] && { log "[ERROR] No mutations in pyclone-vi input"; exit 1; }

# ── 2. Run pyclone-vi fit ─────────────────────────────────────────────────────
H5_OUT="${OUTDIR}/results.h5"
log "Running pyclone-vi fit (beta-binomial, 10 clusters, 8 restarts)..."
"${PYCLONE_VI}" fit \
    -i "${PCVI_INPUT}" \
    -o "${H5_OUT}" \
    -c 10 \
    -d beta-binomial \
    -g 100 \
    -r 8 \
    -t "${THREADS}" 2>&1

log "pyclone-vi fit complete."

# ── 3. Write results to TSV ───────────────────────────────────────────────────
LOCI_RAW="${OUTDIR}/loci_raw.tsv"
log "Writing results to TSV..."
"${PYCLONE_VI}" write-results-file \
    -i "${H5_OUT}" \
    -o "${LOCI_RAW}" 2>&1

log "Results written: ${LOCI_RAW}"

# ── 4. Reformat to loci.tsv / cluster.tsv compatible with step 7 ─────────────
# loci_raw columns: mutation_id, sample_id, cluster_id, cellular_prevalence,
#                   cellular_prevalence_std, variant_allele_frequency
log "Reformatting to loci.tsv + cluster.tsv..."

"${PYTHON}" - << 'PYEOF'
import os, sys
import pandas as pd

results_dir = os.environ["RESULTS_DIR"]
outdir  = f"{results_dir}/pyclone/output"
tables  = f"{outdir}/tables"
os.makedirs(tables, exist_ok=True)

# pyclone-vi columns: mutation_id, sample_id, cluster_id, cellular_prevalence,
#                     cellular_prevalence_std, cluster_assignment_prob
df = pd.read_csv(f"{outdir}/loci_raw.tsv", sep="\t")
print(f"loci_raw columns: {list(df.columns)}", flush=True)

# Pivot: one row per mutation, cellular_prevalence as wide columns
loci_wide = df.pivot_table(
    index  =["mutation_id","cluster_id"],
    columns="sample_id",
    values ="cellular_prevalence"
).reset_index()
loci_wide.columns = [
    c if c in ("mutation_id","cluster_id") else f"{c}_cellular_prevalence"
    for c in loci_wide.columns
]
loci_wide.to_csv(f"{tables}/loci.tsv", sep="\t", index=False)
print(f"loci.tsv: {len(loci_wide)} rows", flush=True)

# cluster.tsv — one row per mutation with its cluster_id
cluster = df[["mutation_id","cluster_id"]].drop_duplicates()
cluster.to_csv(f"{tables}/cluster.tsv", sep="\t", index=False)
print(f"cluster.tsv: {len(cluster)} rows", flush=True)

# Cluster summary
clust_summary = cluster.groupby("cluster_id").size().reset_index(name="n_mutations")
clust_summary.to_csv(f"{tables}/cluster_summary.tsv", sep="\t", index=False)
print(f"cluster_summary.tsv: {len(clust_summary)} clusters", flush=True)

# Also save the raw long-format table for step 7 (cellular prevalence over time)
df.to_csv(f"{tables}/prevalence_long.tsv", sep="\t", index=False)
print(f"prevalence_long.tsv: {len(df)} rows", flush=True)
PYEOF

log "Reformatted tables written to ${TABLES}/"

# ── 5. Upload outputs ─────────────────────────────────────────────────────────
log "Uploading to S3..."
aws s3 sync "${OUTDIR}/" "${S3_OUT}/output/" 2>&1 | tail -5

log "Step 6 complete."
