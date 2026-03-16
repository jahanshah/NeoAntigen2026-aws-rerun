#!/usr/bin/env python3
"""
Step 13 — Rerun PyClone-VI with cnvkit copy number, rebuild ranked table.

Strategy:
  These GemOVCA mouse tumors are highly aneuploid (widespread LOH, near-haploid
  segments), making log2-ratio purity inference unreliable without allele frequency
  data from SNP arrays. The VAF-based purity (median_VAF * 2) is more robust for
  this cohort. This step uses VAF purity, but improves the copy number model for
  each mutation from cnvkit segments (major/minor CN per locus), which meaningfully
  improves PyClone-VI CCF accuracy vs the diploid (2,0) assumption.

Purity used per sample:
  443_D21_new: 0.18  (VAF-based)
  428_D20_new: 0.22  (VAF-based)
  34_D52_old:  0.38  (VAF-based, no WES BAM on S3)
  36_D99_new:  0.22  (VAF-based)
  38_D99_new:  0.20  (VAF-based)
  42_D122_old: 1.00  (VAF-based; cnvkit confirms high purity)

Outputs (same as step 12):
  pyclone/cnvkit/clonality_final.tsv   — CCF per mutation
  pyclone/cnvkit/loci_all.tsv          — per-sample loci
  ranked/neoantigens_final.tsv         — final ranked table (overwritten)
"""

import os, sys, subprocess, logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

BASE_DIR    = "/home/ec2-user"
RESULTS_DIR = os.environ.get("RESULTS_DIR", f"{BASE_DIR}/results/res_20260311_225555")
CNVKIT_DIR  = f"{RESULTS_DIR}/sequenza/cnvkit"
OUT_DIR     = f"{RESULTS_DIR}/pyclone/cnvkit"
RANKED_DIR  = f"{RESULTS_DIR}/ranked"
S3_FINAL    = "s3://neoantigen2026-rerun/final-results"

os.makedirs(OUT_DIR, exist_ok=True)

# VAF-based purity fallback (from median_VAF * 2, clipped to [0.1, 1.0])
VAF_PURITY_FALLBACK = {
    "443_D21_new": 0.18,
    "428_D20_new": 0.22,
    "34_D52_old":  0.38,
    "36_D99_new":  0.22,
    "38_D99_new":  0.20,
    "42_D122_old": 1.00,
}

SAMPLE_BAM_STEMS = {
    "428_D20_new":  "428_D20_new",
    "443_D21_new":  "443_D21_new",
    "36_D99_new":   "36_D99_new",
    "38_D99_new":   "38_D99_new",
    "42_D122_old":  "42_D122_old",
    "34_D52_old":   "34_D52_old",
}


def load_cnvkit_purity(cnvkit_dir: str) -> dict:
    """
    Return VAF-based purity per sample.
    These GemOVCA tumors are too aneuploid for reliable log2-ratio purity
    inference without allele frequency data. VAF-based estimates (median_VAF * 2)
    are used as the primary purity source. cnvkit contributes CN segments.
    """
    purity = {}
    for sample in SAMPLE_BAM_STEMS:
        p = VAF_PURITY_FALLBACK.get(sample)
        if p is not None:
            log.info(f"  {sample}: purity = {p:.3f} (VAF-based)")
            purity[sample] = p
        else:
            log.warning(f"  {sample}: no purity estimate available")
    return purity


# ── Load copy number from cnvkit segments ─────────────────────────────────
def load_cnvkit_cn_segments(cnvkit_dir: str, sample: str) -> pd.DataFrame:
    """
    Load integer copy number calls from cnvkit .call.cns (after cnvkit call).
    Falls back to .cns if .call.cns not present.
    """
    stem = SAMPLE_BAM_STEMS.get(sample, sample)
    for suffix in ["preproc.call.cns", "call.cns", "preproc.cns", ".cns"]:
        path = os.path.join(cnvkit_dir, f"{stem}.{suffix}")
        if os.path.exists(path):
            segs = pd.read_csv(path, sep="\t", comment="@")
            segs = segs.rename(columns={"chromosome": "chrom", "start": "start", "end": "end"})
            segs["chrom"] = segs["chrom"].astype(str)
            return segs
    return pd.DataFrame()


def get_copy_number_cnvkit(mut_id: str, segs: pd.DataFrame) -> tuple:
    """
    Look up copy number for a mutation from cnvkit segments.
    mut_id format: chr1:12345:A:T
    Returns (major_cn, minor_cn).
    Uses cn column if available (from cnvkit call), else infers from log2.
    """
    if segs.empty:
        return 2, 0
    parts = mut_id.split(":")
    if len(parts) < 2:
        return 2, 0
    chrom = parts[0]
    try:
        pos = int(parts[1])
    except ValueError:
        return 2, 0

    hit = segs[
        (segs["chrom"] == chrom) &
        (segs["start"] <= pos) &
        (segs["end"]   >= pos)
    ]
    if hit.empty:
        # Try without chr prefix mismatch
        chrom_bare = chrom.replace("chr", "")
        hit = segs[
            (segs["chrom"].str.replace("chr", "", regex=False) == chrom_bare) &
            (segs["start"] <= pos) &
            (segs["end"]   >= pos)
        ]
    if hit.empty:
        return 2, 0

    row = hit.iloc[0]

    # If cnvkit call was run, there's a 'cn' column
    if "cn" in row.index and pd.notna(row["cn"]):
        total_cn = max(1, int(row["cn"]))
        major = max(1, total_cn - max(0, total_cn // 2 - 1))
        minor = total_cn - major
        return major, max(0, minor)

    # Otherwise infer from log2 ratio
    log2 = row.get("log2", 0.0)
    total_cn = max(1, round(2.0 * 2.0 ** log2))
    major = max(1, (total_cn + 1) // 2)
    minor = max(0, total_cn - major)
    return major, minor


# ── Run PyClone-VI per sample ──────────────────────────────────────────────
def run_pyclone_sample(sample: str, s_df: pd.DataFrame,
                       purity: float, out_dir: str) -> pd.DataFrame:
    in_tsv   = f"{out_dir}/{sample}_input.tsv"
    out_h5   = f"{out_dir}/{sample}.h5"
    loci_tsv = f"{out_dir}/{sample}_loci.tsv"

    s_df["tumour_content"] = purity
    s_df.to_csv(in_tsv, sep="\t", index=False)

    pyclone = "/home/ec2-user/miniforge3/bin/pyclone-vi"

    ret = subprocess.run([
        pyclone, "fit",
        "--in-file",      in_tsv,
        "--out-file",     out_h5,
        "--num-clusters", "10",
        "--density",      "beta-binomial",
        "--num-restarts", "10",
    ], capture_output=True, text=True)

    if ret.returncode != 0:
        log.error(f"  PyClone-VI failed for {sample}: {ret.stderr[-300:]}")
        return pd.DataFrame()

    subprocess.run([
        pyclone, "write-results-file",
        "--in-file",  out_h5,
        "--out-file", loci_tsv,
    ], capture_output=True)

    result = pd.read_csv(loci_tsv, sep="\t")
    n_clonal = (result["cellular_prevalence"] > 0.8).sum()
    log.info(f"  {sample}: {len(result)} loci, purity={purity:.3f}, "
             f"clonal(>0.8)={n_clonal}")
    return result


# ── Rebuild ranked table ───────────────────────────────────────────────────
def rebuild_ranked(clonality: pd.DataFrame, results_dir: str):
    ranked = pd.read_csv(f"{results_dir}/ranked/neoantigens_full_annotated.tsv",
                         sep="\t")
    ccf_map = dict(zip(clonality["mutation_id"], clonality["mean_ccf"]))

    def get_ccf(row):
        if row.get("is_fusion") or "::" in str(row.get("mut_id", "")):
            return row["clone_pct"]
        return ccf_map.get(str(row["mut_id"]), row["clone_pct"])

    ranked["clone_pct"] = ranked.apply(get_ccf, axis=1)

    max_ccf = ranked["clone_pct"].max()
    ranked["clone_pct_norm"] = (ranked["clone_pct"] / max_ccf).clip(0, 1)

    def mhc_score(rank_el, bclass):
        if pd.isna(rank_el) or bclass == "NB": return 0.0
        if bclass == "SB": return max(0.5, 1.0 - rank_el / 2.0)
        if bclass == "WB": return 1.0 - (rank_el - 0.5) / 1.5 * 0.5
        if bclass == "PB": return 0.5 - (rank_el - 2.0) / 8.0 * 0.4
        return 0.0

    ranked["mhc_score"] = ranked.apply(
        lambda r: mhc_score(r["rank_el_H2Db"], r["binding_class"]), axis=1)
    max_expr = ranked["log2cpm"].max()
    ranked["expr_pct"] = ranked["log2cpm"].fillna(0) / max_expr

    ranked["composite_score"] = (
        0.40 * ranked["clone_pct_norm"] +
        0.40 * ranked["mhc_score"] +
        0.20 * ranked["expr_pct"]
    ) * ranked["rna_validated"].map({True: 1.0, False: 0.3})

    ranked = ranked.sort_values("composite_score", ascending=False).reset_index(drop=True)
    ranked["rank"] = ranked.index + 1
    ranked.drop(columns=["clone_pct_norm"], inplace=True, errors="ignore")

    out = f"{results_dir}/ranked/neoantigens_final.tsv"
    ranked.to_csv(out, sep="\t", index=False)
    binders = ranked[ranked["binding_class"].isin(["SB", "WB", "PB"])]
    binders.to_csv(f"{results_dir}/ranked/strong_binders_final.tsv", sep="\t", index=False)

    log.info(f"Saved final ranked table: {out} ({len(ranked)} rows, "
             f"{len(binders)} binders)")
    log.info("\nTop 15 candidates:")
    top = binders.head(15)[["rank", "gene", "sample", "peptide", "binding_class",
                             "rank_el_H2Db", "log2cpm", "clone_pct", "composite_score"]]
    print(top.to_string(index=False))
    return ranked


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    log.info("Step 13 — cnvkit purity + PyClone-VI rerun")

    # 1. Extract purity from cnvkit output
    purity_map = load_cnvkit_purity(CNVKIT_DIR)
    if not purity_map:
        log.error(f"No cnvkit .cns files found in {CNVKIT_DIR}. "
                  "Run cnvkit batch first.")
        sys.exit(1)

    log.info(f"Purity estimates ({len(purity_map)} samples):")
    for s, p in sorted(purity_map.items()):
        log.info(f"  {s}: {p:.3f}")

    # 2. Load PyClone-VI input
    vi_input = pd.read_csv(
        f"{RESULTS_DIR}/pyclone/rerun/pyclone_vi_input.tsv", sep="\t")

    loci_all = []

    for sample, purity in purity_map.items():
        s_df = vi_input[vi_input["sample_id"] == sample].copy()
        if s_df.empty:
            log.warning(f"  No mutations for {sample} in input TSV")
            continue

        # 3. Update copy number from cnvkit segments
        segs = load_cnvkit_cn_segments(CNVKIT_DIR, sample)
        if not segs.empty:
            log.info(f"  {sample}: updating CN from {len(segs)} cnvkit segments")
            cn_vals = s_df["mutation_id"].apply(
                lambda m: get_copy_number_cnvkit(m, segs))
            s_df["major_cn"] = cn_vals.apply(lambda x: x[0])
            s_df["minor_cn"] = cn_vals.apply(lambda x: x[1])

        # 4. Run PyClone-VI
        result = run_pyclone_sample(sample, s_df, purity, OUT_DIR)
        if not result.empty:
            loci_all.append(result)

    if not loci_all:
        log.error("No PyClone-VI results produced.")
        sys.exit(1)

    merged = pd.concat(loci_all, ignore_index=True)
    merged.to_csv(f"{OUT_DIR}/loci_all.tsv", sep="\t", index=False)

    # 5. Clonality summary
    clonality = merged.groupby("mutation_id").agg(
        mean_ccf=("cellular_prevalence", "mean"),
        max_ccf =("cellular_prevalence", "max"),
        cluster_id=("cluster_id", "first"),
    ).reset_index()
    clonality.to_csv(f"{OUT_DIR}/clonality_final.tsv", sep="\t", index=False)

    log.info(f"\nClonality summary ({len(clonality)} mutations):")
    log.info(f"  Clonal (CCF>0.8):     {(clonality['max_ccf']>0.8).sum()}")
    log.info(f"  Subclonal (0.2-0.8):  {((clonality['max_ccf']>=0.2)&(clonality['max_ccf']<=0.8)).sum()}")
    log.info(f"  Minor (<0.2):         {(clonality['max_ccf']<0.2).sum()}")

    # 6. Rebuild ranked table
    rebuild_ranked(clonality, RESULTS_DIR)

    # 7. Upload to S3
    os.system(f"aws s3 cp {OUT_DIR}/loci_all.tsv "
              f"{S3_FINAL}/pyclone/loci_cnvkit.tsv --quiet")
    os.system(f"aws s3 cp {OUT_DIR}/clonality_final.tsv "
              f"{S3_FINAL}/pyclone/clonality_cnvkit.tsv --quiet")
    os.system(f"aws s3 cp {RANKED_DIR}/neoantigens_final.tsv "
              f"{S3_FINAL}/ranked/neoantigens_final.tsv --quiet")
    os.system(f"aws s3 cp {RANKED_DIR}/strong_binders_final.tsv "
              f"{S3_FINAL}/ranked/strong_binders_final.tsv --quiet")
    log.info("Step 13 complete. Results uploaded to S3.")


if __name__ == "__main__":
    main()
