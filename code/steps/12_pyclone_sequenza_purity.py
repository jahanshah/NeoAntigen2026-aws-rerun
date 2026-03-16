#!/usr/bin/env python3
"""
Step 12 — PyClone-VI with Sequenza-derived Tumor Purity

Reads purity estimates from Sequenza output per sample, then reruns
PyClone-VI per sample with the correct tumour_content value.

Also integrates CNV (major/minor copy number) from Sequenza segments
for per-mutation copy number, removing the diploid assumption.

Outputs:
  pyclone/sequenza/clonality_final.tsv     — CCF per mutation
  pyclone/sequenza/loci_all.tsv            — per-sample loci
  ranked/neoantigens_final.tsv             — final ranked table
"""

import os, sys, subprocess, logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

BASE_DIR    = "/home/ec2-user"
RESULTS_DIR = os.environ.get("RESULTS_DIR", f"{BASE_DIR}/results/res_20260311_225555")
SEQ_DIR     = f"{RESULTS_DIR}/sequenza"
OUT_DIR     = f"{RESULTS_DIR}/pyclone/sequenza"
RANKED_DIR  = f"{RESULTS_DIR}/ranked"
S3_FINAL    = "s3://neoantigen2026-rerun/final-results"

os.makedirs(OUT_DIR, exist_ok=True)


# ── Load Sequenza purity results ───────────────────────────────────────────
def load_sequenza_purity(seq_dir: str) -> dict:
    """
    Read *_confints_CP.txt from each sample's Sequenza output directory.
    Returns dict: sample_id → purity (float)
    """
    purity = {}
    for sample in os.listdir(seq_dir):
        sample_dir = os.path.join(seq_dir, sample)
        if not os.path.isdir(sample_dir):
            continue
        # Sequenza writes <sample>_confints_CP.txt with columns:
        #   name, cellularity.mean, cellularity.025, cellularity.975, ploidy.mean, ...
        cp_file = os.path.join(sample_dir, f"{sample}_confints_CP.txt")
        if not os.path.exists(cp_file):
            log.warning(f"  No confints_CP file for {sample}: {cp_file}")
            continue
        cp = pd.read_csv(cp_file, sep="\t")
        # Row with name == "confint" gives the maximum likelihood estimate
        best_row = cp[cp["name"] == "confint"]
        if best_row.empty:
            best_row = cp.head(1)
        p = float(best_row["cellularity.mean"].iloc[0])
        pl = float(best_row["ploidy.mean"].iloc[0])
        purity[sample] = p
        log.info(f"  {sample}: purity={p:.3f}  ploidy={pl:.2f}")
    return purity


# ── Load Sequenza CNV segments ─────────────────────────────────────────────
def load_sequenza_segments(seq_dir: str, sample: str) -> pd.DataFrame:
    """
    Load copy number segments from Sequenza *_segments.txt.
    Columns used: chromosome, start.pos, end.pos, CNt (total CN), A (major), B (minor)
    """
    seg_file = os.path.join(seq_dir, sample, f"{sample}_segments.txt")
    if not os.path.exists(seg_file):
        return pd.DataFrame()
    segs = pd.read_csv(seg_file, sep="\t")
    # Rename columns to standard names
    segs = segs.rename(columns={
        "chromosome": "chrom",
        "start.pos":  "start",
        "end.pos":    "end",
        "A":          "major_cn",
        "B":          "minor_cn",
        "CNt":        "total_cn",
    })
    segs["chrom"] = segs["chrom"].astype(str)
    return segs


def get_copy_number(mut_id: str, segs: pd.DataFrame) -> tuple:
    """
    Look up copy number for a mutation from Sequenza segments.
    mut_id format: chr1:12345:A:T
    Returns (major_cn, minor_cn) or (2, 0) if not found.
    """
    if segs.empty:
        return 2, 0
    parts = mut_id.split(":")
    if len(parts) < 2:
        return 2, 0
    chrom = parts[0].replace("chr", "")
    try:
        pos = int(parts[1])
    except ValueError:
        return 2, 0
    hit = segs[
        (segs["chrom"].str.replace("chr", "") == chrom) &
        (segs["start"] <= pos) &
        (segs["end"]   >= pos)
    ]
    if hit.empty:
        return 2, 0
    row = hit.iloc[0]
    major = max(1, int(row.get("major_cn", 2)))
    minor = max(0, int(row.get("minor_cn", 0)))
    return major, minor


# ── Run PyClone-VI per sample ──────────────────────────────────────────────
def run_pyclone_sample(sample: str, s_df: pd.DataFrame,
                       purity: float, out_dir: str) -> pd.DataFrame:
    in_tsv   = f"{out_dir}/{sample}_input.tsv"
    out_h5   = f"{out_dir}/{sample}.h5"
    loci_tsv = f"{out_dir}/{sample}_loci.tsv"

    s_df["tumour_content"] = purity
    s_df.to_csv(in_tsv, sep="\t", index=False)

    ret = subprocess.run([
        "pyclone-vi", "fit",
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
        "pyclone-vi", "write-results-file",
        "--in-file",  out_h5,
        "--out-file", loci_tsv,
    ], capture_output=True)

    result = pd.read_csv(loci_tsv, sep="\t")
    n_clonal = (result["cellular_prevalence"] > 0.8).sum()
    log.info(f"  {sample}: {len(result)} loci, purity={purity:.3f}, "
             f"clonal(>0.8)={n_clonal}")
    return result


# ── Update ranked table ────────────────────────────────────────────────────
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
    log.info("Step 12 — PyClone-VI with Sequenza purity")

    # Load Sequenza purity estimates
    purity_map = load_sequenza_purity(SEQ_DIR)
    if not purity_map:
        log.error(f"No Sequenza results found in {SEQ_DIR}. "
                  "Run 11_sequenza_purity.sh first.")
        sys.exit(1)

    log.info(f"Loaded purity for {len(purity_map)} samples: {purity_map}")

    # Load PyClone-VI input
    vi_input = pd.read_csv(
        f"{RESULTS_DIR}/pyclone/rerun/pyclone_vi_input.tsv", sep="\t")

    loci_all = []

    for sample, purity in purity_map.items():
        s_df = vi_input[vi_input["sample_id"] == sample].copy()
        if s_df.empty:
            log.warning(f"  No mutations for {sample} in input TSV")
            continue

        # Update copy number from Sequenza segments if available
        segs = load_sequenza_segments(SEQ_DIR, sample)
        if not segs.empty:
            log.info(f"  {sample}: updating CN from {len(segs)} Sequenza segments")
            cn_vals = s_df["mutation_id"].apply(
                lambda m: get_copy_number(m, segs))
            s_df["major_cn"] = cn_vals.apply(lambda x: x[0])
            s_df["minor_cn"] = cn_vals.apply(lambda x: x[1])

        result = run_pyclone_sample(sample, s_df, purity, OUT_DIR)
        if not result.empty:
            loci_all.append(result)

    if not loci_all:
        log.error("No PyClone-VI results produced.")
        sys.exit(1)

    merged = pd.concat(loci_all, ignore_index=True)
    merged.to_csv(f"{OUT_DIR}/loci_all.tsv", sep="\t", index=False)

    # Clonality summary
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

    # Rebuild ranked table
    rebuild_ranked(clonality, RESULTS_DIR)

    # Upload to S3
    os.system(f"aws s3 cp {OUT_DIR}/loci_all.tsv "
              f"{S3_FINAL}/pyclone/loci_sequenza.tsv --quiet")
    os.system(f"aws s3 cp {OUT_DIR}/clonality_final.tsv "
              f"{S3_FINAL}/pyclone/clonality_sequenza.tsv --quiet")
    os.system(f"aws s3 cp {RANKED_DIR}/neoantigens_final.tsv "
              f"{S3_FINAL}/ranked/neoantigens_final.tsv --quiet")
    os.system(f"aws s3 cp {RANKED_DIR}/strong_binders_final.tsv "
              f"{S3_FINAL}/ranked/strong_binders_final.tsv --quiet")
    log.info("Step 12 complete. Results uploaded to S3.")


if __name__ == "__main__":
    main()
