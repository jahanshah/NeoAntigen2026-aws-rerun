#!/usr/bin/env python3
"""
Step 10 — Neoantigen Peptide Ranking

Ranks peptides by a composite score combining three biological signals:
  1. Clonality    — PyClone-VI cellular prevalence (higher = more clonal = better target)
  2. Rank_EL%     — netMHCpan eluted-ligand rank   (lower = stronger binder = better)
  3. Expression   — RNA-seq gene CPM               (higher = more expressed = better)

Each signal is converted to a percentile score [0–1], then combined:
  composite = w_clone * clonality_pct + w_expr * expression_pct + w_mhc * mhc_pct

Usage:
  python3 10_rank_neoantigens.py [--mhc MHC_TSV] [--out OUT_TSV]

  --mhc   Path to merged netMHCpan output TSV (optional; ranking works without it)
  --out   Output path (default: RESULTS_DIR/ranked/ranked_neoantigens.tsv)

netMHCpan TSV expected columns (tab-separated, one row per peptide-allele):
  Peptide, MHC, %Rank_EL  (additional columns ignored)
"""

import os, re, gzip, logging, argparse, sys
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR    = "/home/ec2-user"
RESULTS_WES = os.environ.get("RESULTS_DIR",  f"{BASE_DIR}/results/res_20260311_225555")
S3_WES      = os.environ.get("S3_RESULTS",   "s3://neoantigen2026-rerun/results/res_20260311_225555")
S3_RNA      = "s3://neoantigen2026-rerun/results/res_20260313_071626/rnaseq"
S3_FINAL    = "s3://neoantigen2026-rerun/final-results"

PEPTIDES_TSV  = f"{RESULTS_WES}/peptides/all_peptides.tsv"
PYCLONE_LOCI  = f"{RESULTS_WES}/pyclone/output/tables/loci.tsv"
COUNTS_TXT    = f"{RESULTS_WES}/rnaseq/raw_counts.txt"
GTF_GZ        = f"{BASE_DIR}/.cache/pyensembl/GRCm38/ensembl102/Mus_musculus.GRCm38.102.gtf.gz"
OUT_DIR       = f"{RESULTS_WES}/ranked"

# Weights for composite score (must sum to 1 when MHC available; clonality+expression otherwise)
W_CLONE = 0.40
W_MHC   = 0.40
W_EXPR  = 0.20

TUMOR_SAMPLES = ["443_D21_new", "428_D20_new", "34_D52_old",
                 "36_D99_new",  "38_D99_new",  "42_D122_old"]

# ── Helpers ───────────────────────────────────────────────────────────────────

def pct_rank(series: pd.Series) -> pd.Series:
    """Convert a series to percentile ranks in [0, 1]; NaN stays NaN."""
    return series.rank(pct=True, na_option="keep")


def fetch_s3(s3_path: str, local_path: str):
    """Download from S3 if local file is missing."""
    if not os.path.exists(local_path):
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        ret = os.system(f"aws s3 cp {s3_path} {local_path} --quiet")
        if ret != 0:
            raise FileNotFoundError(f"Failed to fetch {s3_path}")


# ── 1. Load peptides ───────────────────────────────────────────────────────────

def load_peptides() -> pd.DataFrame:
    if not os.path.exists(PEPTIDES_TSV):
        fetch_s3(f"{S3_WES}/peptides/all_peptides.tsv", PEPTIDES_TSV)
    df = pd.read_csv(PEPTIDES_TSV, sep="\t")
    log.info(f"Peptides loaded: {len(df)} rows, {df['peptide'].nunique()} unique sequences")
    return df


# ── 2. Load clonality ─────────────────────────────────────────────────────────

def load_clonality() -> pd.DataFrame:
    """
    Returns DataFrame with columns:
      mut_id, cluster_id, mean_prevalence, max_prevalence,
      n_samples_present, <sample>_cellular_prevalence ...
    """
    if not os.path.exists(PYCLONE_LOCI):
        fetch_s3(f"{S3_WES}/pyclone/output/tables/loci.tsv", PYCLONE_LOCI)

    loci = pd.read_csv(PYCLONE_LOCI, sep="\t")
    loci = loci.rename(columns={"mutation_id": "mut_id"})

    prev_cols = [c for c in loci.columns if c.endswith("_cellular_prevalence")]
    loci["mean_prevalence"] = loci[prev_cols].mean(axis=1)
    loci["max_prevalence"]  = loci[prev_cols].max(axis=1)
    loci["n_samples_present"] = (loci[prev_cols] > 0.05).sum(axis=1)

    log.info(f"PyClone loci: {len(loci)} mutations, clusters: {sorted(loci['cluster_id'].unique())}")
    return loci


# ── 3. Load expression ────────────────────────────────────────────────────────

def _build_ensembl_to_symbol(gtf_gz: str) -> dict:
    """Parse GTF gene features to build Ensembl gene ID → gene symbol map."""
    mapping = {}
    log.info("Building Ensembl ID → gene symbol map from GTF...")
    opener = gzip.open if gtf_gz.endswith(".gz") else open
    gid_re  = re.compile(r'gene_id "([^"]+)"')
    gsym_re = re.compile(r'gene_name "([^"]+)"')
    with opener(gtf_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            gid  = gid_re.search(line)
            gsym = gsym_re.search(line)
            if gid and gsym:
                mapping[gid.group(1)] = gsym.group(1)
    log.info(f"  Mapped {len(mapping)} Ensembl gene IDs")
    return mapping


def load_expression() -> pd.DataFrame:
    """
    Returns long-format DataFrame:
      sample, gene, cpm
    """
    os.makedirs(os.path.dirname(COUNTS_TXT), exist_ok=True)
    if not os.path.exists(COUNTS_TXT):
        fetch_s3(f"{S3_RNA}/counts/raw_counts.txt", COUNTS_TXT)

    counts = pd.read_csv(COUNTS_TXT, sep="\t", index_col=0)
    # Drop STAR summary rows (N_unmapped etc.)
    counts = counts[~counts.index.str.startswith("N_")]
    # Strip Ensembl version suffix (ENSMUSG00000051951.2 → ENSMUSG00000051951)
    counts.index = counts.index.str.split(".").str[0]

    # Build Ensembl → symbol mapping
    ensembl2sym = _build_ensembl_to_symbol(GTF_GZ)
    counts.index = counts.index.map(lambda x: ensembl2sym.get(x, x))
    counts = counts[~counts.index.duplicated(keep="first")]

    # Compute CPM per sample
    lib_sizes = counts.sum(axis=0)
    cpm = counts.divide(lib_sizes, axis=1) * 1e6

    # Pivot to long format, keep only tumor samples present in counts
    tumor_cols = [c for c in TUMOR_SAMPLES if c in cpm.columns]
    missing    = set(TUMOR_SAMPLES) - set(tumor_cols)
    if missing:
        log.warning(f"Expression missing for samples: {missing} — will use NaN")

    cpm.index.name = "gene"
    long = (cpm[tumor_cols]
            .reset_index()
            .melt(id_vars="gene", var_name="sample", value_name="cpm"))
    long["log2cpm"] = np.log2(long["cpm"] + 1)
    log.info(f"Expression loaded: {len(tumor_cols)} samples, {cpm.shape[0]} genes")
    return long


# ── 4. Load netMHCpan (optional) ──────────────────────────────────────────────

def load_mhc(mhc_path: str) -> pd.DataFrame:
    """
    Parse merged netMHCpan output. Expected columns: Peptide, MHC, %Rank_EL.
    Returns best (lowest) %Rank_EL per peptide across all alleles.
    """
    df = pd.read_csv(mhc_path, sep="\t", comment="#")
    df.columns = df.columns.str.strip()

    # Normalise column names across netMHCpan versions
    rename = {}
    for col in df.columns:
        lc = col.lower()
        if "peptide" in lc:
            rename[col] = "peptide"
        elif "rank_el" in lc or "rank" in lc:
            rename[col] = "rank_el"
        elif "mhc" in lc or "allele" in lc:
            rename[col] = "MHC"
    df = df.rename(columns=rename)

    if "rank_el" not in df.columns:
        raise ValueError(f"Could not find %Rank_EL column in {mhc_path}. "
                         f"Columns found: {list(df.columns)}")

    # Best rank across alleles per peptide
    best = (df.groupby("peptide")["rank_el"]
              .min()
              .reset_index()
              .rename(columns={"rank_el": "best_rank_el"}))
    log.info(f"netMHCpan loaded: {len(best)} unique peptides with binding predictions")
    return best


# ── 5. Merge & score ──────────────────────────────────────────────────────────

def build_ranked_table(peptides: pd.DataFrame,
                       clonality: pd.DataFrame,
                       expression: pd.DataFrame,
                       mhc: pd.DataFrame | None) -> pd.DataFrame:

    df = peptides.copy()

    # ── Clonality: join by mut_id ─────────────────────────────────────────────
    prev_cols = [c for c in clonality.columns if c.endswith("_cellular_prevalence")]
    clone_cols = ["mut_id", "cluster_id", "mean_prevalence",
                  "max_prevalence", "n_samples_present"] + prev_cols
    df = df.merge(clonality[clone_cols], on="mut_id", how="left")

    # For mutations not in PyClone (low-confidence filtered out), set prevalence = 0
    df["mean_prevalence"] = df["mean_prevalence"].fillna(0)
    df["max_prevalence"]  = df["max_prevalence"].fillna(0)

    # Per-row: use sample-specific cellular prevalence if available
    prev_cols = [c for c in clonality.columns if c.endswith("_cellular_prevalence")]
    def _sample_prev(row):
        col = f"{row['sample']}_cellular_prevalence"
        val = row.get(col, np.nan)
        if pd.isna(val):
            return row["mean_prevalence"]  # already 0-filled for non-PyClone mutations
        return val
    df["sample_prevalence"] = df.apply(_sample_prev, axis=1).fillna(0)

    # ── Expression: join by (sample, gene) ───────────────────────────────────
    df = df.merge(expression[["sample", "gene", "log2cpm"]],
                  on=["sample", "gene"], how="left")

    # ── MHC binding: join by peptide ─────────────────────────────────────────
    if mhc is not None:
        df = df.merge(mhc, on="peptide", how="left")
        df["mhc_score"] = 1 - pct_rank(df["best_rank_el"])  # invert: lower rank_el → higher score
    else:
        df["best_rank_el"] = np.nan
        df["mhc_score"]    = np.nan

    # ── Percentile ranks ─────────────────────────────────────────────────────
    df["clone_pct"] = pct_rank(df["sample_prevalence"])
    df["expr_pct"]  = pct_rank(df["log2cpm"])

    # ── Composite score ───────────────────────────────────────────────────────
    if mhc is not None:
        df["composite_score"] = (
            W_CLONE * df["clone_pct"] +
            W_MHC   * df["mhc_score"] +
            W_EXPR  * df["expr_pct"]
        )
    else:
        # Redistribute MHC weight equally when binding data unavailable
        w_c = W_CLONE / (W_CLONE + W_EXPR)
        w_e = W_EXPR  / (W_CLONE + W_EXPR)
        df["composite_score"] = w_c * df["clone_pct"] + w_e * df["expr_pct"]
        log.warning("netMHCpan data not provided — composite score uses clonality + expression only")

    df["rank"] = df["composite_score"].rank(ascending=False, method="min").astype("Int64")
    df = df.sort_values("rank")

    return df


# ── 6. Output ─────────────────────────────────────────────────────────────────

def write_outputs(df: pd.DataFrame, mhc_available: bool):
    os.makedirs(OUT_DIR, exist_ok=True)

    # Select output columns
    base_cols = [
        "rank", "sample", "gene", "hgvsp", "effect", "peptide", "length",
        "is_frameshift", "mut_id",
        "cluster_id", "sample_prevalence", "mean_prevalence", "n_samples_present",
        "log2cpm", "composite_score", "clone_pct", "expr_pct",
    ]
    if mhc_available:
        base_cols += ["best_rank_el", "mhc_score"]

    out_cols = [c for c in base_cols if c in df.columns]
    out = df[out_cols].copy()

    # Round floats
    for col in ["sample_prevalence", "mean_prevalence", "log2cpm",
                "composite_score", "clone_pct", "expr_pct", "mhc_score", "best_rank_el"]:
        if col in out.columns:
            out[col] = out[col].round(4)

    ranked_tsv = f"{OUT_DIR}/ranked_neoantigens.tsv"
    out.to_csv(ranked_tsv, sep="\t", index=False)
    os.system(f"aws s3 cp {ranked_tsv} {S3_FINAL}/ranked/ranked_neoantigens.tsv --quiet")
    log.info(f"Full ranked table: {len(out)} rows → {ranked_tsv}")

    # Summary: top peptide per mutation per sample
    summary = (out.sort_values("rank")
                  .drop_duplicates(subset=["sample", "mut_id"])
                  [out_cols])
    summary_tsv = f"{OUT_DIR}/ranked_neoantigens_summary.tsv"
    summary.to_csv(summary_tsv, sep="\t", index=False)
    os.system(f"aws s3 cp {summary_tsv} {S3_FINAL}/ranked/ranked_neoantigens_summary.tsv --quiet")
    log.info(f"Summary (top peptide/mutation): {len(summary)} rows → {summary_tsv}")

    # Strong binders only (if MHC available: Rank_EL < 2%)
    if mhc_available and "best_rank_el" in out.columns:
        strong = out[out["best_rank_el"] < 2.0].copy()
        strong_tsv = f"{OUT_DIR}/ranked_neoantigens_strong_binders.tsv"
        strong.to_csv(strong_tsv, sep="\t", index=False)
        os.system(f"aws s3 cp {strong_tsv} {S3_FINAL}/ranked/ranked_neoantigens_strong_binders.tsv --quiet")
        log.info(f"Strong binders (<2% Rank_EL): {len(strong)} peptides → {strong_tsv}")

    # Print top 20
    print("\n" + "="*90)
    print("TOP 20 RANKED NEOANTIGENS")
    print("="*90)
    display_cols = ["rank", "sample", "gene", "hgvsp", "peptide", "length",
                    "sample_prevalence", "log2cpm", "composite_score"]
    if mhc_available:
        display_cols.insert(-1, "best_rank_el")
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.width", 200)
    print(out[display_cols].head(20).to_string(index=False))
    print("="*90 + "\n")

    return ranked_tsv


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mhc", default=None,
                        help="Path to merged netMHCpan TSV (optional)")
    parser.add_argument("--out", default=None,
                        help="Override output directory")
    args = parser.parse_args()

    global OUT_DIR
    if args.out:
        OUT_DIR = args.out

    log.info("="*60)
    log.info("Step 10 — Neoantigen Ranking")
    log.info(f"  MHC input : {args.mhc or 'not provided (clonality+expression only)'}")
    log.info(f"  Weights   : clonality={W_CLONE}, MHC={W_MHC}, expression={W_EXPR}")
    log.info("="*60)

    peptides   = load_peptides()
    clonality  = load_clonality()
    expression = load_expression()
    mhc        = load_mhc(args.mhc) if args.mhc else None

    df = build_ranked_table(peptides, clonality, expression, mhc)
    write_outputs(df, mhc_available=(mhc is not None))

    log.info("Step 10 complete.")


if __name__ == "__main__":
    main()
