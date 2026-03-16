#!/usr/bin/env python3
"""
Step 14 — Build fully annotated final peptides table.

Joins all pipeline outputs into a single comprehensive table:

  ranked table          → rank, composite_score, MHC predictions, expression
  clonality_final.tsv   → mean_ccf, max_ccf across samples
  loci_all.tsv          → per-sample CCF, CCF_std, cluster_prob (pivoted)
  per-sample inputs     → per-sample VAF, ref/alt counts, major/minor CN, purity

Output columns (one row per peptide):
  --- Identity ---
  rank, gene, sample, transcript, mut_id, hgvsp, peptide, length
  variant_type, is_frameshift, is_fusion

  --- MHC Binding ---
  best_allele, binding_class, rank_el_H2Db, mhc_score

  --- Expression ---
  log2cpm, expr_pct, rna_validated

  --- Clonality (per-sample CCF) ---
  tumor_purity          ← VAF-based per-sample purity used in PyClone-VI
  ccf_mean, ccf_max     ← mean/max across samples (mutations present in >1 get averaged)
  ccf_std               ← std of CCF across samples where mutation detected
  clonality_class       ← 'clonal'/'subclonal'/'minor' based on max_ccf threshold
  cluster_id            ← PyClone-VI cluster assignment
  n_samples_pyclone     ← number of samples where PyClone-VI ran for this mutation
  ccf_443_D21_new, ccf_428_D20_new, ccf_34_D52_old,
  ccf_36_D99_new, ccf_38_D99_new, ccf_42_D122_old   ← per-sample CCF (NaN if absent)
  ccf_std_443_D21_new, ... (per-sample CCF std)

  --- Variant Read Counts (per sample where present) ---
  ref_counts_{sample}, alt_counts_{sample}, vaf_{sample}
  major_cn_{sample}, minor_cn_{sample}, total_cn_{sample}

  --- Cohort Prevalence ---
  sample_prevalence     ← fraction of WES samples with this mutation
  n_samples_present     ← number of WES samples where this mutation was called

  --- Score ---
  composite_score

Output: results/ranked/final_peptides_{timestamp}.tsv
         s3://neoantigen2026-rerun/final-results/ranked/final_peptides_{timestamp}.tsv
"""

import os, sys, logging
from datetime import datetime
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

BASE_DIR    = "/home/ec2-user"
RESULTS_DIR = os.environ.get("RESULTS_DIR", f"{BASE_DIR}/results/res_20260311_225555")
CNVKIT_DIR  = f"{RESULTS_DIR}/pyclone/cnvkit"
RANKED_DIR  = f"{RESULTS_DIR}/ranked"
S3_FINAL    = "s3://neoantigen2026-rerun/final-results"

SAMPLE_PURITY = {
    "443_D21_new": 0.18,
    "428_D20_new": 0.22,
    "34_D52_old":  0.38,
    "36_D99_new":  0.22,
    "38_D99_new":  0.20,
    "42_D122_old": 1.00,
}
SAMPLES = sorted(SAMPLE_PURITY.keys())

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")


def load_ranked() -> pd.DataFrame:
    path = f"{RANKED_DIR}/neoantigens_final.tsv"
    df = pd.read_csv(path, sep="\t")
    log.info(f"Loaded ranked table: {len(df)} rows")
    return df


def load_clonality() -> pd.DataFrame:
    path = f"{CNVKIT_DIR}/clonality_final.tsv"
    clon = pd.read_csv(path, sep="\t")
    clon = clon.rename(columns={
        "mutation_id": "mut_id",
        "mean_ccf":    "ccf_mean",
        "max_ccf":     "ccf_max",
        "cluster_id":  "cluster_id_pyclone",
    })
    log.info(f"Loaded clonality: {len(clon)} mutations")
    return clon


def load_per_sample_ccf() -> pd.DataFrame:
    """
    Pivot loci_all.tsv to get per-sample CCF columns.
    Returns wide DataFrame indexed by mutation_id.
    Note: summary stats (mean, max, std) come from clonality_final.tsv — not duplicated here.
    """
    path = f"{CNVKIT_DIR}/loci_all.tsv"
    loci = pd.read_csv(path, sep="\t")

    # Count samples per mutation
    n_samples = loci.groupby("mutation_id")["sample_id"].nunique().rename("n_samples_pyclone")

    # Std of CCF across samples (std is not in clonality_final, add here)
    ccf_std = loci.groupby("mutation_id")["cellular_prevalence"].std().rename("ccf_std")

    # Per-sample CCF
    ccf_wide = loci.pivot_table(
        index="mutation_id",
        columns="sample_id",
        values="cellular_prevalence",
        aggfunc="first",
    )
    ccf_wide.columns = [f"ccf_{c}" for c in ccf_wide.columns]

    # Per-sample CCF measurement std (model uncertainty per sample)
    std_wide = loci.pivot_table(
        index="mutation_id",
        columns="sample_id",
        values="cellular_prevalence_std",
        aggfunc="first",
    )
    std_wide.columns = [f"ccf_mstd_{c}" for c in std_wide.columns]

    # Per-sample cluster assignment probability
    prob_wide = loci.pivot_table(
        index="mutation_id",
        columns="sample_id",
        values="cluster_assignment_prob",
        aggfunc="first",
    )
    prob_wide.columns = [f"cluster_prob_{c}" for c in prob_wide.columns]

    result = pd.concat([n_samples, ccf_std, ccf_wide, std_wide, prob_wide], axis=1).reset_index()
    result = result.rename(columns={"mutation_id": "mut_id"})
    log.info(f"Per-sample CCF wide table: {len(result)} mutations, "
             f"{len(result.columns)} columns")
    return result


def load_per_sample_counts() -> pd.DataFrame:
    """
    Load per-sample read counts, VAF, and copy number from individual input TSVs.
    Returns wide DataFrame: one row per mutation_id.
    """
    dfs = []
    for sample in SAMPLES:
        path = f"{CNVKIT_DIR}/{sample}_input.tsv"
        if not os.path.exists(path):
            log.warning(f"  No input TSV for {sample}")
            continue
        s = pd.read_csv(path, sep="\t")
        s = s.rename(columns={"mutation_id": "mut_id"})
        s["vaf"] = s["alt_counts"] / (s["ref_counts"] + s["alt_counts"] + 1e-9)
        s["vaf"] = s["vaf"].round(4)
        s["total_cn"] = s["major_cn"] + s["minor_cn"]
        # keep only what we want per sample
        keep = ["mut_id", "ref_counts", "alt_counts", "vaf",
                "major_cn", "minor_cn", "total_cn", "tumour_content"]
        s = s[keep].rename(columns={
            "ref_counts":    f"ref_counts_{sample}",
            "alt_counts":    f"alt_counts_{sample}",
            "vaf":           f"vaf_{sample}",
            "major_cn":      f"major_cn_{sample}",
            "minor_cn":      f"minor_cn_{sample}",
            "total_cn":      f"total_cn_{sample}",
            "tumour_content": f"tumor_purity_{sample}",
        })
        dfs.append(s)
        log.info(f"  Loaded counts for {sample}: {len(s)} mutations")

    if not dfs:
        return pd.DataFrame()

    # Merge all samples on mut_id (outer join — not all mutations present in all samples)
    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="mut_id", how="outer")
    log.info(f"Per-sample counts wide: {len(merged)} mutations")
    return merged


def classify_clonality(max_ccf: float) -> str:
    if pd.isna(max_ccf):
        return "unassigned"
    if max_ccf > 0.8:
        return "clonal"
    if max_ccf >= 0.2:
        return "subclonal"
    return "minor"


def build_purity_column(df: pd.DataFrame) -> pd.Series:
    """Add tumor_purity from the sample column."""
    return df["sample"].map(SAMPLE_PURITY)


def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Reorder for readability."""
    # Define preferred order
    priority = [
        # Identity
        "rank", "composite_score", "gene", "sample", "tumor_purity",
        "transcript", "mut_id", "hgvsp", "peptide", "length",
        "variant_type", "is_frameshift", "is_fusion",
        # MHC
        "best_allele", "binding_class", "rank_el_H2Db", "mhc_score",
        # Expression
        "log2cpm", "expr_pct", "rna_validated",
        # Clonality summary
        "clonality_class", "ccf_mean", "ccf_max", "ccf_std",
        "n_samples_pyclone", "cluster_id",
        # Cohort prevalence
        "sample_prevalence", "n_samples_present",
    ]
    # Per-sample CCF columns (exclude summary columns already in priority list)
    _summary_ccf = {"ccf_mean", "ccf_max", "ccf_std"}
    ccf_cols  = sorted([c for c in df.columns
                        if c.startswith("ccf_") and not c.startswith("ccf_mstd_")
                        and c not in _summary_ccf])
    ccf_std_cols = sorted([c for c in df.columns if c.startswith("ccf_mstd_")])
    prob_cols = sorted([c for c in df.columns if c.startswith("cluster_prob_")])
    # Per-sample counts/CN
    ref_cols   = sorted([c for c in df.columns if c.startswith("ref_counts_")])
    alt_cols   = sorted([c for c in df.columns if c.startswith("alt_counts_")])
    vaf_cols   = sorted([c for c in df.columns if c.startswith("vaf_")])
    maj_cols   = sorted([c for c in df.columns if c.startswith("major_cn_")])
    min_cols   = sorted([c for c in df.columns if c.startswith("minor_cn_")])
    tot_cols   = sorted([c for c in df.columns if c.startswith("total_cn_")])
    pur_cols   = sorted([c for c in df.columns if c.startswith("tumor_purity_")])

    ordered = (priority + ccf_cols + ccf_std_cols + prob_cols +
               ref_cols + alt_cols + vaf_cols +
               maj_cols + min_cols + tot_cols + pur_cols)

    # Add any remaining columns not in ordered list
    remaining = [c for c in df.columns if c not in ordered]
    final_order = [c for c in ordered if c in df.columns] + remaining
    return df[final_order]


def main():
    log.info("Step 14 — Building fully annotated final peptides table")

    # 1. Base ranked table
    ranked = load_ranked()

    # 2. Add per-sample purity from sample column
    ranked["tumor_purity"] = build_purity_column(ranked)

    # 3. Clonality summary (mean_ccf, max_ccf per mutation across samples)
    clon = load_clonality()
    ranked = ranked.merge(clon, on="mut_id", how="left")

    # Prefer cluster_id from clonality file (drop old one)
    if "cluster_id_pyclone" in ranked.columns:
        ranked["cluster_id"] = ranked["cluster_id_pyclone"].combine_first(ranked["cluster_id"])
        ranked.drop(columns=["cluster_id_pyclone"], inplace=True)

    # Add clonality class based on max_ccf
    ranked["clonality_class"] = ranked["ccf_max"].apply(classify_clonality)
    # For fusion peptides, use clone_pct as proxy
    fusion_mask = ranked["is_fusion"] == True
    ranked.loc[fusion_mask & ranked["ccf_max"].isna(), "clonality_class"] = \
        ranked.loc[fusion_mask & ranked["ccf_max"].isna(), "clone_pct"].apply(classify_clonality)

    # 4. Per-sample CCF (pivoted)
    ccf_wide = load_per_sample_ccf()
    if not ccf_wide.empty:
        ranked = ranked.merge(ccf_wide, on="mut_id", how="left")

    # 5. Per-sample read counts + CN
    counts_wide = load_per_sample_counts()
    if not counts_wide.empty:
        ranked = ranked.merge(counts_wide, on="mut_id", how="left")

    # 6. Ensure all per-sample columns exist (fill NaN for absent samples)
    for sample in SAMPLES:
        for prefix in ["ccf_", "ccf_mstd_", "cluster_prob_",
                       "ref_counts_", "alt_counts_", "vaf_",
                       "major_cn_", "minor_cn_", "total_cn_", "tumor_purity_"]:
            col = f"{prefix}{sample}"
            if col not in ranked.columns:
                ranked[col] = np.nan

    # Drop legacy clone_pct (superseded by ccf_mean/ccf_max from PyClone output)
    # Keep it renamed for traceability
    if "clone_pct" in ranked.columns:
        ranked = ranked.rename(columns={"clone_pct": "ccf_legacy"})

    # Round numeric columns
    float_cols = ranked.select_dtypes(include=[float]).columns
    ranked[float_cols] = ranked[float_cols].round(6)

    # 7. Reorder columns
    ranked = reorder_columns(ranked)

    # 8. Save
    out_path = f"{RANKED_DIR}/final_peptides_{TIMESTAMP}.tsv"
    ranked.to_csv(out_path, sep="\t", index=False)
    log.info(f"Saved: {out_path}")
    log.info(f"  Shape: {ranked.shape[0]} rows × {ranked.shape[1]} columns")
    log.info(f"  Columns: {list(ranked.columns)}")

    # Summary stats
    binders = ranked[ranked["binding_class"].isin(["SB", "WB", "PB"])]
    log.info(f"\nSummary:")
    log.info(f"  Total peptides:   {len(ranked)}")
    log.info(f"  Binders (SB+WB+PB): {len(binders)}")
    log.info(f"  Strong binders (SB): {(ranked['binding_class']=='SB').sum()}")
    log.info(f"  Clonal (max_ccf>0.8): {(ranked['clonality_class']=='clonal').sum()}")
    log.info(f"  Clonal binders:   {((ranked['clonality_class']=='clonal') & ranked['binding_class'].isin(['SB','WB','PB'])).sum()}")
    log.info(f"  RNA-validated:    {ranked['rna_validated'].sum()}")
    log.info(f"\n  Purity per sample:")
    for s, p in sorted(SAMPLE_PURITY.items()):
        log.info(f"    {s}: {p:.2f}")
    log.info(f"\nTop 10 (all columns):")
    top = binders.head(10)[[
        "rank", "gene", "sample", "tumor_purity", "peptide", "binding_class",
        "rank_el_H2Db", "log2cpm", "clonality_class", "ccf_mean", "ccf_max",
        "ccf_std", "rna_validated", "composite_score"
    ]]
    print(top.to_string(index=False))

    # 9. Upload to S3
    s3_path = f"{S3_FINAL}/ranked/final_peptides_{TIMESTAMP}.tsv"
    ret = os.system(f"aws s3 cp '{out_path}' '{s3_path}' --quiet")
    if ret == 0:
        log.info(f"\nUploaded to S3: {s3_path}")
    else:
        log.error(f"S3 upload failed (exit {ret})")

    return out_path


if __name__ == "__main__":
    main()
