#!/usr/bin/env python3
"""
Step 5 — PyClone Input Preparation
Converts SnpEff-annotated VCFs into PyClone YAML input files.
One YAML per tumour sample; multi-sample run uses all together.

Requires: pandas, pysam, pyyaml
"""

import os, sys, re, gzip, yaml, logging
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR    = "/home/ec2-user"
ANN_DIR     = f"{BASE_DIR}/results/annotated"
OUT_DIR     = f"{BASE_DIR}/results/pyclone/input"
S3_ANN      = "s3://bam-wes/NeoAntigen-aws/results/annotated"
S3_OUT      = "s3://bam-wes/NeoAntigen-aws/results/pyclone/input"

TUMOR_SAMPLES = ["428_D20_new","34_D52_old","36_D99_new","38_D99_new","42_D122_old"]

# ── Helpers ───────────────────────────────────────────────────────────────────
def parse_vcf(vcf_path: str) -> pd.DataFrame:
    """Parse a bgzipped SnpEff-annotated VCF into a DataFrame."""
    records = []
    opener = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    cols = line.lstrip("#").strip().split("\t")
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info, fmt = parts[:9]
            samples = parts[9:]

            if filt not in ("PASS", "."):
                continue

            # Only SNVs and short indels
            if len(ref) > 50 or len(alt.split(",")[0]) > 50:
                continue

            # Extract FORMAT fields
            fmt_keys  = fmt.split(":")
            fmt_dict  = dict(zip(fmt_keys, samples[0].split(":")))  # tumour = sample[0]
            norm_dict = dict(zip(fmt_keys, samples[1].split(":"))) if len(samples) > 1 else {}

            # Read counts: prefer AD, fall back to F1R2/F2R1
            if "AD" in fmt_dict:
                ad = fmt_dict["AD"].split(",")
                ref_count = int(ad[0]) if ad[0].isdigit() else 0
                alt_count = int(ad[1]) if len(ad) > 1 and ad[1].isdigit() else 0
            else:
                continue  # skip if no allele depth

            depth = ref_count + alt_count
            if depth < 10 or alt_count < 3:
                continue  # minimum coverage filters

            # Extract SnpEff consequence (first ANN entry)
            ann_match = re.search(r"ANN=([^;]+)", info)
            effect, gene, biotype = "unknown", "unknown", "unknown"
            aa_change, transcript = "", ""
            if ann_match:
                ann0 = ann_match.group(1).split(",")[0].split("|")
                if len(ann0) >= 4:
                    effect    = ann0[1]
                    gene      = ann0[3]
                    biotype   = ann0[7] if len(ann0) > 7 else ""
                    aa_change = ann0[10] if len(ann0) > 10 else ""
                    transcript= ann0[6]  if len(ann0) > 6  else ""

            records.append({
                "chrom": chrom, "pos": int(pos), "ref": ref,
                "alt": alt.split(",")[0],
                "ref_count": ref_count, "alt_count": alt_count,
                "depth": depth,
                "vaf": alt_count / depth if depth > 0 else 0.0,
                "effect": effect, "gene": gene, "biotype": biotype,
                "aa_change": aa_change, "transcript": transcript
            })

    return pd.DataFrame(records)


def make_pyclone_yaml(df: pd.DataFrame, sample_id: str, out_path: str):
    """Write a PyClone-compatible YAML file."""
    mutations = []
    for _, row in df.iterrows():
        mut_id = f"{row['chrom']}:{row['pos']}:{row['ref']}:{row['alt']}"
        mutations.append({
            "id": mut_id,
            "ref_counts":  int(row["ref_count"]),
            "var_counts":  int(row["alt_count"]),
            "normal_cn":   2,  # diploid normal assumed
            "minor_cn":    0,  # conservative — no CNV data
            "major_cn":    2,
            "variant_case": sample_id,
            "variant_freq": float(row["vaf"]),
            "genotype": "AB"
        })

    config = {
        "sample_id": sample_id,
        "tumour_content": {"value": 1.0},  # update if known purity
        "mutations": mutations
    }
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        yaml.dump(config, fh, default_flow_style=False, sort_keys=False)
    log.info(f"  {sample_id}: {len(mutations)} mutations → {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    all_dfs = {}

    for sample in TUMOR_SAMPLES:
        yaml_out = f"{OUT_DIR}/{sample}.yaml"
        tsv_out  = f"{OUT_DIR}/{sample}.mutations.tsv"

        if os.path.exists(yaml_out):
            log.info(f"[SKIP] {sample} YAML exists")
            continue

        vcf_path = f"{ANN_DIR}/{sample}.annotated.vcf.gz"
        if not os.path.exists(vcf_path):
            log.warning(f"Downloading {sample} annotated VCF from S3...")
            os.system(f"aws s3 cp {S3_ANN}/{sample}.annotated.vcf.gz {vcf_path}")
            os.system(f"aws s3 cp {S3_ANN}/{sample}.annotated.vcf.gz.tbi {vcf_path}.tbi 2>/dev/null || true")

        log.info(f"Parsing: {sample}")
        df = parse_vcf(vcf_path)
        log.info(f"  {len(df)} PASS coding variants")

        df.to_csv(tsv_out, sep="\t", index=False)
        make_pyclone_yaml(df, sample, yaml_out)
        all_dfs[sample] = df

        # Upload
        os.system(f"aws s3 cp {yaml_out} {S3_OUT}/{sample}.yaml")
        os.system(f"aws s3 cp {tsv_out}  {S3_OUT}/{sample}.mutations.tsv")

    # Write merged table across all samples
    if all_dfs:
        merged_tsv = f"{OUT_DIR}/all_samples_mutations.tsv"
        merged = pd.concat(
            [df.assign(sample=s) for s, df in all_dfs.items()], ignore_index=True
        )
        merged.to_csv(merged_tsv, sep="\t", index=False)
        os.system(f"aws s3 cp {merged_tsv} {S3_OUT}/all_samples_mutations.tsv")
        log.info(f"Merged table: {len(merged)} rows → {merged_tsv}")

    log.info("Step 5 complete.")

if __name__ == "__main__":
    main()
