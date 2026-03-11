#!/usr/bin/env python3
"""
Step 8b — Fusion Peptide Extraction from Arriba Fusion Calls
Extracts 8/9/10-mer junction peptides from Arriba fusion predictions.

Arriba output columns used:
  #gene1, gene2         — fusion partner gene names
  confidence            — high/medium/low
  reading_frame         — in-frame / out-of-frame / .
  peptide_sequence      — fusion protein AA sequence; | marks junction, ___ is truncation

Outputs:
  results/peptides/fusions_all.tsv     — per-fusion metadata + peptides
  results/peptides/peptides_{8,9,10}mer.fasta  — APPENDED (not overwritten)
"""

import os, re, sys, logging
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ── Config ─────────────────────────────────────────────────────────────────────
BASE_DIR  = "/home/ec2-user/results"
OUT_DIR   = f"{BASE_DIR}/peptides"
S3_ARRIBA = "s3://bam-wes/NeoAntigen-aws/data/FusionCalling/arriba_files"
S3_OUT    = "s3://bam-wes/NeoAntigen-aws/results/peptides"
LOCAL_DIR = f"{BASE_DIR}/../tmp/arriba"
PEPTIDE_LENS = [8, 9, 10]
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOCAL_DIR, exist_ok=True)

# Map Arriba file suffixes → WES/RNASeq sample IDs
SAMPLE_MAP = {
    "D20_new":    "428_D20_new",
    "D52_old":    "34_D52_old",
    "D99_new_36": "36_D99_new",
    "D99_new_38": "38_D99_new",
    "D122_old":   "42_D122_old",
    "D88_old":    "D88_old",    # RNASeq-only timepoint
    "D99_old":    "D99_old",    # RNASeq-only timepoint
    "D109_new":   "D109_new",   # RNASeq-only timepoint
}

# ── Helpers ────────────────────────────────────────────────────────────────────
def sliding_peptides(seq: str, lengths: list, window_around_junction: int) -> list:
    """Return k-mers that span the junction (marked by | position)."""
    jpos = seq.find("|")
    if jpos < 0:
        return []
    clean = seq.replace("|", "").replace("_", "")
    peptides = []
    for k in lengths:
        # Only report k-mers that overlap the junction
        start = max(0, jpos - k + 1)
        end   = min(jpos + 1, len(clean) - k + 1)
        for i in range(start, end):
            pep = clean[i:i+k]
            if len(pep) == k and "*" not in pep and "X" not in pep and "." not in pep:
                peptides.append(pep)
    return list(set(peptides))


# ── Main ───────────────────────────────────────────────────────────────────────
def process_arriba(arriba_file: str, sample_id: str) -> list:
    """Parse one Arriba file and return list of peptide records."""
    records = []
    df = pd.read_csv(arriba_file, sep="\t", comment=None)
    # Rename #gene1 → gene1
    df.columns = [c.lstrip("#") for c in df.columns]

    # Filter
    df = df[df["confidence"].isin(["high", "medium"])]
    log.info(f"  {sample_id}: {len(df)} high/medium confidence fusions")

    for _, row in df.iterrows():
        pep_seq = str(row.get("peptide_sequence", "") or "")
        if not pep_seq or pep_seq == "." or "|" not in pep_seq:
            continue

        # Determine mutation type: in-frame = missense-like, out-of-frame = frameshift
        rf = str(row.get("reading_frame", "."))
        if rf == ".":
            rf = "unknown"

        gene1 = str(row.get("gene1", "")).split("(")[0]
        gene2 = str(row.get("gene2", "")).split("(")[0]
        fusion_id = f"{gene1}::{gene2}"
        bp1  = str(row.get("breakpoint1", ""))
        bp2  = str(row.get("breakpoint2", ""))
        conf = str(row.get("confidence", ""))

        peps = sliding_peptides(pep_seq, PEPTIDE_LENS, window_around_junction=15)
        for pep in peps:
            records.append({
                "sample":       sample_id,
                "mut_id":       f"{bp1}::{bp2}",
                "gene":         fusion_id,
                "transcript":   "",
                "effect":       f"fusion_{rf}",
                "hgvsp":        "",
                "peptide":      pep,
                "length":       len(pep),
                "is_frameshift": rf == "out-of-frame",
                "confidence":   conf,
            })
    return records


def main():
    all_records = []

    for suffix, sample_id in SAMPLE_MAP.items():
        local_f = f"{LOCAL_DIR}/Arriba_{suffix}.txt"
        s3_path = f"{S3_ARRIBA}/Arriba_{suffix}.txt"

        if not os.path.exists(local_f):
            ret = os.system(f"aws s3 cp '{s3_path}' '{local_f}' 2>/dev/null")
            if ret != 0:
                log.warning(f"  [SKIP] Not found on S3: {suffix}")
                continue

        log.info(f"Processing: Arriba_{suffix}.txt → {sample_id}")
        recs = process_arriba(local_f, sample_id)
        log.info(f"  {sample_id}: {len(recs)} fusion junction peptides")
        all_records.extend(recs)

    if not all_records:
        log.warning("No fusion peptides extracted")
        return

    df = pd.DataFrame(all_records).drop_duplicates(subset=["sample", "peptide"])
    out_tsv = f"{OUT_DIR}/fusions_all.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    os.system(f"aws s3 cp '{out_tsv}' '{S3_OUT}/fusions_all.tsv'")
    log.info(f"Fusion metadata: {len(df)} rows → {out_tsv}")

    # Append junction peptides to per-length FASTAs (create if missing)
    for plen in PEPTIDE_LENS:
        subset = df[df["length"] == plen]["peptide"].drop_duplicates()
        fasta = f"{OUT_DIR}/peptides_{plen}mer.fasta"
        mode  = "a" if os.path.exists(fasta) else "w"
        n_written = 0
        with open(fasta, mode) as fh:
            for i, pep in enumerate(subset):
                fh.write(f">fusion_{plen}_{i}\n{pep}\n")
                n_written += 1
        os.system(f"aws s3 cp '{fasta}' '{S3_OUT}/peptides_{plen}mer.fasta'")
        log.info(f"  {plen}-mer FASTA: {n_written} fusion peptides appended → {fasta}")

    log.info(f"Step 8b complete. Total fusion peptides: {len(df)}")


if __name__ == "__main__":
    main()
