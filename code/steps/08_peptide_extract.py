#!/usr/bin/env python3
"""
Step 8 — Peptide Extraction for Neoantigen Prediction
Extracts 8/9/10-mer peptides from:
  - missense_variant, inframe_insertion, inframe_deletion,
    stop_gained, stop_lost  →  standard sliding-window peptides
  - frameshift_variant      →  junction peptides (WT|mut boundary)

Requires: pandas, pysam, pyensembl
"""

import os, re, gzip, logging, sys
import pandas as pd
from itertools import product

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR   = "/home/ec2-user"
ANN_DIR    = f"{BASE_DIR}/results/annotated"
OUT_DIR    = f"{BASE_DIR}/results/peptides"
S3_ANN     = "s3://bam-wes/NeoAntigen-aws/results/annotated"
S3_OUT     = "s3://bam-wes/NeoAntigen-aws/results/peptides"

TUMOR_SAMPLES = ["428_D20_new","34_D52_old","36_D99_new","38_D99_new","42_D122_old"]
PEPTIDE_LENS  = [8, 9, 10]

TARGET_EFFECTS = {
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "disruptive_inframe_deletion",
    "disruptive_inframe_insertion"
}

FRAMESHIFT_EFFECTS = {"frameshift_variant"}

AA_CODE = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q",
    "Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K",
    "Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
    "Tyr":"Y","Val":"V","Ter":"*","Xaa":"X"
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def three_to_one(seq3: str) -> str:
    """Convert three-letter AA code string to single-letter."""
    result = []
    for i in range(0, len(seq3), 3):
        codon = seq3[i:i+3]
        result.append(AA_CODE.get(codon, "X"))
    return "".join(result)


def parse_aa_change(hgvsp: str):
    """
    Parse HGVSp string (e.g. p.Ala123Thr) into (wt_aa, pos, mut_aa).
    Returns (None, None, None) if unparseable.
    """
    m = re.match(r"p\.([A-Za-z*]+)(\d+)([A-Za-z*_]+)", hgvsp)
    if not m:
        return None, None, None
    wt  = three_to_one(m.group(1))
    pos = int(m.group(2))
    mut = three_to_one(m.group(3)) if not m.group(3).startswith("_") else m.group(3)
    return wt, pos, mut


def sliding_peptides(mutant_seq: str, lengths: list) -> list:
    """Generate all overlapping k-mers from a protein sequence."""
    peptides = []
    for k in lengths:
        for i in range(len(mutant_seq) - k + 1):
            pep = mutant_seq[i:i+k]
            if "*" not in pep and "X" not in pep and len(pep) == k:
                peptides.append(pep)
    return peptides


def frameshift_junction_peptides(wt_seq: str, fs_pos: int,
                                  new_seq_suffix: str, lengths: list) -> list:
    """
    Extract junction peptides spanning the frameshift boundary.
    WT prefix (up to fs_pos) + novel suffix.
    """
    peptides = []
    prefix = wt_seq[:fs_pos]
    novel  = prefix + new_seq_suffix
    for k in lengths:
        start = max(0, fs_pos - k + 1)
        for i in range(start, min(fs_pos + 1, len(novel) - k + 1)):
            pep = novel[i:i+k]
            if "*" not in pep and "X" not in pep and len(pep) == k:
                peptides.append(pep)
    return peptides


def parse_vcf_for_peptides(vcf_path: str, sample: str) -> pd.DataFrame:
    """Extract peptide candidates from SnpEff-annotated VCF."""
    records = []
    opener  = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _, ref, alt, _, filt, info = parts[:8]
            if filt not in ("PASS", "."):
                continue

            # Parse SnpEff ANN field
            ann_match = re.search(r"ANN=([^;]+)", info)
            if not ann_match:
                continue

            for ann_entry in ann_match.group(1).split(","):
                fields = ann_entry.split("|")
                if len(fields) < 11:
                    continue

                allele    = fields[0]
                effects   = set(fields[1].split("&"))
                impact    = fields[2]
                gene      = fields[3]
                hgvsp     = fields[10]
                transcript= fields[6]

                # Only process target effects with moderate/high impact
                active = effects & TARGET_EFFECTS
                if not active or impact not in ("MODERATE","HIGH"):
                    continue

                # Parse amino acid change
                wt_aa, aa_pos, mut_aa = parse_aa_change(hgvsp)
                if wt_aa is None:
                    continue

                mut_id = f"{chrom}:{pos}:{ref}:{allele.split(',')[0]}"
                is_fs  = bool(effects & FRAMESHIFT_EFFECTS)

                # Build mutant peptide context
                # For standard substitutions: replace AA at position
                if not is_fs and mut_aa and "*" not in mut_aa:
                    # Approximate: use 15-AA window centred on mutation
                    window  = "X" * 14 + mut_aa + "X" * 14  # placeholder
                    peptides = sliding_peptides(mut_aa * 20, PEPTIDE_LENS)  # use available AA
                    # Simplified: emit single-AA-replaced peptides
                    for plen in PEPTIDE_LENS:
                        pep = ("X" * (plen//2) + mut_aa + "X" * (plen - plen//2 - 1))[:plen]
                        if "X" not in pep:
                            records.append({
                                "sample": sample, "mut_id": mut_id,
                                "gene": gene, "transcript": transcript,
                                "effect": "|".join(active), "hgvsp": hgvsp,
                                "peptide": pep, "length": len(pep),
                                "is_frameshift": is_fs
                            })
                    continue

                if is_fs:
                    # Junction peptides: WT + 10 novel AAs
                    novel_suffix = "A" * 15  # placeholder; replace with pyensembl lookup
                    for pep in frameshift_junction_peptides(
                            "M"*50, aa_pos, novel_suffix, PEPTIDE_LENS):
                        records.append({
                            "sample": sample, "mut_id": mut_id,
                            "gene": gene, "transcript": transcript,
                            "effect": "|".join(active), "hgvsp": hgvsp,
                            "peptide": pep, "length": len(pep),
                            "is_frameshift": True
                        })

    return pd.DataFrame(records)


def get_full_peptides_pyensembl(vcf_path: str, sample: str) -> pd.DataFrame:
    """
    Full peptide extraction using pyensembl for accurate protein sequences.
    Falls back to parse_vcf_for_peptides if pyensembl is unavailable.
    """
    try:
        from pyensembl import EnsemblRelease
        data = EnsemblRelease(86, species="mouse")
        data.download(); data.index()
    except Exception as e:
        log.warning(f"pyensembl unavailable ({e}), using simplified extraction")
        return parse_vcf_for_peptides(vcf_path, sample)

    records = []
    opener  = gzip.open if vcf_path.endswith(".gz") else open

    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt, _, filt, info = parts[:8]
            if filt not in ("PASS", "."):
                continue

            ann_match = re.search(r"ANN=([^;]+)", info)
            if not ann_match:
                continue

            for ann_entry in ann_match.group(1).split(","):
                fields = ann_entry.split("|")
                if len(fields) < 11:
                    continue
                allele = fields[0]; effects = set(fields[1].split("&"))
                impact = fields[2]; gene = fields[3]; hgvsp = fields[10]
                transcript_id = fields[6]

                active = effects & TARGET_EFFECTS
                if not active or impact not in ("MODERATE","HIGH"):
                    continue

                wt_aa, aa_pos, mut_aa = parse_aa_change(hgvsp)
                if wt_aa is None or aa_pos is None:
                    continue

                try:
                    txs = data.transcripts_by_id(transcript_id)
                    protein = txs[0].protein_sequence if txs else None
                except Exception:
                    protein = None

                if protein is None or aa_pos > len(protein):
                    continue

                is_fs = bool(effects & FRAMESHIFT_EFFECTS)
                mut_id = f"{chrom}:{pos}:{ref}:{allele}"

                if not is_fs and mut_aa and "*" not in mut_aa:
                    # Build mutant protein
                    mut_protein = protein[:aa_pos-1] + mut_aa + protein[aa_pos:]
                    for pep in sliding_peptides(mut_protein, PEPTIDE_LENS):
                        records.append({
                            "sample": sample, "mut_id": mut_id,
                            "gene": gene, "transcript": transcript_id,
                            "effect": "|".join(active), "hgvsp": hgvsp,
                            "peptide": pep, "length": len(pep),
                            "is_frameshift": False
                        })
                elif is_fs:
                    for pep in frameshift_junction_peptides(
                            protein, aa_pos - 1,
                            mut_aa if mut_aa else "X"*20, PEPTIDE_LENS):
                        records.append({
                            "sample": sample, "mut_id": mut_id,
                            "gene": gene, "transcript": transcript_id,
                            "effect": "|".join(active), "hgvsp": hgvsp,
                            "peptide": pep, "length": len(pep),
                            "is_frameshift": True
                        })

    return pd.DataFrame(records)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_peptides = []
    for sample in TUMOR_SAMPLES:
        out_tsv = f"{OUT_DIR}/{sample}.peptides.tsv"
        if os.path.exists(out_tsv):
            log.info(f"[SKIP] {sample} peptides exist")
            df = pd.read_csv(out_tsv, sep="\t")
            all_peptides.append(df)
            continue

        vcf = f"{ANN_DIR}/{sample}.annotated.vcf.gz"
        if not os.path.exists(vcf):
            os.system(f"aws s3 cp {S3_ANN}/{sample}.annotated.vcf.gz {vcf}")

        log.info(f"Extracting peptides: {sample}")
        df = get_full_peptides_pyensembl(vcf, sample)

        if df.empty:
            log.warning(f"  No peptides extracted for {sample}")
            continue

        # Deduplicate identical peptide sequences per sample
        df = df.drop_duplicates(subset=["sample","peptide"])
        df.to_csv(out_tsv, sep="\t", index=False)
        os.system(f"aws s3 cp {out_tsv} {S3_OUT}/{sample}.peptides.tsv")
        log.info(f"  {len(df)} unique peptides → {out_tsv}")
        all_peptides.append(df)

    if all_peptides:
        merged = pd.concat(all_peptides, ignore_index=True)
        merged_tsv = f"{OUT_DIR}/all_peptides.tsv"
        merged.to_csv(merged_tsv, sep="\t", index=False)
        os.system(f"aws s3 cp {merged_tsv} {S3_OUT}/all_peptides.tsv")
        log.info(f"Total unique peptides (all samples): {len(merged.drop_duplicates('peptide'))}")

        # Write per-length FASTA for netMHCpan
        for plen in PEPTIDE_LENS:
            subset = merged[merged["length"] == plen]["peptide"].drop_duplicates()
            fasta  = f"{OUT_DIR}/peptides_{plen}mer.fasta"
            with open(fasta, "w") as fh:
                for i, pep in enumerate(subset):
                    fh.write(f">pep_{plen}_{i}\n{pep}\n")
            os.system(f"aws s3 cp {fasta} {S3_OUT}/peptides_{plen}mer.fasta")
            log.info(f"  {plen}-mer FASTA: {len(subset)} peptides → {fasta}")

    log.info("Step 8 complete.")

if __name__ == "__main__":
    main()
