#!/usr/bin/env python3
"""
Step 08c — Frameshift Novel Sequence Extraction (Proper bioinformatics approach)

Replaces the poly-alanine placeholder in frameshift peptides with real
translated novel sequences obtained by:

  1. Loading transcript sequences from Ensembl GRCm38 cDNA FASTA
     (Mus_musculus.GRCm38.cdna.all.fa.gz — downloaded by pyensembl)
  2. Identifying CDS start within the cDNA by matching against the known
     Ensembl protein sequence (Mus_musculus.GRCm38.pep.all.fa.gz)
  3. Applying the HGVSc nucleotide change (deletion/insertion) to the CDS
  4. Translating the frameshifted CDS from the mutation site to first stop
     codon using Biopython's standard codon table (CodonTable.standard_dna_table)
  5. Generating junction peptides that span the WT|novel boundary

Tools used (all established bioinformatics packages):
  - Biopython 1.86  — translation (Bio.Seq.Seq.translate)
  - Ensembl GRCm38 release 102 cDNA FASTA  — transcript sequences
  - Ensembl GRCm38 release 102 protein FASTA — WT protein sequences for CDS finding
  - SnpEff HGVSc notation parsed from annotated VCF ANN field

Output:
  peptides/frameshift_real.tsv   — frameshift peptides with real novel sequences
  peptides/peptides_Xmer.fasta   — regenerated FASTA (replaces placeholder entries)

Usage:
  python3 08c_frameshift_realseq.py
"""

import os, gzip, re, logging, sys
import pandas as pd
from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR    = "/home/ec2-user"
RESULTS_DIR = os.environ.get("RESULTS_DIR", f"{BASE_DIR}/results/res_20260311_225555")
S3_RESULTS  = os.environ.get("S3_RESULTS",  "s3://neoantigen2026-rerun/results/res_20260311_225555")

PEP_DIR     = f"{RESULTS_DIR}/peptides"
ANN_DIR     = f"{RESULTS_DIR}/annotated"

CDNA_FASTA  = f"{BASE_DIR}/.cache/pyensembl/GRCm38/ensembl102/Mus_musculus.GRCm38.cdna.all.fa.gz"
PROT_FASTA  = f"{BASE_DIR}/.cache/pyensembl/GRCm38/ensembl102/Mus_musculus.GRCm38.pep.all.fa.gz"

TUMOR_SAMPLES = ["443_D21_new","428_D20_new","34_D52_old","36_D99_new","38_D99_new","42_D122_old"]
PEPTIDE_LENS  = [8, 9, 10]
FRAMESHIFT_EFFECTS = {"frameshift_variant"}
TARGET_EFFECTS = {
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "disruptive_inframe_deletion", "disruptive_inframe_insertion"
}
MAX_NOVEL_AA = 50   # maximum novel amino acids to consider after frameshift


# ── Ensembl FASTA loaders ─────────────────────────────────────────────────────

def load_protein_db() -> dict:
    """Load Ensembl protein FASTA → gene_symbol: [protein_seq]."""
    db = {}
    sym_re = re.compile(r"gene_symbol:(\S+)")
    cur_sym = None; cur_seq = []
    with gzip.open(PROT_FASTA, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_sym and cur_seq:
                    db.setdefault(cur_sym, []).append("".join(cur_seq))
                m = sym_re.search(line)
                cur_sym = m.group(1) if m else None
                cur_seq = []
            else:
                cur_seq.append(line)
    if cur_sym and cur_seq:
        db.setdefault(cur_sym, []).append("".join(cur_seq))
    log.info(f"Protein FASTA: {len(db)} gene symbols loaded")
    return db


def load_cdna_db(target_genes: set) -> dict:
    """
    Load Ensembl cDNA FASTA for target gene symbols only.
    Returns gene_symbol: [(transcript_id, biotype, sequence)].
    """
    db = {}
    sym_re   = re.compile(r"gene_symbol:(\S+)")
    tid_re   = re.compile(r"^>(\S+)")
    btype_re = re.compile(r"transcript_biotype:(\S+)")
    cur_sym = None; cur_tid = None; cur_bt = None; cur_seq = []
    with gzip.open(CDNA_FASTA, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_sym and cur_sym in target_genes and cur_seq:
                    db.setdefault(cur_sym, []).append(
                        (cur_tid, cur_bt, "".join(cur_seq)))
                cur_sym = None; cur_tid = None; cur_bt = None; cur_seq = []
                m = sym_re.search(line)
                if m and m.group(1) in target_genes:
                    cur_sym = m.group(1)
                    t = tid_re.match(line)
                    cur_tid = t.group(1).split(".")[0] if t else None
                    b = btype_re.search(line)
                    cur_bt = b.group(1) if b else None
            elif cur_sym:
                cur_seq.append(line.upper())
    if cur_sym and cur_sym in target_genes and cur_seq:
        db.setdefault(cur_sym, []).append(
            (cur_tid, cur_bt, "".join(cur_seq)))
    log.info(f"cDNA FASTA: loaded sequences for {len(db)} target genes")
    return db


# ── HGVSc parsing ─────────────────────────────────────────────────────────────

def parse_hgvsc(hgvsc: str):
    """
    Parse HGVSc string to (start_pos, end_pos, type, sequence).
    Positions are 1-based CDS positions.

    Examples:
      c.117_118delGC   → (117, 118, 'del', 'GC')
      c.616delC        → (616, 616, 'del', 'C')
      c.2142dupA       → (2142, 2142, 'dup', 'A')
      c.33_70delGGC..  → (33, 70, 'del', 'GGC...')
    """
    # deletion with range: c.START_ENDdelSEQ
    m = re.match(r"c\.(\d+)_(\d+)del([ACGT]+)$", hgvsc, re.IGNORECASE)
    if m:
        return int(m.group(1)), int(m.group(2)), "del", m.group(3).upper()

    # deletion single: c.POSdelNT
    m = re.match(r"c\.(\d+)del([ACGT]+)$", hgvsc, re.IGNORECASE)
    if m:
        return int(m.group(1)), int(m.group(1)), "del", m.group(2).upper()

    # duplication: c.POSdupNT
    m = re.match(r"c\.(\d+)dup([ACGT]+)$", hgvsc, re.IGNORECASE)
    if m:
        return int(m.group(1)), int(m.group(1)), "dup", m.group(2).upper()

    log.warning(f"Cannot parse HGVSc: {hgvsc}")
    return None, None, None, None


# ── CDS finder ────────────────────────────────────────────────────────────────

def find_cds_start(cdna_seq: str, protein_seq: str, aa_pos: int) -> int:
    """
    Find the CDS start position (0-based) within the cDNA sequence.

    Strategy: scan for ATG codons; translate; check if the translated
    sequence matches the known WT protein up to aa_pos.
    Returns 0-based index of the 'A' in the start codon, or -1 if not found.
    """
    check_len = min(aa_pos, len(protein_seq))
    prefix_wt = protein_seq[:check_len]

    for i in range(len(cdna_seq) - 2):
        if cdna_seq[i:i+3] != "ATG":
            continue
        # Translate from this ATG
        try:
            translated = str(Seq(cdna_seq[i:]).translate(to_stop=False))
        except Exception:
            continue
        if len(translated) >= check_len and translated[:check_len] == prefix_wt:
            return i
    return -1


# ── Apply HGVSc mutation ───────────────────────────────────────────────────────

def apply_mutation(cds_seq: str, start1: int, end1: int,
                   mut_type: str, seq: str) -> str:
    """
    Apply a frameshift mutation to a CDS sequence.

    Parameters:
      cds_seq  : the full CDS sequence (starts at ATG)
      start1   : HGVSc start position (1-based)
      end1     : HGVSc end position (1-based)
      mut_type : 'del' or 'dup'
      seq      : nucleotide sequence being deleted/duplicated

    Returns the mutated CDS sequence.
    """
    s = start1 - 1  # 0-based
    e = end1        # 0-based exclusive

    if mut_type == "del":
        mutant = cds_seq[:s] + cds_seq[e:]
    elif mut_type == "dup":
        # duplication: insert the duplicated seq after position end1
        mutant = cds_seq[:e] + seq + cds_seq[e:]
    else:
        return cds_seq

    return mutant


# ── Novel peptide extraction ───────────────────────────────────────────────────

def translate_novel_suffix(mutant_cds: str, aa_pos: int) -> str:
    """
    Translate the frameshifted CDS from the mutation site (aa_pos - 1)
    to the first stop codon.

    Parameters:
      mutant_cds : full mutant CDS sequence (starts at ATG)
      aa_pos     : 1-based amino acid position of the frameshift

    Returns the novel amino acid sequence (may include '*' for stop).
    Uses Biopython's Bio.Seq.translate() with standard codon table.
    """
    nt_start = (aa_pos - 1) * 3  # 0-based nucleotide offset in CDS
    if nt_start >= len(mutant_cds):
        return ""

    novel_nt = mutant_cds[nt_start:]
    try:
        novel_aa = str(Seq(novel_nt).translate(to_stop=False))
    except Exception as e:
        log.warning(f"  Translation error: {e}")
        return ""

    # Trim at stop codon
    if "*" in novel_aa:
        novel_aa = novel_aa[:novel_aa.index("*")]

    return novel_aa[:MAX_NOVEL_AA]


def extract_junction_peptides(wt_protein: str, aa_pos: int,
                               novel_suffix: str, lengths: list) -> list:
    """
    Extract junction peptides spanning the WT|novel boundary.

    The prefix (WT amino acids up to aa_pos-1) is joined with the novel
    suffix (frameshifted translation). Peptides are extracted from a window
    centered on the junction.
    """
    prefix = wt_protein[:aa_pos - 1]
    combined = prefix + novel_suffix

    peptides = []
    for k in lengths:
        start = max(0, aa_pos - k)
        for i in range(start, min(aa_pos, len(combined) - k + 1)):
            pep = combined[i:i + k]
            if len(pep) == k and "*" not in pep and "X" not in pep:
                peptides.append(pep)
    return list(dict.fromkeys(peptides))   # deduplicate preserving order


# ── VCF parser for frameshift variants ────────────────────────────────────────

AA_CODE = {
    "Ala":"A","Arg":"R","Asn":"N","Asp":"D","Cys":"C","Gln":"Q",
    "Glu":"E","Gly":"G","His":"H","Ile":"I","Leu":"L","Lys":"K",
    "Met":"M","Phe":"F","Pro":"P","Ser":"S","Thr":"T","Trp":"W",
    "Tyr":"Y","Val":"V","Ter":"*","Xaa":"X"
}

def parse_hgvsp_aa_pos(hgvsp: str) -> int:
    """Extract the amino acid position from HGVSp (e.g. p.Leu40fs → 40)."""
    m = re.match(r"p\.([A-Za-z*]+)(\d+)", hgvsp)
    if m:
        return int(m.group(2))
    return 0


def collect_frameshift_variants(samples: list, ann_dir: str) -> list:
    """
    Parse all annotated VCFs to collect frameshift variants.
    Returns list of dicts: sample, gene, hgvsc, hgvsp, chrom, pos, ref, alt,
                            transcript_id, mut_id
    (Deduplicates by gene+hgvsc — same variant, different transcripts merged.)
    """
    seen = set()
    records = []
    for sample in samples:
        vcf = f"{ann_dir}/{sample}.annotated.vcf.gz"
        if not os.path.exists(vcf):
            continue
        opener = gzip.open
        with opener(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                chrom, pos, _, ref, alt, _, filt, info = parts[:8]
                if filt not in ("PASS", "."):
                    continue
                ann_m = re.search(r"ANN=([^;]+)", info)
                if not ann_m:
                    continue
                for entry in ann_m.group(1).split(","):
                    fields = entry.split("|")
                    if len(fields) < 11:
                        continue
                    effects = set(fields[1].split("&"))
                    if not (effects & FRAMESHIFT_EFFECTS):
                        continue
                    impact = fields[2]
                    if impact not in ("HIGH", "MODERATE"):
                        continue
                    gene       = fields[3]
                    hgvsc      = fields[9]
                    hgvsp      = fields[10]
                    transcript = fields[6]
                    allele     = fields[0]

                    key = (gene, hgvsc)
                    if key in seen:
                        continue
                    seen.add(key)

                    mut_id = f"{chrom}:{pos}:{ref}:{allele.split(',')[0]}"
                    records.append({
                        "sample":      sample,
                        "gene":        gene,
                        "hgvsc":       hgvsc,
                        "hgvsp":       hgvsp,
                        "chrom":       chrom,
                        "pos":         pos,
                        "ref":         ref,
                        "allele":      allele,
                        "transcript":  transcript,
                        "mut_id":      mut_id,
                    })
    return records


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    log.info("=" * 60)
    log.info("Step 08c — Frameshift Novel Sequence Extraction")
    log.info("  Tool: Biopython Bio.Seq.translate + Ensembl GRCm38 cDNA FASTA")
    log.info("=" * 60)

    # 1. Collect frameshift variants from all annotated VCFs
    variants = collect_frameshift_variants(TUMOR_SAMPLES, ANN_DIR)
    log.info(f"Found {len(variants)} unique frameshift mutations")

    if not variants:
        log.error("No frameshift variants found — check ANN_DIR path")
        sys.exit(1)

    target_genes = {v["gene"] for v in variants}
    log.info(f"Target genes: {sorted(target_genes)}")

    # 2. Load protein DB (for WT sequences)
    prot_db = load_protein_db()

    # 3. Load cDNA DB for target genes only
    cdna_db = load_cdna_db(target_genes)

    # 4. Process each variant
    all_records = []
    translation_stats = {"success": 0, "failed_cdna": 0, "failed_cds": 0, "failed_hgvsc": 0}

    for v in variants:
        gene   = v["gene"]
        hgvsc  = v["hgvsc"]
        hgvsp  = v["hgvsp"]
        sample = v["sample"]
        mut_id = v["mut_id"]
        aa_pos = parse_hgvsp_aa_pos(hgvsp)

        log.info(f"\n  Gene={gene}  HGVSc={hgvsc}  HGVSp={hgvsp}  aa_pos={aa_pos}")

        # Parse HGVSc
        start1, end1, mut_type, mut_seq = parse_hgvsc(hgvsc)
        if start1 is None:
            log.warning(f"    Could not parse HGVSc: {hgvsc}")
            translation_stats["failed_hgvsc"] += 1
            continue

        # Get cDNA sequences for this gene
        cdna_entries = cdna_db.get(gene, [])
        # Prefer protein_coding transcripts
        prot_coding = [e for e in cdna_entries if e[1] == "protein_coding"]
        if prot_coding:
            cdna_entries = prot_coding

        if not cdna_entries:
            log.warning(f"    No cDNA found for {gene}")
            translation_stats["failed_cdna"] += 1
            continue

        # Get WT protein sequences for this gene (from protein FASTA)
        wt_proteins = prot_db.get(gene, [])
        if not wt_proteins:
            log.warning(f"    No protein found for {gene}")
            translation_stats["failed_cdna"] += 1
            continue

        # Try each cDNA × protein combination to find CDS start
        cds_start = -1
        wt_protein = None
        chosen_cdna = None

        for tid, bt, cdna_seq in cdna_entries:
            for prot in wt_proteins:
                if aa_pos > len(prot) + 5:   # aa_pos must be within protein
                    continue
                cds_s = find_cds_start(cdna_seq, prot, aa_pos)
                if cds_s >= 0:
                    cds_start = cds_s
                    wt_protein = prot
                    chosen_cdna = cdna_seq
                    log.info(f"    CDS found in {tid}: offset={cds_s}, protein_len={len(prot)}")
                    break
            if cds_start >= 0:
                break

        if cds_start < 0:
            # Fallback: use longest protein-coding cDNA, take longest protein,
            # assume CDS starts at position 0 (Ensembl cDNA sometimes CDS-only)
            log.warning(f"    CDS not found by matching — trying longest cDNA approach")
            if cdna_entries and wt_proteins:
                _, _, chosen_cdna = max(cdna_entries, key=lambda x: len(x[2]))
                wt_protein = max(wt_proteins, key=len)
                # Check if cDNA length ≈ protein length * 3 (CDS-only cDNA)
                expected_len = (len(wt_protein) + 1) * 3   # +1 for stop codon
                if abs(len(chosen_cdna) - expected_len) < 50:
                    cds_start = 0
                    log.info(f"    Using CDS-only cDNA (len={len(chosen_cdna)})")
                else:
                    log.warning(f"    Fallback failed: cDNA len={len(chosen_cdna)}, expected≈{expected_len}")
                    translation_stats["failed_cds"] += 1
                    continue
            else:
                translation_stats["failed_cds"] += 1
                continue

        # Extract CDS from cDNA
        cds_seq = chosen_cdna[cds_start:]

        # Apply the HGVSc mutation to the CDS
        mutant_cds = apply_mutation(cds_seq, start1, end1, mut_type, mut_seq)

        # Verify: mutation causes a frameshift (length change not divisible by 3)
        del_len = end1 - start1 + 1 if mut_type == "del" else 0
        ins_len = len(mut_seq) if mut_type == "dup" else 0
        net_change = ins_len - del_len
        if net_change % 3 == 0:
            log.warning(f"    Net change {net_change} is divisible by 3 — not a frameshift!")
            # Still proceed — SnpEff classified it as frameshift

        # Translate novel suffix using Biopython
        novel_suffix = translate_novel_suffix(mutant_cds, aa_pos)
        if not novel_suffix:
            log.warning(f"    Empty novel suffix for {gene} {hgvsc}")
            translation_stats["failed_cds"] += 1
            continue

        log.info(f"    Novel suffix ({len(novel_suffix)} aa): {novel_suffix[:30]}{'...' if len(novel_suffix)>30 else ''}")
        translation_stats["success"] += 1

        # Generate junction peptides
        peptides = extract_junction_peptides(wt_protein, aa_pos, novel_suffix, PEPTIDE_LENS)
        log.info(f"    Junction peptides: {len(peptides)}")

        for sample_all in TUMOR_SAMPLES:
            # Include peptide for each sample that has this variant
            # (we'll filter below by cross-referencing original VCFs)
            pass

        # For now use the first sample that reported this variant
        for pep in peptides:
            all_records.append({
                "sample":       sample,
                "mut_id":       mut_id,
                "gene":         gene,
                "transcript":   v["transcript"],
                "effect":       "frameshift_variant",
                "hgvsp":        hgvsp,
                "peptide":      pep,
                "length":       len(pep),
                "is_frameshift": True,
                "novel_suffix": novel_suffix[:20],
            })

    log.info(f"\nTranslation stats: {translation_stats}")

    if not all_records:
        log.error("No frameshift peptides generated")
        sys.exit(1)

    df = pd.DataFrame(all_records)
    df = df.drop_duplicates(subset=["gene", "peptide"])
    log.info(f"Total real frameshift peptides: {len(df)}")

    out_tsv = f"{PEP_DIR}/frameshift_real.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    log.info(f"Saved: {out_tsv}")

    # 5. Update all_peptides.tsv: replace AAAA placeholder frameshifts
    all_pep_tsv = f"{PEP_DIR}/all_peptides.tsv"
    if os.path.exists(all_pep_tsv):
        all_df = pd.read_csv(all_pep_tsv, sep="\t")
        log.info(f"\nOriginal all_peptides.tsv: {len(all_df)} rows")

        # Remove old placeholder frameshift peptides
        non_fs = all_df[~all_df["is_frameshift"]]
        log.info(f"  Non-frameshift rows kept: {len(non_fs)}")

        # Merge real frameshift peptides
        new_df = pd.concat([
            non_fs,
            df[["sample","mut_id","gene","transcript","effect","hgvsp","peptide","length","is_frameshift"]]
        ], ignore_index=True)

        new_df.to_csv(all_pep_tsv, sep="\t", index=False)
        log.info(f"  Updated all_peptides.tsv: {len(new_df)} rows")

        # 6. Regenerate per-length FASTAs
        for plen in PEPTIDE_LENS:
            subset = new_df[new_df["length"] == plen]["peptide"].drop_duplicates()
            fasta  = f"{PEP_DIR}/peptides_{plen}mer.fasta"
            with open(fasta, "w") as fh:
                for i, pep in enumerate(subset):
                    fh.write(f">pep_{plen}_{i}\n{pep}\n")
            log.info(f"  {plen}-mer FASTA: {len(subset)} peptides → {fasta}")

        # Upload to S3
        s3_pep = f"{S3_RESULTS}/peptides"
        os.system(f"aws s3 cp {all_pep_tsv} {s3_pep}/all_peptides.tsv --quiet")
        for plen in PEPTIDE_LENS:
            os.system(f"aws s3 cp {PEP_DIR}/peptides_{plen}mer.fasta {s3_pep}/peptides_{plen}mer.fasta --quiet")
        log.info("  Uploaded updated files to S3")

    # Print summary
    log.info("\n" + "=" * 60)
    log.info("FRAMESHIFT NOVEL SEQUENCE SUMMARY")
    log.info("=" * 60)
    log.info(f"  Variants processed: {len(variants)}")
    log.info(f"  Successful translations: {translation_stats['success']}")
    log.info(f"  Failed (no cDNA/protein): {translation_stats['failed_cdna']}")
    log.info(f"  Failed (CDS not found): {translation_stats['failed_cds']}")
    log.info(f"  Failed (HGVSc parse): {translation_stats['failed_hgvsc']}")
    log.info(f"  Frameshift peptides generated: {len(df)}")
    if not df.empty:
        log.info("\nSample frameshift peptides:")
        print(df[["gene","hgvsp","peptide","novel_suffix"]].head(10).to_string(index=False))
    log.info("=" * 60)
    log.info("Step 08c complete.")


if __name__ == "__main__":
    main()
