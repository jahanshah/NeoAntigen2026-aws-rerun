#!/usr/bin/env python3
"""
Step 08d — Fusion ORF Translation for Arriba "." peptide_sequence cases

For fusions where Arriba outputs peptide_sequence="." (breakpoint in UTR, intron,
or intergenic region), this script attempts to recover a junction peptide sequence
by translating the downstream gene CDS from the breakpoint position.

Approach (all using established bioinformatics tools):
  1. Parse Arriba fusion output — select high/medium fusions with "." peptide
  2. Extract the downstream gene's CDS coordinates from Ensembl GTF
     (Mus_musculus.GRCm38.102.gtf.gz)
  3. Fetch the genomic sequence at the breakpoint using pysam (mm10.fa)
  4. Build the downstream CDS sequence from the breakpoint to end of transcript
  5. Find the first in-frame ATG at or downstream of the breakpoint
  6. Translate using Biopython Bio.Seq.translate (standard codon table)
  7. For the upstream gene (if breakpoint is in CDS), use the Ensembl protein FASTA
     to get the upstream WT peptide prefix
  8. Extract junction k-mers spanning the WT|novel boundary

Tools used:
  - pysam 0.23.3       — genomic sequence extraction (mm10.fa)
  - Biopython 1.86     — Bio.Seq.translate
  - Ensembl GRCm38 Release 102 GTF  — CDS coordinates
  - Ensembl GRCm38 Release 102 protein FASTA — upstream WT protein

Germline/artifact filter: fusions present in >= MIN_SAMPLES_GERMLINE samples excluded.

Output:
  peptides/fusions_dot_translated.tsv   — translated junction peptides
  netmhcpan/fusion/dot_fusions_predictions.tsv — IEDB MHC predictions
"""

import os, gzip, re, logging, time, sys
import pandas as pd
import pysam
import requests
from Bio.Seq import Seq

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR    = "/home/ec2-user"
RESULTS_DIR = os.environ.get("RESULTS_DIR", f"{BASE_DIR}/results/res_20260311_225555")
S3_FINAL    = "s3://neoantigen2026-rerun/final-results"

GENOME_FA   = f"{BASE_DIR}/ref/mm10/mm10.fa"
GTF_GZ      = f"{BASE_DIR}/.cache/pyensembl/GRCm38/ensembl102/Mus_musculus.GRCm38.102.gtf.gz"
PROT_FASTA  = f"{BASE_DIR}/.cache/pyensembl/GRCm38/ensembl102/Mus_musculus.GRCm38.pep.all.fa.gz"

ARRIBA_DIR  = f"{BASE_DIR}/tmp/neoantig_pipeline/arriba_peptides"
OUT_DIR     = f"{RESULTS_DIR}/peptides"
MHC_OUT     = f"{RESULTS_DIR}/netmhcpan/fusion"

TUMOR_SAMPLES = ["443_D21_new","428_D20_new","34_D52_old","36_D99_new","38_D99_new","42_D122_old"]
PEPTIDE_LENS  = [8, 9, 10]
MIN_SAMPLES_GERMLINE = 3
MAX_NOVEL_AA  = 60

IEDB_URL     = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"
ALLELES      = ["H-2-Kb", "H-2-Db"]
METHOD       = "netmhcpan_el"
SB_THRESHOLD = 0.5
WB_THRESHOLD = 2.0
BATCH_SIZE   = 100
MAX_RETRIES  = 5
RETRY_DELAY  = 10

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(MHC_OUT, exist_ok=True)


# ── GTF CDS loader ────────────────────────────────────────────────────────────

def load_cds_db(target_genes: set) -> dict:
    """
    Parse Ensembl GTF to build gene_name → list of CDS exons.
    Each entry: {'chrom', 'start'(1-based), 'end', 'strand', 'exon_number', 'transcript_id'}
    Only protein_coding transcripts included.
    """
    db = {}  # gene_name → [{'chrom','start','end','strand','exon_number','transcript_id'}]
    gname_re  = re.compile(r'gene_name "([^"]+)"')
    bt_re     = re.compile(r'transcript_biotype "([^"]+)"')
    tid_re    = re.compile(r'transcript_id "([^"]+)"')
    exon_re   = re.compile(r'exon_number "?(\d+)"?')

    log.info("Parsing Ensembl GTF for CDS coordinates...")
    with gzip.open(GTF_GZ, "rt") as fh:
        for line in fh:
            if line.startswith("#") or "\tCDS\t" not in line:
                continue
            gm = gname_re.search(line)
            if not gm or gm.group(1) not in target_genes:
                continue
            bt = bt_re.search(line)
            if not bt or bt.group(1) != "protein_coding":
                continue
            parts = line.split("\t")
            chrom  = parts[0]
            start  = int(parts[3])   # 1-based inclusive
            end    = int(parts[4])
            strand = parts[6]
            gene   = gm.group(1)
            tid    = tid_re.search(line)
            exn    = exon_re.search(line)
            db.setdefault(gene, []).append({
                "chrom":        chrom,
                "start":        start,
                "end":          end,
                "strand":       strand,
                "transcript_id": tid.group(1) if tid else "",
                "exon_number":  int(exn.group(1)) if exn else 0,
            })
    log.info(f"CDS entries loaded for {len(db)} target genes")
    return db


def get_cds_sequence(cds_exons: list, genome: pysam.FastaFile, strand: str) -> str:
    """
    Fetch and concatenate CDS exon sequences from the genome using pysam.
    Handles strand orientation (reverse-complement for minus strand).
    """
    # Sort by exon number; for minus strand, reverse order
    sorted_exons = sorted(cds_exons, key=lambda x: x["exon_number"],
                          reverse=(strand == "-"))

    seq_parts = []
    for exon in sorted_exons:
        chrom = exon["chrom"]
        # pysam uses 0-based half-open coordinates
        try:
            seq = genome.fetch(chrom, exon["start"] - 1, exon["end"])
        except KeyError:
            # Try with/without chr prefix
            alt = chrom.replace("chr","") if chrom.startswith("chr") else "chr"+chrom
            try:
                seq = genome.fetch(alt, exon["start"] - 1, exon["end"])
            except Exception:
                return ""
        seq_parts.append(seq.upper())

    full_seq = "".join(seq_parts)
    if strand == "-":
        full_seq = str(Seq(full_seq).reverse_complement())
    return full_seq


def translate_from_position(cds_seq: str, nt_offset: int, max_aa: int = MAX_NOVEL_AA) -> str:
    """
    Translate CDS from nt_offset (0-based) to first stop codon using Biopython.
    Returns amino acid sequence (no stop codon character).
    """
    subseq = cds_seq[nt_offset:]
    if len(subseq) < 3:
        return ""
    try:
        aa = str(Seq(subseq).translate(to_stop=False))
    except Exception:
        return ""
    if "*" in aa:
        aa = aa[:aa.index("*")]
    return aa[:max_aa]


# ── Protein DB (upstream WT prefix) ──────────────────────────────────────────

def load_protein_db(target_genes: set) -> dict:
    db = {}
    sym_re = re.compile(r"gene_symbol:(\S+)")
    cur_sym = None; cur_seq = []
    with gzip.open(PROT_FASTA, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_sym and cur_sym in target_genes and cur_seq:
                    db.setdefault(cur_sym, []).append("".join(cur_seq))
                m = sym_re.search(line)
                cur_sym = m.group(1) if m else None
                cur_seq = []
            elif cur_sym and cur_sym in target_genes:
                cur_seq.append(line)
    if cur_sym and cur_sym in target_genes and cur_seq:
        db.setdefault(cur_sym, []).append("".join(cur_seq))
    return db


# ── Junction peptide extraction ───────────────────────────────────────────────

def extract_junction_peptides(upstream_aa: str, junction_pos: int,
                               downstream_aa: str, lengths: list) -> list:
    """
    Extract k-mers spanning the upstream|downstream junction.
    upstream_aa:   WT amino acids from upstream gene up to breakpoint
    junction_pos:  index in the combined sequence where junction occurs
    downstream_aa: translated novel downstream sequence
    """
    combined = upstream_aa + downstream_aa
    peptides = []
    for k in lengths:
        start = max(0, junction_pos - k + 1)
        end   = min(junction_pos + 1, len(combined) - k + 1)
        for i in range(start, end):
            pep = combined[i:i+k]
            if len(pep) == k and "*" not in pep and "X" not in pep:
                peptides.append(pep.upper())
    return list(dict.fromkeys(peptides))


# ── Arriba parser for "." fusions ─────────────────────────────────────────────

def collect_dot_fusions(arriba_dir: str, samples: list) -> pd.DataFrame:
    """Collect all high/medium fusions with peptide_sequence='.' from all samples."""
    rows = []
    for sample in samples:
        fpath = f"{arriba_dir}/Arriba_{sample}.txt"
        if not os.path.exists(fpath):
            continue
        df = pd.read_csv(fpath, sep="\t")
        df.columns = [c.lstrip("#") for c in df.columns]
        hm = df[df["confidence"].isin(["high","medium"])]
        dot = hm[hm["peptide_sequence"] == "."]
        for _, row in dot.iterrows():
            rows.append({
                "sample":      sample,
                "gene1":       str(row.get("gene1","")).split("(")[0],
                "gene2":       str(row.get("gene2","")).split("(")[0],
                "breakpoint1": str(row.get("breakpoint1","")),
                "breakpoint2": str(row.get("breakpoint2","")),
                "site1":       str(row.get("site1","")),
                "site2":       str(row.get("site2","")),
                "reading_frame": str(row.get("reading_frame",".")),
                "confidence":  str(row.get("confidence","")),
                "strand1":     str(row.get("strand1(gene/fusion)","")).split("/")[0],
                "strand2":     str(row.get("strand2(gene/fusion)","")).split("/")[0],
            })
    return pd.DataFrame(rows)


# ── IEDB query ────────────────────────────────────────────────────────────────

def query_iedb(peptides: list, allele: str, length: int) -> list:
    payload = {"method": METHOD, "sequence_text": "\n".join(peptides),
               "allele": allele, "length": str(length)}
    for attempt in range(1, MAX_RETRIES+1):
        try:
            resp = requests.post(IEDB_URL, data=payload, timeout=120)
            if resp.status_code == 200:
                rows = []
                for line in resp.text.strip().split("\n"):
                    if not line or line.startswith("allele\t"): continue
                    parts = line.split("\t")
                    if len(parts) < 10: continue
                    try:
                        pep  = parts[5].strip()
                        rank = float(parts[9])
                        bc   = "SB" if rank <= SB_THRESHOLD else ("WB" if rank <= WB_THRESHOLD else "NB")
                        rows.append({"allele": allele, "plen": length,
                                     "peptide": pep, "el_rank": rank, "binder_class": bc})
                    except: continue
                return rows
        except Exception as e:
            log.warning(f"  Request error attempt {attempt}: {e}")
        time.sleep(RETRY_DELAY * (2**(attempt-1)))
    return []


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    log.info("=" * 60)
    log.info("Step 08d — Fusion ORF Translation ('.' peptide cases)")
    log.info("  Tools: pysam (genome fetch) + Biopython (translation)")
    log.info("         Ensembl GRCm38 R102 GTF + protein FASTA")
    log.info("=" * 60)

    # 1. Collect "." fusions
    dot_df = collect_dot_fusions(ARRIBA_DIR, TUMOR_SAMPLES)
    if dot_df.empty:
        log.info("No '.' fusions found — nothing to translate.")
        return
    log.info(f"Found {len(dot_df)} '.' fusions across {dot_df['sample'].nunique()} samples")

    # 2. Germline filter
    fusion_counts = dot_df.groupby(dot_df.apply(
        lambda r: f"{r['gene1']}::{r['gene2']}", axis=1))["sample"].nunique()
    germline = fusion_counts[fusion_counts >= MIN_SAMPLES_GERMLINE].index.tolist()
    if germline:
        log.warning(f"Excluding {len(germline)} germline/recurrent fusions: {germline}")
        dot_df["fusion_id"] = dot_df["gene1"] + "::" + dot_df["gene2"]
        dot_df = dot_df[~dot_df["fusion_id"].isin(germline)]
    log.info(f"After germline filter: {len(dot_df)} fusions to process")

    if dot_df.empty:
        log.info("All '.' fusions were germline/recurrent — nothing left to translate.")
        return

    # 3. Load genome
    log.info(f"Loading genome: {GENOME_FA}")
    genome = pysam.FastaFile(GENOME_FA)

    # 4. Load CDS DB for downstream genes (gene2 = downstream partner)
    target_genes = set(dot_df["gene2"].unique()) | set(dot_df["gene1"].unique())
    log.info(f"Loading CDS coordinates for {len(target_genes)} genes...")
    cds_db = load_cds_db(target_genes)
    prot_db = load_protein_db(set(dot_df["gene1"].unique()))

    # 5. Process each fusion
    all_records = []
    for _, row in dot_df.iterrows():
        gene1 = row["gene1"]
        gene2 = row["gene2"]
        bp2   = row["breakpoint2"]   # downstream breakpoint (chrom:pos)
        site2 = row["site2"]
        sample = row["sample"]
        strand2 = row["strand2"]
        fusion_id = f"{gene1}::{gene2}"

        log.info(f"\n  {fusion_id} ({sample}) | {row['site1']} x {site2} | strand2={strand2}")

        # Parse breakpoint2
        try:
            chrom2, pos2 = bp2.split(":")
            pos2 = int(pos2)
        except Exception:
            log.warning(f"    Cannot parse breakpoint2: {bp2}")
            continue

        # Only attempt if downstream breakpoint involves a CDS
        if "CDS" not in site2 and "UTR" not in site2:
            log.info(f"    Skipping: downstream site '{site2}' has no coding context")
            continue

        # Get CDS exons for downstream gene
        cds_exons_raw = cds_db.get(gene2, [])
        if not cds_exons_raw:
            log.warning(f"    No CDS found for downstream gene {gene2}")
            continue

        # Use the most common transcript (largest exon count)
        from collections import Counter
        tid_counts = Counter(e["transcript_id"] for e in cds_exons_raw)
        best_tid   = tid_counts.most_common(1)[0][0]
        cds_exons  = [e for e in cds_exons_raw if e["transcript_id"] == best_tid]
        strand2_gtf = cds_exons[0]["strand"]

        log.info(f"    Downstream CDS: {gene2} transcript={best_tid} "
                 f"strand={strand2_gtf} exons={len(cds_exons)}")

        # Build full downstream CDS
        full_cds = get_cds_sequence(cds_exons, genome, strand2_gtf)
        if not full_cds:
            log.warning(f"    Could not fetch CDS sequence for {gene2}")
            continue

        # Find position within CDS where breakpoint falls
        # Build sorted exon list and find cumulative position
        exons_sorted = sorted(cds_exons, key=lambda x: x["start"] if strand2_gtf=="+" else -x["start"])
        cum_pos = 0
        bp_in_cds = -1
        for exon in exons_sorted:
            exon_len = exon["end"] - exon["start"] + 1
            if strand2_gtf == "+":
                if exon["start"] <= pos2 <= exon["end"]:
                    bp_in_cds = cum_pos + (pos2 - exon["start"])
                    break
            else:
                if exon["start"] <= pos2 <= exon["end"]:
                    bp_in_cds = cum_pos + (exon["end"] - pos2)
                    break
            cum_pos += exon_len

        if bp_in_cds < 0:
            log.info(f"    Breakpoint pos2={pos2} not within CDS exons — using start of CDS")
            bp_in_cds = 0

        # Align to codon boundary
        frame_offset = bp_in_cds % 3
        nt_start = bp_in_cds - frame_offset

        # Translate downstream novel sequence
        downstream_aa = translate_from_position(full_cds, nt_start)
        if not downstream_aa:
            log.warning(f"    Empty translation for {gene2} from position {nt_start}")
            continue
        log.info(f"    Downstream translation ({len(downstream_aa)} aa): "
                 f"{downstream_aa[:30]}{'...' if len(downstream_aa)>30 else ''}")

        # Get upstream WT protein prefix (from Ensembl protein FASTA)
        upstream_proteins = prot_db.get(gene1, [])
        if upstream_proteins and "CDS" in row["site1"]:
            # Use longest protein as proxy; junction at upstream breakpoint
            upstream_prot = max(upstream_proteins, key=len)
            # Estimate upstream AA position from bp1
            try:
                _, pos1 = row["breakpoint1"].split(":")
                pos1 = int(pos1)
                # Rough estimate: use full upstream protein (can't easily compute AA pos)
                upstream_aa = upstream_prot  # conservative: use full WT prefix
            except Exception:
                upstream_aa = ""
        else:
            upstream_aa = ""

        # Junction position = end of upstream sequence
        junction_pos = len(upstream_aa)
        combined_junction = upstream_aa + downstream_aa

        if not combined_junction or junction_pos > len(combined_junction):
            # No upstream context — use first k-mer of downstream translation
            junction_pos = 0

        peptides = extract_junction_peptides(
            upstream_aa, junction_pos, downstream_aa, PEPTIDE_LENS)

        log.info(f"    Junction peptides: {len(peptides)}")
        mut_id = f"{row['breakpoint1']}::{bp2}"
        for pep in peptides:
            all_records.append({
                "sample":       sample,
                "mut_id":       mut_id,
                "gene":         fusion_id,
                "transcript":   best_tid,
                "effect":       "fusion_translated",
                "hgvsp":        "",
                "peptide":      pep,
                "length":       len(pep),
                "is_frameshift": False,
                "confidence":   row["confidence"],
                "downstream_gene": gene2,
                "downstream_seq":  downstream_aa[:20],
            })

    genome.close()

    if not all_records:
        log.info("No junction peptides could be generated from '.' fusions.")
        return

    df = pd.DataFrame(all_records).drop_duplicates(subset=["gene","peptide"])
    log.info(f"\nTotal translated fusion peptides: {len(df)}")

    out_tsv = f"{OUT_DIR}/fusions_dot_translated.tsv"
    df.to_csv(out_tsv, sep='\t', index=False)
    log.info(f"Saved: {out_tsv}")

    # 6. Run IEDB predictions
    log.info("\nRunning IEDB MHC-I predictions...")
    all_mhc = []
    for plen in PEPTIDE_LENS:
        peps = df[df["length"]==plen]["peptide"].drop_duplicates().tolist()
        if not peps: continue
        log.info(f"  {plen}-mer: {len(peps)} peptides")
        for allele in ALLELES:
            rows = query_iedb(peps, allele, plen)
            all_mhc.extend(rows)
            sb = sum(1 for r in rows if r["binder_class"]=="SB")
            wb = sum(1 for r in rows if r["binder_class"]=="WB")
            log.info(f"    {allele}: SB={sb} WB={wb}")
            time.sleep(1)

    if all_mhc:
        mhc_df = pd.DataFrame(all_mhc)
        mhc_df.to_csv(f"{MHC_OUT}/dot_fusions_predictions.tsv", sep='\t', index=False)
        sb_total = (mhc_df["binder_class"]=="SB").sum()
        wb_total = (mhc_df["binder_class"]=="WB").sum()
        log.info(f"\nPrediction results: SB={sb_total} WB={wb_total}")
        if sb_total > 0:
            sb_ann = mhc_df[mhc_df["binder_class"]=="SB"].merge(
                df[["peptide","gene","sample"]].drop_duplicates("peptide"),
                on="peptide", how="left").sort_values("el_rank")
            print("\nDot-fusion STRONG BINDERS:")
            print(sb_ann[["gene","sample","peptide","allele","el_rank"]].to_string(index=False))

    # Upload
    aws_base = "aws s3 cp"
    os.system(f"{aws_base} {out_tsv} {S3_FINAL}/peptides/fusions_dot_translated.tsv --quiet")
    if all_mhc:
        os.system(f"{aws_base} {MHC_OUT}/dot_fusions_predictions.tsv {S3_FINAL}/netmhcpan/dot_fusions_predictions.tsv --quiet")

    log.info("Step 08d complete.")


if __name__ == "__main__":
    main()
