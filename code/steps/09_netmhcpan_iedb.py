#!/usr/bin/env python3
"""
Step 9 — MHC-I Binding Prediction via IEDB REST API

Predicts H-2-Kb and H-2-Db binding for 8/9/10-mer peptides using
netMHCpan_el via the IEDB public API (no registration required).

Batches peptides (100 per request) to respect API limits.
Retries on failure with exponential back-off.

Output:
  netmhcpan/all_predictions.tsv    — all SB/WB hits merged
  netmhcpan/<len>mer_<allele>.tsv  — per-length per-allele
  S3: <S3_RESULTS>/netmhcpan/

Usage:
  python3 09_netmhcpan_iedb.py
"""

import os, time, logging, json, sys
import requests
import pandas as pd

logging.basicConfig(level=logging.INFO,
                    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger(__name__)

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR    = "/home/ec2-user"
RESULTS_DIR = os.environ.get("RESULTS_DIR", f"{BASE_DIR}/results/res_20260311_225555")
S3_RESULTS  = os.environ.get("S3_RESULTS",  "s3://neoantigen2026-rerun/results/res_20260311_225555")
S3_FINAL    = "s3://neoantigen2026-rerun/final-results"

PEP_DIR  = f"{RESULTS_DIR}/peptides"
OUT_DIR  = f"{RESULTS_DIR}/netmhcpan"
S3_OUT   = f"{S3_RESULTS}/netmhcpan"

ALLELES      = ["H-2-Kb", "H-2-Db"]
PEPTIDE_LENS = [8, 9, 10]
METHOD       = "netmhcpan_el"
IEDB_URL     = "http://tools-cluster-interface.iedb.org/tools_api/mhci/"

SB_THRESHOLD = 0.5   # % rank — strong binder
WB_THRESHOLD = 2.0   # % rank — weak binder

BATCH_SIZE   = 100   # peptides per API request
MAX_RETRIES  = 5
RETRY_DELAY  = 10    # seconds base delay (exponential back-off)


# ── Helpers ───────────────────────────────────────────────────────────────────

def read_fasta(fasta_path: str) -> list[str]:
    peptides = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(">"):
                peptides.append(line)
    return peptides


def query_iedb(peptides: list[str], allele: str, length: int) -> list[dict]:
    """Submit a batch of peptides to IEDB API and return parsed rows."""
    seq_text = "\n".join(peptides)
    payload = {
        "method":        METHOD,
        "sequence_text": seq_text,
        "allele":        allele,
        "length":        str(length),
    }
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.post(IEDB_URL, data=payload, timeout=120)
            if resp.status_code == 200:
                return parse_iedb_response(resp.text, allele, length)
            else:
                log.warning(f"  HTTP {resp.status_code} on attempt {attempt}: {resp.text[:200]}")
        except requests.RequestException as e:
            log.warning(f"  Request error attempt {attempt}: {e}")
        wait = RETRY_DELAY * (2 ** (attempt - 1))
        log.info(f"  Retrying in {wait}s...")
        time.sleep(wait)
    log.error(f"  All {MAX_RETRIES} attempts failed for {allele} {length}-mer batch")
    return []


def parse_iedb_response(text: str, allele: str, length: int) -> list[dict]:
    """
    IEDB returns tab-separated text. Typical columns:
      allele  seq_num  start  end  length  peptide  core  icore
      score   rank     ... (affinity columns vary by method)

    We extract: peptide, el_rank (% rank EL), classify as SB/WB/NB.
    """
    rows = []
    # Response columns: allele seq_num start end length peptide core icore score percentile_rank
    lines = [l for l in text.strip().split("\n") if l and not l.startswith("allele\t")]
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 10:
            continue
        try:
            peptide = parts[5].strip()
            rank_el = float(parts[9])   # percentile_rank column
            binder = ("SB" if rank_el <= SB_THRESHOLD
                      else "WB" if rank_el <= WB_THRESHOLD
                      else "NB")
            rows.append({
                "allele":       allele,
                "plen":         length,
                "peptide":      peptide,
                "el_rank":      rank_el,
                "binder_class": binder,
            })
        except (IndexError, ValueError):
            continue
    return rows


def run_predictions(peptides: list[str], allele: str, length: int) -> pd.DataFrame:
    """Run all batches for one allele/length combination."""
    all_rows = []
    n_batches = (len(peptides) + BATCH_SIZE - 1) // BATCH_SIZE
    log.info(f"  {allele} {length}-mer: {len(peptides)} peptides in {n_batches} batches")

    for i in range(0, len(peptides), BATCH_SIZE):
        batch  = peptides[i:i + BATCH_SIZE]
        batch_n = i // BATCH_SIZE + 1
        rows   = query_iedb(batch, allele, length)
        all_rows.extend(rows)
        log.info(f"    Batch {batch_n}/{n_batches}: {len(rows)} results")
        time.sleep(1)   # be polite to the API

    return pd.DataFrame(all_rows)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    log.info("=" * 60)
    log.info("Step 9 — MHC-I Binding Prediction (IEDB API / netMHCpan_el)")
    log.info(f"  Alleles : {ALLELES}")
    log.info(f"  Lengths : {PEPTIDE_LENS}")
    log.info(f"  Method  : {METHOD}")
    log.info(f"  SB/WB   : <{SB_THRESHOLD}% / <{WB_THRESHOLD}% rank_EL")
    log.info("=" * 60)

    all_results = []

    for length in PEPTIDE_LENS:
        fasta = f"{PEP_DIR}/peptides_{length}mer.fasta"
        if not os.path.exists(fasta):
            os.system(f"aws s3 cp {S3_RESULTS}/peptides/peptides_{length}mer.fasta {fasta} --quiet")
        if not os.path.exists(fasta):
            log.warning(f"No {length}-mer FASTA found — skipping")
            continue

        peptides = read_fasta(fasta)
        log.info(f"\n{'='*40}")
        log.info(f"{length}-mer peptides: {len(peptides)}")
        log.info(f"{'='*40}")

        for allele in ALLELES:
            allele_safe = allele.replace("-", "_")
            out_tsv = f"{OUT_DIR}/{length}mer_{allele_safe}.tsv"

            # Skip if already done
            if os.path.exists(out_tsv) and os.path.getsize(out_tsv) > 50:
                log.info(f"[SKIP] {length}mer {allele} — already exists")
                df = pd.read_csv(out_tsv, sep="\t")
                all_results.append(df)
                continue

            df = run_predictions(peptides, allele, length)

            if df.empty:
                log.warning(f"  No results for {allele} {length}-mer")
                continue

            sb = (df["binder_class"] == "SB").sum()
            wb = (df["binder_class"] == "WB").sum()
            log.info(f"  {allele} {length}-mer → SB={sb}  WB={wb}  total={len(df)}")

            df.to_csv(out_tsv, sep="\t", index=False)
            os.system(f"aws s3 cp {out_tsv} {S3_OUT}/{length}mer_{allele_safe}.tsv --quiet")
            all_results.append(df)

    if not all_results:
        log.error("No predictions obtained. Check IEDB API connectivity.")
        sys.exit(1)

    # Merge all results
    merged = pd.concat(all_results, ignore_index=True)

    # Keep only binders (SB + WB) — full table
    binders = merged[merged["binder_class"] != "NB"].copy()
    binders_tsv = f"{OUT_DIR}/all_predictions.tsv"
    binders.to_csv(binders_tsv, sep="\t", index=False)
    os.system(f"aws s3 cp {binders_tsv} {S3_OUT}/all_predictions.tsv --quiet")
    os.system(f"aws s3 cp {binders_tsv} {S3_FINAL}/netmhcpan/all_predictions.tsv --quiet")

    total_sb = (binders["binder_class"] == "SB").sum()
    total_wb = (binders["binder_class"] == "WB").sum()

    log.info("\n" + "=" * 60)
    log.info("BINDING PREDICTION SUMMARY")
    log.info("=" * 60)
    log.info(f"  Strong binders (Rank_EL < 0.5%): {total_sb}")
    log.info(f"  Weak binders   (Rank_EL < 2.0%): {total_wb}")
    log.info(f"  Total peptides with binding:      {len(binders)}")

    # Summary per allele
    for allele in ALLELES:
        a_df = binders[binders["allele"] == allele]
        sb = (a_df["binder_class"] == "SB").sum()
        wb = (a_df["binder_class"] == "WB").sum()
        log.info(f"  {allele}: SB={sb}  WB={wb}")

    log.info(f"  Saved: {binders_tsv}")
    log.info("=" * 60)
    log.info("Step 9 complete.")

    # Print top 20 strong binders
    if total_sb > 0:
        top = (binders[binders["binder_class"] == "SB"]
               .sort_values("el_rank")
               .head(20))
        print("\nTop 20 Strong Binders (lowest Rank_EL%):")
        print(top[["allele", "plen", "peptide", "el_rank", "binder_class"]].to_string(index=False))


if __name__ == "__main__":
    main()
