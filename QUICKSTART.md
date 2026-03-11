# NeoAntigen 2026 — WES Pipeline Quickstart

Mouse (mm10) whole-exome somatic neoantigen calling pipeline.
Tumor/normal longitudinal samples. Runs on AWS EC2 with S3 storage.

---

## Requirements

- AWS EC2 instance (recommended: `r5.2xlarge` or larger, 64 GB RAM)
- AWS CLI configured with access to `s3://neoantigen2026-rerun`
- [Miniforge3](https://github.com/conda-forge/miniforge) installed at `~/miniforge3/`
- Git, Bash 4+, R 4+

---

## 1. Clone the Repository

```bash
git clone https://github.com/jahanshah/NeoAntigen2026-aws-rerun.git
cd NeoAntigen2026-aws-rerun
```

---

## 2. Configure Your Environment

Copy the environment template and edit for your setup:

```bash
cp env.sh ~/pipeline_env.sh
nano ~/pipeline_env.sh   # adjust BASE_DIR, CONDA_BIN, THREADS if needed
source ~/pipeline_env.sh
```

Then update `code/config.sh` if your paths differ from defaults (see comments in that file).

---

## 3. Install Dependencies & Download References

Run once per new instance:

```bash
bash code/setup_pipeline.sh
```

This installs: `samtools`, `bcftools`, `gatk4`, `picard`, `snpEff`, `netMHCpan`, `PyClone`, and Python packages.
It also downloads the mm10 reference FASTA, exon BED, and MGP SNPs for BQSR — **or pulls them from S3** if already uploaded.

> **Note:** `netMHCpan` requires manual registration and installation from
> https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/

---

## 4. S3 Data Layout

```
s3://neoantigen2026-rerun/
├── data/
│   ├── bam/wes/            ← raw WES BAM + BAI files (input)
│   └── reference/wes/      ← mm10.fa, .fai, .dict, exons.bed, mgp_snps.vcf.gz
└── results/
    └── res_YYYYMMDD_HHMMSS/ ← one folder per pipeline run (auto-generated)
```

---

## 5. Samples

| Sample | Role | Timepoint |
|---|---|---|
| 423_D0_old | Normal (primary) | D0 |
| 424_D0__old | Normal (PoN) | D0 |
| 428_D20_new | Tumor | D20 |
| 34_D52_old | Tumor | D52 |
| 36_D99_new | Tumor | D99 |
| 38_D99_new | Tumor | D99 |
| 42_D122_old | Tumor | D122 |

---

## 6. Run the Full Pipeline

```bash
bash code/run_pipeline.sh
```

Each run automatically creates a timestamped results folder:
- **Local:** `~/results/res_YYYYMMDD_HHMMSS/`
- **S3:** `s3://neoantigen2026-rerun/results/res_YYYYMMDD_HHMMSS/`

### Run a subset of steps

```bash
bash code/run_pipeline.sh 5      # run steps 5 through 12
bash code/run_pipeline.sh 3 7    # run steps 3 to 7 only
```

### Run a single step

```bash
bash code/steps/01_bam_preprocess.sh
```

> When running a single step directly, a new `RUN_ID` is auto-generated.
> To continue an existing run, export the same ID first:
> ```bash
> export RUN_ID=res_20260311_184352
> bash code/steps/02_mutect2.sh
> ```

---

## 7. Pipeline Steps

| Step | Script | Description |
|---|---|---|
| 1 | `01_bam_preprocess.sh` | Sort → MarkDuplicates → BQSR → ValidateSamFile |
| 2 | `02_mutect2.sh` | Somatic variant calling (Mutect2 + Panel of Normals) |
| 3 | `03_filter_vcf.sh` | FilterMutectCalls + SelectVariants (exons only) |
| 4 | `04_annotate.sh` | SnpEff functional annotation (GRCm38.86) |
| 5 | `05_pyclone_prep.py` | Convert annotated VCFs → PyClone YAML input |
| 6 | `06_pyclone_run.sh` | PyClone MCMC clonal inference (beta-binomial) |
| 7 | `07_pyclone_tables.R` | PyClone cluster/loci tables + figures |
| 8 | `08_peptide_extract.py` | Extract 8/9/10-mer missense/frameshift peptides |
| 9 | `08b_fusion_peptides.py` | Extract fusion junction peptides (Arriba) |
| 10 | `09_netmhcpan.sh` | netMHCpan H-2Kb/H-2Db binding prediction |
| 11 | `11_rnaseq_expression.R` | RNASeq batch correction + expression matrix |
| 12 | `10_filter_results.R` | Rank candidates (binding + clonality + expression) |

### Skip / Reprocess Behaviour (Step 1)

When Step 1 runs, it checks S3 for previously preprocessed BAMs for the current `RUN_ID`.
For any sample already found, it will interactively ask:

```
Sample '423_D0_old' was previously preprocessed on S3 (2026-03-11 18:43).
Skip or reprocess? [s=skip / r=reprocess]:
```

---

## 8. Monitor a Running Job

```bash
# Follow the preprocessing log
tail -f ~/results/bam_preprocessing.log

# Follow the latest run log (if piped)
tail -f ~/results/res_<RUN_ID>/pipeline.log

# Check what's been uploaded to S3 for this run
aws s3 ls s3://neoantigen2026-rerun/results/res_<RUN_ID>/ --recursive
```

---

## 9. Find Results

### List all runs
```bash
aws s3 ls s3://neoantigen2026-rerun/results/
```

### Browse a specific run
```bash
aws s3 ls s3://neoantigen2026-rerun/results/res_20260311_184352/ --recursive
```

### Download results locally
```bash
RUN_ID=res_20260311_184352
aws s3 sync s3://neoantigen2026-rerun/results/${RUN_ID}/ ~/results/${RUN_ID}/
```

### Key output files

| File | Description |
|---|---|
| `preprocessed_bams/*.preproc.bam` | BQSR-corrected BAMs |
| `vcf/*.filtered.vcf.gz` | Filtered somatic variants (per tumor) |
| `vcf/*.annotated.vcf` | SnpEff-annotated variants |
| `pyclone/output/` | Clonal inference tables + figures |
| `peptides/*.fasta` | Candidate neoantigen peptides |
| `netmhcpan/*.xls` | MHC binding predictions |
| `final_candidates.tsv` | Ranked neoantigen candidates |

---

## 10. Evaluate BAMs Before Running

To check all BAMs are suitable for Mutect2 and VarScan:

```bash
bash code/evaluate_bams_for_variant_calling.sh
# Results: ~/results/bam_evaluation/variant_calling_readiness.tsv
```

---

## 11. Push Results to GitHub

After a run completes, push reports, tables, figures, and HTML outputs to the repo:

```bash
# Auto-detect latest local run
bash push_results.sh

# Push a specific run
bash push_results.sh res_20260311_184352

# Pull from S3 first, then push (useful on a fresh instance)
bash push_results.sh res_20260311_184352 --from-s3
```

**What gets committed:**

| Included | Excluded |
|---|---|
| `.html`, `.pdf`, `.png`, `.svg` | `.bam`, `.bai`, `.vcf.gz` |
| `.tsv`, `.csv`, `.txt`, `.xlsx` | `.gz`, `.zip`, `.tar.gz` |
| `.vcf` (filtered/annotated) | `preprocessed_bams/` |
| `.md`, `.R`, `.py`, `.sh` | `*.recal.table`, `*.tmp` |

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `S3 access denied` | Run `aws configure` and check IAM permissions for `s3://neoantigen2026-rerun` |
| `GATK not found` | Run `bash code/setup_pipeline.sh` |
| `netMHCpan not found` | Manual install required — see Step 3 note above |
| `BQSR warning: no known SNPs` | Check `~/ref/mm10/mgp_snps.vcf.gz` exists; re-run `setup_pipeline.sh` |
| Step fails mid-run | Re-run with `bash code/run_pipeline.sh <STEP>` after fixing the issue |
