# NeoAntigen Pipeline — Commands Log

## Environment

| Item | Value |
|------|-------|
| Instance | EC2 (ca-central-1) |
| OS | Amazon Linux 2023 |
| S3 Bucket | `s3://bam-wes/NeoAntigen-aws/` |
| Working Dir | `/home/ec2-user/` |
| Code Dir | `/home/ec2-user/code/` |
| Results Dir | `/home/ec2-user/results/` |
| Ref Dir | `/home/ec2-user/ref/mm10/` |

---

## Tool Versions

| Tool | Version | Install Method |
|------|---------|---------------|
| samtools | 1.23 | `mamba install -c bioconda samtools` |
| bcftools | 1.23 | `mamba install -c bioconda bcftools` |
| GATK4 | 4.6.2.0 | `mamba install -c bioconda gatk4` |
| R | 4.3.2 | `sudo dnf install -y R` |
| Java (OpenJDK) | 25.0.2 | pre-installed |
| BWA-MEM | 0.7.17 | used for original alignment (reference) |

---

## Initial Setup

### 1. Install Miniforge + bioinformatics tools
```bash
# Install Miniforge
curl -sLO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3

# Install samtools and bcftools
~/miniforge3/bin/mamba install -y -c bioconda samtools bcftools

# Install GATK4
~/miniforge3/bin/mamba install -y -c bioconda gatk4

# Fix GATK python shebang (Amazon Linux 2023 uses python3, not python)
sed -i 's|#!/usr/bin/env python|#!/usr/bin/env python3|' ~/miniforge3/bin/gatk

# Set GATK jar path (add to ~/.bashrc for persistence)
export GATK_LOCAL_JAR=~/miniforge3/share/gatk4-4.6.2.0-0/gatk-package-4.6.2.0-local.jar
```

### 2. Install R and system dependencies
```bash
sudo dnf install -y R
sudo dnf install -y harfbuzz-devel fribidi-devel libcurl-devel openssl-devel libxml2-devel
```

> R packages (ggplot2, dplyr, tidyr, readr, patchwork, RColorBrewer, scales, ggrepel)
> are installed automatically on first run into `~/R/library`.

---

## Step 1 — BAM Validation

**Script:** `code/validate_bams.sh`
**R report:** `code/summarize_bam_qc.R`
**Input:** `s3://bam-wes/NeoAntigen-aws/data/bam/*.bam` (7 files, ~48.9 GB total)
**Output:** `results/bam_validation/`

```bash
bash /home/ec2-user/code/validate_bams.sh
```

Checks performed per BAM:

| Check | Tool | Method |
|-------|------|--------|
| BGZF EOF integrity | AWS CLI + xxd | `aws s3api get-object --range` last 28 bytes |
| BAM header / sort order | samtools view -H | First 10 MB via `aws s3api get-object --range` |
| Reference genome build | samtools view -H | Parse `@PG`/`@SQ` tags for mm10/hg38/hg19 |
| Read groups | samtools view -H | Count `@RG` lines in header |
| BAM index (.bai) | AWS CLI | `aws s3 ls *.bam.bai` |
| Read stats | samtools flagstat | Full BAM streamed: `aws s3 cp ... - \| samtools flagstat -` |

### Results (2026-03-10)

| Sample | Size | Mapped % | Sort | Reference | BAI | QC | Notes |
|--------|------|----------|------|-----------|-----|----|-------|
| 34_D52_old | 4.25 GB | 100% | coordinate | mm10/GRCm38 | YES | PASS | Valid BGZF; indexed |
| 36_D99_new | 9.88 GB | 100% | coordinate | mm10/GRCm38 | YES | PASS | Valid BGZF; indexed |
| 38_D99_new | 9.39 GB | 100% | coordinate | mm10/GRCm38 | YES | PASS | Valid BGZF; indexed |
| 423_D0_old | 4.71 GB | 100% | coordinate | mm10/GRCm38 | YES | PASS | Valid BGZF; indexed |
| 424_D0__old | 4.07 GB | 100% | coordinate | mm10/GRCm38 | YES | PASS | Valid BGZF; indexed |
| 428_D20_new | 9.49 GB | 100% | coordinate | mm10/GRCm38 | YES | PASS | Valid BGZF; indexed |
| 42_D122_old | 3.78 GB | 100% | coordinate | mm10/GRCm38 | YES | PASS | Valid BGZF; indexed |

**Key findings:**
- All 7 BAMs aligned with BWA-MEM 0.7.17 to mm10/GRCm38 (UCSC chr-prefix), coordinate-sorted, 100% mapping rate.
- Initial BGZF EOF false-negatives caused by S3 regional redirect — resolved by using `aws s3api` with explicit bucket/key.
- No `.bai` index files existed — created in Step 2.

### Re-run R report only
```bash
Rscript /home/ec2-user/code/summarize_bam_qc.R \
  /home/ec2-user/results/bam_validation/bam_validation_results.tsv \
  /home/ec2-user/results/bam_validation
```

### Output files

| File | Description |
|------|-------------|
| `results/bam_validation/bam_validation_results.tsv` | Raw per-sample QC table |
| `results/bam_validation/bam_qc_summary_table.tsv` | Formatted summary table |
| `results/bam_validation/bam_qc_summary_table.html` | Colour-coded HTML table |
| `results/bam_validation/bam_qc_report.pdf` | QC figures (PDF) |
| `results/bam_validation/bam_qc_report.png` | QC figures (PNG) |
| `results/bam_validation/flagstat/*.flagstat.txt` | Per-sample flagstat output |

---

## Step 2 — BAM Indexing

**Script:** `code/index_bams.sh`
**Input:** `s3://bam-wes/NeoAntigen-aws/data/bam/*.bam`
**Output:** `.bai` files uploaded back to `s3://bam-wes/NeoAntigen-aws/data/bam/`

```bash
bash /home/ec2-user/code/index_bams.sh
```

Strategy: download one BAM → `samtools index -@ 4` → upload `.bai` to S3 → delete local BAM → repeat.

### Results (2026-03-10)

| Sample | Size | Index Time | Status |
|--------|------|-----------|--------|
| 34_D52_old | 4.25 GB | 0.5 min | PASS |
| 36_D99_new | 9.88 GB | 1.1 min | PASS |
| 38_D99_new | 9.39 GB | 1.1 min | PASS |
| 423_D0_old | 4.71 GB | 0.6 min | PASS |
| 424_D0__old | 4.07 GB | 0.5 min | PASS |
| 428_D20_new | 9.49 GB | 1.1 min | PASS |
| 42_D122_old | 3.78 GB | 0.5 min | PASS |

---

## Step 3 — Variant Calling (GATK4 Mutect2)

**Scripts:** `code/setup_variant_calling.sh`, `code/variant_calling.sh`
**R report:** `code/summarize_variants.R`
**Reference:** mm10/GRCm38 — UCSC (chr-prefixed), downloaded from `hgdownload.soe.ucsc.edu`
**Output:** `results/variant_calling/`  +  `s3://bam-wes/NeoAntigen-aws/results/variant_calling/`

### Tumor-Normal Pairing

| Role | Sample | SM Tag | Timepoint |
|------|--------|--------|-----------|
| Normal (matched) | 423_D0_old | D0_old | D0 (baseline) |
| PoN normal 1 | 423_D0_old | D0_old | D0 |
| PoN normal 2 | 424_D0__old | D0_old | D0 |
| Tumor | 428_D20_new | D20_new | D20 |
| Tumor | 34_D52_old | D52_old | D52 |
| Tumor | 36_D99_new | D99_new | D99 |
| Tumor | 38_D99_new | D99_new | D99 |
| Tumor | 42_D122_old | D122_old | D122 |

### Step 3a — Setup (run once)

```bash
bash /home/ec2-user/code/setup_variant_calling.sh
```

What it does:
1. Installs GATK4 via mamba (`-c bioconda gatk4`)
2. Fixes GATK python shebang: `sed -i 's|python|python3|' ~/miniforge3/bin/gatk`
3. Sets `GATK_LOCAL_JAR` environment variable
4. Downloads mm10.fa.gz (~800 MB) from UCSC and decompresses with `gunzip`
5. Downloads mm10.dict and mm10.fa.fai from S3

**Reference files after setup:**
```
/home/ec2-user/ref/mm10/
├── mm10.fa      (2.6 GB — UCSC mm10 genome)
├── mm10.fa.fai  (FASTA index)
└── mm10.dict    (sequence dictionary)
```

> **Note:** Use `gunzip` (not `bgzip`) to decompress mm10.fa.gz on Amazon Linux 2023.
> `bgzip` is available at `~/miniforge3/bin/bgzip` if needed for other steps.

### Step 3b — Run Variant Calling

```bash
bash /home/ec2-user/code/variant_calling.sh
```

Pipeline steps per sample:

| Step | GATK Tool | Purpose |
|------|-----------|---------|
| 3a | `Mutect2` (tumor-only) | Call variants on each normal for PoN |
| 3a | `GenomicsDBImport` | Consolidate normal VCFs into genomics DB |
| 3a | `CreateSomaticPanelOfNormals` | Build PoN VCF from normal DB |
| 3b | `Mutect2` (tumor-normal) | Somatic SNV/indel calling with PoN |
| 3c | `FilterMutectCalls` | Apply somatic filters → PASS variants |
| 3d | `bcftools stats` | Count SNVs vs indels in PASS calls |
| 3d | `summarize_variants.R` | Summary table + figures |

Disk management: each tumor BAM is downloaded, processed, and deleted before the next.
VCFs and stats are uploaded to S3 immediately after each sample.

### Re-run R report only
```bash
Rscript /home/ec2-user/code/summarize_variants.R \
  /home/ec2-user/results/variant_calling/variant_calling_results.tsv \
  /home/ec2-user/results/variant_calling
```

### Output files

| File | Description |
|------|-------------|
| `results/variant_calling/variant_calling_results.tsv` | Per-sample variant counts (raw) |
| `results/variant_calling/variant_summary_table.tsv` | Formatted summary table |
| `results/variant_calling/variant_summary_table.html` | Colour-coded HTML table |
| `results/variant_calling/variant_calling_report.pdf` | QC figures — 4 panels (PDF) |
| `results/variant_calling/variant_calling_report.png` | QC figures — 4 panels (PNG) |
| `results/variant_calling/pon/pon.vcf.gz` | Panel of Normals VCF |
| `results/variant_calling/vcf_raw/*.vcf.gz` | Raw Mutect2 VCFs (also on S3) |
| `results/variant_calling/vcf_filtered/*.filtered.vcf.gz` | Filtered PASS VCFs (also on S3) |
| `results/variant_calling/stats/*.mutect2.stats` | Mutect2 stats per sample |

---

## Step 4 — SnpEff Annotation

**Script:** `code/steps/04_annotate.sh`
**Input:** `results/variant_calling/vcf_filtered/*.filtered.vcf.gz`
**Output:** `results/annotated/*.annotated.vcf.gz`  +  S3

```bash
bash /home/ec2-user/code/steps/04_annotate.sh
```

What it does:
- Downloads filtered VCFs from S3 (if not local)
- Runs SnpEff `GRCm38.86` database (`-cancer` flag, HTML summary stats)
- bgzips + tabix-indexes output VCFs
- Uploads annotated VCFs and HTML stats to `s3://bam-wes/NeoAntigen-aws/results/annotated/`

SnpEff setup (run once):
```bash
mamba install -y -c bioconda snpeff
snpEff download GRCm38.86
```

---

## Step 5 — PyClone Input Preparation

**Script:** `code/steps/05_pyclone_prep.py`
**Input:** `results/annotated/*.annotated.vcf.gz`
**Output:** `results/pyclone/input/*.yaml` + `*.mutations.tsv`  +  S3

```bash
python3 /home/ec2-user/code/steps/05_pyclone_prep.py
```

What it does:
- Parses SnpEff-annotated VCFs (PASS variants only)
- Filters: depth ≥ 10, alt_count ≥ 3, SNVs/short indels only
- Writes one PyClone YAML per tumour sample (ref_counts, var_counts, normal_cn=2)
- Writes merged `all_samples_mutations.tsv`

---

## Step 6 — PyClone MCMC

**Script:** `code/steps/06_pyclone_run.sh`
**Input:** `results/pyclone/input/*.yaml`
**Output:** `results/pyclone/output/tables/loci.tsv`, `cluster.tsv`  +  S3

```bash
bash /home/ec2-user/code/steps/06_pyclone_run.sh
```

PyClone setup (run once):
```bash
pip install pyclone
```

Configuration:
- Model: beta-binomial
- Iterations: 10,000 (burn-in: 1,000, thin: 1)
- Init method: disconnected
- tumour_content: 1.0 (update if purity estimates available)

---

## Step 7 — PyClone Tables & Figures

**Script:** `code/steps/07_pyclone_tables.R`
**Input:** `results/pyclone/output/tables/loci.tsv`, `cluster.tsv`
**Output:** `results/pyclone/tables/`  +  S3

```bash
Rscript /home/ec2-user/code/steps/07_pyclone_tables.R
```

Outputs:
| File | Description |
|------|-------------|
| `cluster_summary.tsv` | Mutations per cluster |
| `cluster_prevalence.tsv` | Mean CP ± SD per cluster per timepoint |
| `loci_clusters.tsv` | Mutation → cluster assignments |
| `pyclone_report.pdf/png` | 3-panel figure: sizes, evolution, VAF distributions |

---

## Mutect2 Fix — Exon Interval List

The initial Mutect2 run was killed (exit 144) because it ran over the full mouse genome without an exon interval restriction.

**Fix applied (2026-03-10):**
```bash
# 1. Generate merged exon BED from GRCm38.86 GTF (streaming — no full download)
curl -s "https://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.chr.gtf.gz" \
  | zcat \
  | awk '$3=="exon" {
      chrom = ($1 ~ /^[0-9XYMm]/) ? "chr"$1 : $1
      start = $4 - 1
      end   = $5
      if (end > start) print chrom"\t"start"\t"end
  }' \
  | sort -k1,1 -k2,2n \
  | ~/miniforge3/bin/bedtools merge -i - \
  > /home/ec2-user/ref/mm10/mm10_exons.bed
# Result: 261,746 merged exon intervals (6 MB)

# 2. Remove truncated PoN VCF
rm -f /home/ec2-user/results/variant_calling/pon/423_D0_old.vcf.gz
```

`02_mutect2.sh` was updated to:
- Add `-L ${EXON_BED} --interval-padding 100` to all Mutect2 and GenomicsDBImport calls
- Fall back to raw S3 BAMs if Step 1 preprocessed BAMs are not available

---

## Step 8 — Peptide Extraction

**Script:** `code/steps/08_peptide_extract.py`
**Input:** `results/annotated/*.annotated.vcf.gz`
**Output:** `results/peptides/*.peptides.tsv`, `peptides_{8,9,10}mer.fasta`  +  S3

```bash
python3 /home/ec2-user/code/steps/08_peptide_extract.py
```

Variant effects processed:
- `missense_variant`, `inframe_insertion`, `inframe_deletion`, `stop_gained`,
  `stop_lost`, `start_lost`, `disruptive_inframe_deletion/insertion`  → sliding-window k-mers
- `frameshift_variant`  → junction peptides (WT prefix + novel suffix)

Uses pyensembl (GRCm38 release 86) for accurate protein sequences; falls back to
simplified extraction if pyensembl is unavailable.

---

## Step 9 — netMHCpan Binding Prediction (previously Step 9)

**Script:** `code/steps/09_netmhcpan.sh`
**Input:** `results/peptides/peptides_{8,9,10}mer.fasta`
**Output:** `results/netmhcpan/`  +  S3

```bash
bash /home/ec2-user/code/steps/09_netmhcpan.sh
```

Alleles: `H-2-Kb`, `H-2-Db`  (C57BL/6 H-2b haplotype)
Lengths: 8-mer, 9-mer, 10-mer

Thresholds:
| Class | EL Rank |
|-------|---------|
| Strong Binder (SB) | < 0.5% |
| Weak Binder (WB) | < 2.0% |
| Non-Binder (NB) | ≥ 2.0% (excluded) |

netMHCpan setup (run once):
```bash
mamba install -y -c bioconda netmhcpan
```

Outputs:
| File | Description |
|------|-------------|
| `{plen}mer_{allele}.tsv` | Per-length per-allele binders (SB+WB only) |
| `all_predictions.tsv` | Merged table — all alleles, all lengths |

---

## Step 8b — Fusion Peptide Extraction

**Script:** `code/steps/08b_fusion_peptides.py`
**Input:** `s3://bam-wes/NeoAntigen-aws/data/FusionCalling/arriba_files/Arriba_*.txt`
**Output:** appended to `results/peptides/peptides_{8,9,10}mer.fasta` + `results/peptides/fusions_all.tsv`

```bash
python3 /home/ec2-user/code/steps/08b_fusion_peptides.py
```

Data sources:
- 8 Arriba files per timepoint: D20_new, D52_old, D88_old, D99_old, D99_new_36, D99_new_38, D109_new, D122_old
- Filters: `confidence` ∈ {high, medium}
- `peptide_sequence` column — `|` marks fusion junction, `___` marks stop/frameshift
- Extracts 8/9/10-mer junction-spanning peptides (all k-mers with junction in window)
- In-frame fusions → missense-like; out-of-frame → frameshift-like
- Appends to existing per-length FASTAs so netMHCpan sees both WES + fusion peptides

---

## Step 11 — RNASeq Expression Integration

**Script:** `code/steps/11_rnaseq_expression.R`
**Input:** `s3://bam-wes/NeoAntigen-aws/data/RNASeq/count_files/`
**Output:** `results/rnaseq/`  +  S3

```bash
Rscript /home/ec2-user/code/steps/11_rnaseq_expression.R
```

RNASeq samples (10 timepoints):

| Column | Sample ID | Timepoint | Batch |
|--------|-----------|-----------|-------|
| 0 Day_0_Old | 423_D0_old | D0 | old |
| 20 S428_Day_20_New | 428_D20_new | D20 | new |
| 21 S443_Day_21_New | 443_D21_new | D21 | new |
| 52 S34_Day_52_Old | 34_D52_old | D52 | old |
| 88 S80_Day_88_Old | D88_old | D88 | old |
| 99 S32_Day_99_Old-1 | D99_old | D99 | old |
| 99 S36_Day_99_New-2 | 36_D99_new | D99 | new |
| 99 S38_Day_99_New-3 | 38_D99_new | D99 | new |
| 109 S2661_Day_109_New | D109_new | D109 | new |
| 122 S42_Day_122_Old | 42_D122_old | D122 | old |

Batch correction: `limma::removeBatchEffect(log2(counts + 0.5), batch = old/new)`

Outputs:
| File | Description |
|------|-------------|
| `expression_matrix.tsv` | Wide: genes × samples, batch-corrected log2 |
| `expression_by_sample.tsv` | Long format for joining to candidates |
| `rnaseq_report.pdf/png` | PCA before/after, density, top-50 variable gene heatmap |

---

## Step 10 — Neoantigen Candidate Filtering & Ranking

**Script:** `code/steps/10_filter_results.R`
**Input:** `results/netmhcpan/all_predictions.tsv`,
           `results/peptides/all_peptides.tsv`,
           `results/peptides/fusions_all.tsv`,
           `results/pyclone/tables/loci_clusters.tsv`,
           `results/rnaseq/expression_matrix.tsv`
**Output:** `results/final/`  +  S3

```bash
Rscript /home/ec2-user/code/steps/10_filter_results.R
```

Priority scoring (max = 7):

| Component | Points | Criteria |
|-----------|--------|---------|
| Binding SB | 3 | EL Rank < 0.5% |
| Binding WB | 1 | EL Rank < 2.0% |
| Missense | 2 | SNV amino acid substitution |
| Frameshift/fusion | 1 | Novel junction peptide |
| Clonal | 1 | max_CP ≥ 0.5 (PyClone) |
| Expressed | 1 | median log2 expr > 1.0 in tumour samples |

Outputs:
| File | Description |
|------|-------------|
| `neoantigen_candidates.tsv` | All SB+WB binders with full annotation |
| `neoantigen_top_ranked.tsv` | SB only, ranked by priority score |
| `binding_summary.tsv` | Counts by allele/length/class/mut_type |
| `neoantigen_report.pdf/png` | 6-panel: type breakdown, EL ranks, clonality×expression, scatter, top-20, priority distribution |

---

## Run Full Pipeline

```bash
# Run all 12 steps
bash /home/ec2-user/code/run_pipeline.sh

# Run from a specific step
bash /home/ec2-user/code/run_pipeline.sh 5       # steps 5-12

# Run a range
bash /home/ec2-user/code/run_pipeline.sh 3 7     # steps 3-7
```

Each step checks S3 for existing outputs and skips if already complete.

---

## GitHub Repository

```
https://github.com/jahanshah/NeoAntigen2026-aws
```

### Push updated code and results
```bash
cd /home/ec2-user/code
git add -A
git commit -m "describe changes"
git push origin main
```
