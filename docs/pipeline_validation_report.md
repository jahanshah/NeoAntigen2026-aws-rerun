# NeoAntigen2026 Pipeline Validation Report

**Analysis date:** 2026-03-16 (updated)
**Run ID:** res_20260311_225555
**Platform:** AWS EC2 (GRCm38/mm10, 6 tumor samples)
**Tumor samples:** 443_D21_new, 428_D20_new, 34_D52_old, 36_D99_new, 38_D99_new, 42_D122_old
**Normal sample:** 423_D0_old (matched); additional PoN normal: 451_D0_old
**MHC allele:** H-2-Db (IEDB REST, netMHCpan_EL, 8/9/10-mer)
**Binding thresholds:** SB < 0.5% rank_EL · WB < 2.0% · PB < 10.0%

---

## Section 1: Sequence Mapping and Preprocessing

### 1a. WES — Preparing BAM Files for Variant Calling

**Script:** `code/steps/01_bam_preprocess.sh`
**Purpose:** Convert raw WES FASTQs to analysis-ready BAMs suitable for GATK4 Mutect2 somatic calling.

**Pipeline:**
1. **Sort:** `samtools sort` (coordinate order)
2. **Mark Duplicates:** GATK4 `MarkDuplicates` (PCR/optical duplicates flagged, not removed)
3. **Base Quality Score Recalibration (BQSR):** GATK4 `BaseRecalibrator` + `ApplyBQSR` using known mouse SNPs
4. **Validation:** GATK4 `ValidateSamFile`

**Key commands:**
```bash
# Mark duplicates
gatk MarkDuplicates \
    -I sample.sorted.bam \
    -O sample.dedup.bam \
    -M sample.metrics.txt \
    --REMOVE_DUPLICATES false \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

# BQSR (scattered per chromosome for speed, then merged)
gatk BaseRecalibrator \
    -R mm10.fa -I sample.dedup.bam \
    -L chrN --known-sites mm10_mgp_snps.vcf.gz \
    -O sample.chrN.recal.table

gatk ApplyBQSR \
    -R mm10.fa -I sample.dedup.bam \
    --bqsr-recal-file sample.recal.table \
    -O sample.preproc.bam
```

**Tools and versions:**
| Tool | Version | Purpose |
|------|---------|---------|
| GATK4 | 4.x | MarkDuplicates, BaseRecalibrator, ApplyBQSR |
| samtools | 1.x | Sort, index |
| Reference | mm10 (GRCm38) | Genome reference |
| Known sites | Mouse Genome Project SNPs (GRCm38) | BQSR known variants |

**Output:** `*.preproc.bam` — stored permanently at `s3://neoantigen2026-rerun/data/bam/wes/preprocessed/`

---

### 1b. RNA-seq — Preparing BAM Files for Expression and Fusion Calling

**Scripts:** `code/rnaseq/steps/R02_align.sh` (expression), `code/rnaseq/steps/run_arriba.sh` (fusion)
**Purpose:** Align RNA-seq reads to mm10 for (a) gene expression quantification and (b) fusion gene detection.

**Expression alignment (R02):**
```bash
# STAR 2-pass alignment with GeneCounts
STAR \
    --runThreadN 8 \
    --genomeDir /ref/mm10/star_index \
    --readFilesIn R1.fastq.gz R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --twopassMode Basic \
    --quantMode GeneCounts \
    --outFileNamePrefix sample/

samtools index sample/Aligned.sortedByCoord.out.bam
```

**Fusion alignment (Arriba requires chimeric read mode):**
```bash
# STAR with chimeric read flags (required by Arriba 2.5.1)
STAR-avx2 \
    --runThreadN 8 \
    --genomeDir /ref/mm10/star_index \
    --readFilesIn R1.fastq.gz R2.fastq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --chimSegmentMin 10 \
    --chimOutType WithinBAM SoftClip \
    --chimJunctionOverhangMin 10 \
    --chimScoreMin 1 \
    --chimMultimapNmax 50 \
    --peOverlapNbasesMin 10 \
    --outFileNamePrefix sample_arriba/
```

**Tools and versions:**
| Tool | Version | Purpose |
|------|---------|---------|
| STAR | 2.7.x (STAR-avx2) | RNA-seq alignment (2-pass) |
| samtools | 1.x | BAM indexing |
| Reference | mm10.fa + mm10_chr.gtf | Genome + annotation |
| STAR index | GRCm38/mm10 | Pre-built genome index |

**Samples processed:** 10 RNA-seq samples (6 tumor timepoints + 4 additional)
**Output:** `*.Aligned.sortedByCoord.out.bam` + `ReadsPerGene.out.tab` per sample

---

## Section 2: Variant Calling

**Scripts:** `code/steps/02_mutect2.sh`, `code/steps/03_filter_vcf.sh`
**Purpose:** Detect somatic SNVs and indels in each tumor sample relative to the matched normal.

### 2a. Panel of Normals (PoN)

```bash
# Step 1: Tumor-only calls on each normal
gatk Mutect2 \
    -R mm10.fa -I normal.preproc.bam \
    --tumor-sample NORMAL_ID \
    -O normal.pon.vcf.gz

# Step 2: Combine normals into GenomicsDB
gatk GenomicsDBImport \
    -V 423_D0_old.pon.vcf.gz -V 451_D0_old.pon.vcf.gz \
    --genomicsdb-workspace-path pon_db \
    -L exon_targets.bed

# Step 3: Build PoN
gatk CreateSomaticPanelOfNormals \
    -V gendb://pon_db \
    -O panel_of_normals.vcf.gz
```

### 2b. Somatic Calling (Tumor vs. Normal)

```bash
# Mutect2 tumor-vs-normal (with orientation model)
gatk Mutect2 \
    -R mm10.fa \
    -I tumor.preproc.bam -tumor TUMOR_ID \
    -I normal.preproc.bam -normal NORMAL_ID \
    --panel-of-normals panel_of_normals.vcf.gz \
    --germline-resource gnomad_mm10.vcf.gz \
    -L exon_targets.bed --interval-padding 100 \
    --f1r2-tar-gz tumor.f1r2.tar.gz \
    -O tumor.raw.vcf.gz

# Learn orientation model (OxoG / FFPE bias)
gatk LearnReadOrientationModel \
    -I tumor.f1r2.tar.gz \
    -O tumor.read-orientation-model.tar.gz

# Contamination estimate
gatk GetPileupSummaries \
    -I tumor.preproc.bam -V gnomad_mm10.vcf.gz \
    -L common_biallelic.vcf.gz \
    -O tumor.pileup.table
gatk CalculateContamination \
    -I tumor.pileup.table -matched normal.pileup.table \
    -O tumor.contamination.table
```

### 2c. Filtering

```bash
# FilterMutectCalls (applies all Mutect2 filters)
gatk FilterMutectCalls \
    -R mm10.fa -V tumor.raw.vcf.gz \
    --stats tumor.stats \
    --orientation-bias-artifact-priors tumor.read-orientation-model.tar.gz \
    --contamination-table tumor.contamination.table \
    -O tumor.filtered.vcf.gz

# Restrict to exon targets (PASS calls only)
gatk SelectVariants \
    -V tumor.filtered.vcf.gz \
    -L exon_targets.bed --interval-padding 100 \
    --select-type-to-include SNP \
    --select-type-to-include INDEL \
    -select "FILTER == 'PASS'" \
    -O tumor.filtered.exon.vcf.gz
```

### 2d. Annotation

**Script:** `code/steps/04_annotate.sh`

```bash
# SnpEff functional annotation (Ensembl 102 / GRCm38.86)
snpeff \
    -Xmx8g -v GRCm38.86 \
    -stats tumor.snpeff_stats.html \
    -cancer -noLog \
    tumor.filtered.exon.vcf.gz \
    | bgzip -c > tumor.annotated.vcf.gz

bcftools index -t tumor.annotated.vcf.gz
```

**Tools and versions:**
| Tool | Version | Purpose |
|------|---------|---------|
| GATK4 | 4.x | Mutect2, FilterMutectCalls, PoN creation |
| SnpEff | 5.x | Functional annotation (ANN field) |
| bcftools | 1.x | VCF indexing, manipulation |
| Database | GRCm38.86 (SnpEff) | Variant effect annotation |

**Results:**
- 6 tumor samples called against matched normal (423_D0_old)
- PoN built from 2 normal samples (423_D0_old, 451_D0_old)
- FilterMutectCalls applied read-orientation artifact correction
- PASS + exon-restricted VCFs used for downstream analysis

---

## Section 3: Count Estimation (Gene Expression)

**Script:** `code/rnaseq/steps/R03_quantify.sh`
**Purpose:** Generate a gene-level count matrix from RNA-seq BAMs for differential expression and neoantigen scoring.

**Strandedness confirmation:**
STAR `--quantMode GeneCounts` was used during alignment and produced three count columns per sample: unstranded, forward-stranded, and reverse-stranded. Comparing columns for sample 443_D21_new confirmed forward-stranded library:
- Forward-strand counts: ~8.8M assigned reads
- Reverse-strand counts: ~527K assigned reads (~6% of forward)

This confirmed `-s 1` (forward-stranded) must be used in featureCounts. The previous pipeline run used `-s 2` (reverse-stranded), which assigned only ~6% of reads — causing featureCounts to fail silently and fall back to STAR GeneCounts.

**Primary method: featureCounts with `-s 1` (forward-stranded)**
```bash
featureCounts \
    -a mm10_chr.gtf \
    -o raw_counts_s1.txt \
    -T 8 \
    -p --countReadPairs \
    -s 1 \
    -t exon \
    -g gene_id \
    sample1.bam sample2.bam ... sample6.bam
```

Output: `raw_counts_s1.txt` — 48,709 genes × 6 tumor samples; assignment rate 36.8–42.1% (expected for RNA-seq with intronic reads excluded).

**Previous fallback (no longer used): STAR GeneCounts merge**
The original R03_quantify.sh fell back to merging `ReadsPerGene.out.tab` files using column 4 (parts[3] = reverse-stranded) — incorrect for this forward-stranded library. Manual regeneration in the interim used the unstranded column (~2% difference from forward-stranded). Both have been superseded by featureCounts `-s 1`.

**Tools and versions:**
| Tool | Version | Purpose |
|------|---------|---------|
| STAR | 2.7.x | `--quantMode GeneCounts` (fallback) |
| featureCounts (subread) | 2.x | Primary count estimation |
| GTF | mm10_chr.gtf (Ensembl 102) | Gene annotation |
| Python/pandas | 3.13 / 2.x | Count matrix merging |

**Results:**
- 48,709 genes × 10 samples
- CPM normalization applied in Step 10 for neoantigen scoring
- Ensembl gene IDs mapped to gene symbols via `Mus_musculus.GRCm38.102.gtf.gz` (55,487 genes)

---

## Section 4: Fusion Calling

**Script:** `code/rnaseq/steps/run_arriba.sh`
**Purpose:** Detect gene fusion transcripts from RNA-seq using Arriba (requires chimeric read-aware STAR alignment).

**Note:** The original R02 alignment omitted chimeric read flags (`--chimSegmentMin` etc.) required by Arriba. The dedicated `run_arriba.sh` script re-aligns with correct flags before calling Arriba.

```bash
# Arriba fusion calling
arriba \
    -x sample/Aligned.sortedByCoord.out.bam \
    -a mm10.fa \
    -g mm10_chr.gtf \
    -b blacklist_mm10_GRCm38_v2.5.1.tsv.gz \
    -k known_fusions_mm10_GRCm38_v2.5.1.tsv.gz \
    -o Arriba_sample.txt \
    -O Arriba_sample.discarded.txt
```

**Tools and versions:**
| Tool | Version | Purpose |
|------|---------|---------|
| Arriba | 2.5.1 | Fusion detection |
| STAR-avx2 | 2.7.x | Chimeric-mode alignment (required by Arriba) |
| Blacklist | blacklist_mm10_GRCm38_v2.5.1.tsv.gz | Arriba artifact filter |
| Known fusions | known_fusions_mm10_GRCm38_v2.5.1.tsv.gz | Arriba confidence boost |

**Three-stage fusion peptide recovery:**

**Stage 1 — Arriba junction peptides (08b_fusion_peptides.py):**
- All 6 WES-matched tumor samples + 3 RNA-only timepoints (D88_old, D99_old, D109_new) processed
- Fusions filtered to confidence = high or medium
- 35 high/medium-confidence fusions across 6 WES samples
- 28/35 fusions: `peptide_sequence = "."` (see Stage 3)
- 7/35 fusions: 102 junction peptides extracted
- **RNA-only samples (D88_old, D99_old, D109_new):** Re-aligned from FASTQs with STAR chimeric mode (`--chimSegmentMin 10 --chimOutType WithinBAM SoftClip`); Arriba run completed — no novel somatic fusions found (Adgrf1::Adgrf5 germline SV + Ig VDJ recombination artifacts only)
- Arriba lowercase convention: lowercase letters = out-of-frame translation after junction. Peptides converted to uppercase before IEDB submission (biologically valid — these are real amino acids from the frameshifted reading frame; IEDB requires uppercase input)
- Germline filter: Adgrf1::Adgrf5 found in 7/9 tumor samples (6 WES + 3 RNA-only) — excluded as germline structural variant (threshold: ≥ 3 samples)

**Stage 2 — Dot-fusion ORF translation (08d_fusion_orf.py):**
For fusions where Arriba outputs `"."` with a downstream CDS breakpoint, the downstream gene CDS was reconstructed from the Ensembl GTF and mm10 genome FASTA and translated from the breakpoint using Biopython:
```python
# pysam.FastaFile → genomic sequence at breakpoint
# Ensembl GRCm38 R102 GTF → CDS exon coordinates
# Bio.Seq.translate → protein from breakpoint to stop codon
```
5 non-germline "." fusions processed → 12 translated peptides (4 fusions × 3 lengths)

**Stage 3 — STAR-Fusion 1.13.0 (installed, planned):**
STAR-Fusion 1.13.0 installed as complementary fusion caller. The existing STAR alignment BAMs were not generated with `--chimSegmentMin` flags, so STAR-Fusion would require re-alignment from FASTQs (available at `s3://neoantigen2026-rerun/data/fastq/rnaseq/`). Not yet executed.

**Fusion IEDB results (H-2-Db, netMHCpan_EL):**
| Fusion | Stage | Peptide | rank_EL% | Class |
|--------|-------|---------|----------|-------|
| Fxr1::Zfp704 | 08b | AFYKNSMKV | 0.98% | WB |
| Fxr1::Zfp704 | 08b | GAFYKNSMKV | 0.98% | WB |
| Fxr1::Zfp704 | 08b | AFYKNSMKVM | 2.00% | WB |
| Zfp740::Csad | 08d | MADSKPLRTL | 1.50% | WB |
| Fxr1::Zfp704 | 08b | YKNSMKVMF | 4.20% | PB |
| Zfp740::Csad | 08d | MADSKPLRT | 9.80% | PB |
| Trp53::Sat2 | 08d | ELAEFEKL | — | H-2-Kb WB 1.7% |

---

## Section 5: Peptide Identification

**Scripts:** `code/steps/08_peptide_extract.py`, `code/steps/08b_fusion_peptides.py`, `code/steps/08c_frameshift_realseq.py`

### 5a. Somatic Mutation Peptides

**Approach:** Extract mutant peptides from SnpEff-annotated VCFs using Ensembl GRCm38 release 102 protein sequences.

```python
# code/steps/08_peptide_extract.py

# Load Ensembl protein FASTA (gene_symbol → protein sequences)
# File: ~/.cache/pyensembl/GRCm38/ensembl102/Mus_musculus.GRCm38.pep.all.fa.gz
protein_db = load_protein_db()  # 22,547 gene symbols

# For each PASS+HIGH/MODERATE variant in annotated VCF:
# Parse ANN field → gene, HGVSp, effect type
# Missense/inframe: apply amino acid substitution to WT protein
#   → extract k-mers (k=8,9,10) spanning the mutated position
# Frameshift: see Step 5b below (08c_frameshift_realseq.py)

# Fallback: pyensembl 2.3.13 incompatible with polars 0.19.19 on Python 3.13
# (Rust panic: "Invalid strand: 2" — uncatchable by except Exception)
# Solution: except BaseException → fallback to protein FASTA lookup
try:
    from pyensembl import EnsemblRelease
    data = EnsemblRelease(102, species="mouse")
    data.download(); data.index()
except BaseException as e:
    # Fallback to Ensembl protein FASTA extraction
    return parse_vcf_for_peptides(vcf_path, sample)
```

**Mutation types included:**
- `missense_variant` — amino acid substitution
- `inframe_insertion` / `inframe_deletion` / `disruptive_inframe_*` — in-frame indels
- `stop_gained` / `stop_lost` / `start_lost`
- `frameshift_variant` — see 5b

**Filter criteria:** FILTER=PASS, impact=MODERATE or HIGH

### 5b. Frameshift Novel Sequence Translation (proper bioinformatics approach)

**Script:** `code/steps/08c_frameshift_realseq.py`
**Tools used:** Biopython 1.86 (`Bio.Seq.translate`), Ensembl cDNA FASTA

This step replaces the placeholder `AAAAAAAAAAAAAAAA` suffix used in the initial extraction with actual translated novel sequences:

```python
# code/steps/08c_frameshift_realseq.py

from Bio.Seq import Seq

# 1. Load Ensembl cDNA FASTA for target genes
#    File: Mus_musculus.GRCm38.cdna.all.fa.gz (pyensembl cache)
cdna_db = load_cdna_db(target_genes)  # gene_symbol → [(transcript_id, biotype, sequence)]

# 2. For each frameshift variant (gene, HGVSc, HGVSp from SnpEff ANN):
#    Example: Hmgcs2, c.616delC, p.Arg206fs

# 3. Find CDS start within cDNA by matching against known WT protein
#    (Ensembl protein FASTA used as ground truth)
cds_start = find_cds_start(cdna_seq, wt_protein, aa_pos)
# find_cds_start: scans ATG codons in cDNA, translates with Bio.Seq,
# checks if translation matches WT protein prefix up to aa_pos

# 4. Apply HGVSc nucleotide change
#    c.616delC → delete position 616 in CDS
mutant_cds = apply_mutation(cds_seq, start1=616, end1=616, mut_type='del', seq='C')

# 5. Translate from frameshift position using Biopython
novel_suffix = str(Seq(mutant_cds[(aa_pos-1)*3:]).translate(to_stop=False))
novel_suffix = novel_suffix[:novel_suffix.index('*')]  # truncate at stop

# 6. Extract junction peptides spanning WT|novel boundary
peptides = extract_junction_peptides(wt_protein, aa_pos, novel_suffix, [8,9,10])
```

**Results for 6 frameshift genes:**
| Gene | HGVSc | HGVSp | Novel suffix (first 20 aa) |
|------|-------|-------|---------------------------|
| Kcnq2 | c.117_118delGC | p.Leu40fs | THRGLRGPQARQRFEQAADG |
| Mtmr14 | c.33_70del | p.Ala12fs | GARALGAAGGVLPDSVPCQG |
| Ifi211 | c.293delA | p.Asn98fs | MVKKQVLQHLHQLQATC |
| Ifi211 | c.287delA | p.Asn96fs | EKMVKKQVLQHLHQLQATC |
| Gpr88 | c.592_598del | p.Ala198fs | RCLRRRRCCCTATWASCAAC |
| Hmgcs2 | c.616delC | p.Arg206fs | APQVVLGLWQC |
| Ttc14 | c.2298dupA | p.Ter767fs | MIC |
| Ttc14 | c.2142dupA | p.Ter715fs | MIC |

### 5c. Fusion Junction Peptides

**Script:** `code/steps/08b_fusion_peptides.py`
**Source:** Arriba output `peptide_sequence` column

```python
# code/steps/08b_fusion_peptides.py

# Read Arriba fusion output (all 6 samples)
# Filter to high/medium confidence fusions with predicted peptide_sequence
# Arriba uses: UPPERCASE = canonical reading frame, lowercase = out-of-frame
# '|' marks the fusion junction

# Convert to uppercase (lowercase = real amino acids in alternative frame)
pep_seq = fusion['peptide_sequence'].replace('|','').upper()

# Extract sliding-window k-mers spanning the junction
junction_pos = peptide_sequence.index('|')
for k in [8, 9, 10]:
    start = max(0, junction_pos - k)
    for i in range(start, min(junction_pos + 1, len(seq) - k + 1)):
        peptides.append(seq[i:i+k])
```

**Tools and versions:**
| Tool | Version | Purpose |
|------|---------|---------|
| pyensembl | 2.3.13 (fallback to protein FASTA) | WT protein sequence lookup |
| Biopython | 1.86 | `Bio.Seq.translate` for frameshift novel sequences |
| Ensembl GRCm38 protein FASTA | Release 102 | 22,547 gene symbols |
| Ensembl GRCm38 cDNA FASTA | Release 102 | Transcript sequences for CDS extraction |
| Arriba | 2.5.1 | Fusion junction peptide sequences |

**MHC-I binding prediction:**
```python
# code/steps/09_netmhcpan_iedb.py
# Reference:   GRCm38/mm10
# Allele:      H-2-Db (primary); H-2-Kb also queried
# Lengths:     8-mer, 9-mer, 10-mer
# Method:      netmhcpan_el (eluted ligand) via IEDB REST API
# Endpoint:    http://tools-cluster-interface.iedb.org/tools_api/mhci/
# Batch size:  100 peptides per request; exponential back-off on failure

payload = {
    "method":        "netmhcpan_el",
    "sequence_text": "\n".join(peptides),
    "allele":        "H-2-Db",   # primary allele
    "length":        "9",
}
# Parse: col 5 = peptide, col 9 = percentile_rank_EL
# Thresholds: SB < 0.5% · WB < 2.0% · PB < 10.0%
```

**Final peptide counts (H-2-Db, updated):**
| Category | Count |
|----------|-------|
| Total peptides extracted | 7,892 |
| — Missense/inframe | 7,778 |
| — Frameshift (real translated sequences) | 150 |
| — Fusion junction (Arriba 08b) | 102 |
| — Dot-fusion ORF translated (08d) | 12 |
| Peptides with H-2-Db predictions | 4,404 |
| Strong binders SB (< 0.5%) | 41 |
| Weak binders WB (< 2.0%) | 108 |
| Possible binders PB (< 10.0%) | 478 |
| **Total binders (SB+WB+PB)** | **627** |
| Frameshift binders | 1 (Hmgcs2 `AVYPSGNAA`, PB 6.0%) |
| Fusion WB binders | 4 (Fxr1::Zfp704 ×3, Zfp740::Csad ×1) |

---

## Section 6: Clonality Analysis

**Scripts:** `code/steps/05_pyclone_prep.py`, `code/steps/06_pyclone_run.sh`, `code/steps/07_pyclone_tables.R`
**Purpose:** Estimate cancer cell fraction (cellular prevalence) for each somatic mutation to identify clonal vs. subclonal neoantigens.

```bash
# Step 6: PyClone-VI (variational inference, beta-binomial model)
pyclone-vi fit \
    --in-file pyclone_vi_input.tsv \
    --out-file results.h5 \
    --num-clusters 10 \
    --density beta-binomial \
    --num-restarts 10

pyclone-vi write-results-file \
    --in-file results.h5 \
    --out-file tables/loci.tsv

# Input TSV format (per mutation, per sample):
# mutation_id  sample_id  ref_counts  alt_counts  normal_cn  major_cn  minor_cn
```

**PyClone-VI input preparation (Step 5):**
```python
# code/steps/05_pyclone_prep.py
# For each PASS somatic mutation, extract:
# - alt_counts: AD field (allelic depth)
# - ref_counts: DP - alt_counts
# - copy_number: assumed diploid (major_cn=1, minor_cn=1, normal_cn=2)
# Multi-sample: mutations absent in a sample get alt=0, ref=median_depth
```

**Tools and versions:**
| Tool | Version | Purpose |
|------|---------|---------|
| PyClone-VI | latest | Clonal inference (variational Bayes, beta-binomial) |
| Python/pandas | 3.13 / 2.x | Input preparation |
| R | 4.x | Result table processing |

**Results:**
- 3 mutations met PyClone-VI clustering criteria (sufficient depth + VAF signal)
- 2 clonal clusters identified
- The majority of somatic mutations had insufficient evidence for PyClone-VI (low VAF, borderline depth)
- Mutations without PyClone assignment were given `sample_prevalence = 0` in ranking

---

## Section 7: Neoantigen Ranking

**Script:** `code/steps/10_rank_neoantigens.py`

**Composite score formula:**
```
composite_score = 0.40 × clone_pct_norm + 0.40 × mhc_score + 0.20 × expr_pct
                × rna_penalty  [×1.0 if RNA-validated, ×0.3 if REF_ONLY in RNA-seq]

where:
  clone_pct_norm = PyClone-VI per-sample CCF normalized 0–1 (max CCF in cohort)
  mhc_score      = SB: max(0.5, 1 − rank/2); WB: 1 − (rank−0.5)/1.5 × 0.5; PB: 0.5 − (rank−2)/8 × 0.4
  expr_pct       = log2(CPM+1) / max(log2(CPM+1)) across all peptides
  rna_penalty    = 0.3 applied to genes where 0 alt reads detected in RNA-seq BAM despite coverage
```

**PyClone-VI (per-sample mode, VAF purity + cnvkit CN):** 1,330 unique mutations across 6 samples; CCF distribution:
- Clonal (CCF>0.8): **138 mutations** — identified after purity correction (VAF-based per-sample purity)
- Subclonal (0.2–0.8): 884 mutations
- Minor clone (<0.2): 308 mutations

**Purity correction impact:** Initial run with all-sample purity=1.0 → 0 clonal mutations. After per-sample VAF-based purity correction (0.18–1.0) + cnvkit copy number per locus: 138 clonal mutations identified — a qualitative improvement enabling meaningful tumor vaccine target prioritization.

**Copy number update:** cnvkit 0.9.13 run on 5 tumor BAMs vs. normal (423_D0_old) using autobin targets from WES exon BED. Per-mutation major/minor CN assigned from cnvkit .cns segments. These GemOVCA tumors are highly aneuploid (widespread LOH; 62-122 deletion segments per sample) — log2-ratio purity inference is unreliable without allele frequencies. VAF-based purity retained; cnvkit CN improves absolute CCF accuracy.

**Per-sample purity used:**
| Sample | Purity | Source |
|--------|--------|--------|
| 443_D21_new | 0.18 | VAF-based (median_VAF × 2) |
| 428_D20_new | 0.22 | VAF-based |
| 34_D52_old | 0.38 | VAF-based (no WES BAM on S3) |
| 36_D99_new | 0.22 | VAF-based |
| 38_D99_new | 0.20 | VAF-based |
| 42_D122_old | 1.00 | VAF-based (cnvkit confirms high purity) |

**Normal proteome filter:** 1,078 somatic peptides (13.9%) matched canonical GRCm38 R102 proteins exactly and were removed. Remaining 6,700 somatic + 93 fusion peptides are novel.

**RNA-seq mutant read validation (pysam pileup):**
| Gene | Peptide | depth | alt reads | VAF_RNA | Status |
|------|---------|-------|-----------|---------|--------|
| Nudt8 | QSLRPNPEEV | 103 | 60 | 58.3% | ✓ expressed |
| Tcirg1 | SHRLLLETL | 274 | 168 | 61.3% | ✓ expressed |
| Chka | YGIFPQGRL | 125 | 82 | 65.6% | ✓ expressed |
| Msr1 | SNVEMRFTI | 231 | 125 | 54.1% | ✓ expressed |
| Nisch | RSAAIPYWL | 272 | 0 | 0% | ✗ REF_ONLY (penalized) |
| Vps13b | TSIPGTPVL | 56 | 0 | 0% | ✗ REF_ONLY (penalized) |
| Plxnb2 | NLPEFIVTF | 297 | 0 | 0% | ✗ REF_ONLY (penalized) |

**Top 15 neoantigens after purity correction (H-2-Db, composite score):**
| Rank | Sample | Gene | Peptide | rank_EL% | log2CPM | Class | CCF | Score |
|------|--------|------|---------|----------|---------|-------|-----|-------|
| 1 | 36_D99_new | Birc6 | VTIEQSDEL | 0.17 | 9.76 | **SB** | 0.999 | 0.932 |
| 2 | 38_D99_new | Tent4b | SCMGNGVTL | 0.03 | 5.97 | **SB** | 1.000 | 0.895 |
| 3 | 36_D99_new | Tent4b | SCMGNGVTL | 0.03 | 5.65 | **SB** | 1.000 | 0.890 |
| 4 | 34_D52_old | Ubtd2 | GALTDCYDEL | 0.52 | 3.57 | WB | 0.997 | 0.857 |
| 5 | 428_D20_new | Slc15a2 | SLISTFITPM | 0.94 | 6.61 | WB | 0.999 | 0.853 |
| 6 | 38_D99_new | Tent4b | PSCMGNGVTL | 0.90 | 5.97 | WB | 1.000 | 0.848 |
| 7 | 38_D99_new | Tent4b | SCMGNGVTLI | 0.27 | 5.97 | **SB** | 1.000 | 0.847 |
| 10 | 34_D52_old | Ubtd2 | VALGDNQPL | 0.17 | 3.57 | **SB** | 0.997 | 0.826 |
| 12 | 36_D99_new | Slc24a2 | SAVFNILFF | 0.03 | 0.79 | **SB** | 0.999 | 0.807 |
| 13 | 428_D20_new | Slc9c1 | HQLPHTEYL | 0.03 | 0.32 | **SB** | 1.000 | 0.800 |

Note: Rankings changed dramatically after purity correction. Birc6 VTIEQSDEL (rank 1, CCF=0.999, SB 0.17%) and Tent4b SCMGNGVTL (ranks 2-3, SB 0.03%) are now the top candidates — both clonal (CCF>0.99) and strong binders in highly expressed genes. Tent4b has 3 overlapping peptides (8/9/10-mer) all in top 10.

**Expression lookup:** gene symbols mapped from Ensembl IDs via `Mus_musculus.GRCm38.102.gtf.gz`; CPM computed from featureCounts `-s 1` counts. 6,793/6,793 peptides have expression assigned (Nudt8 gap patched via direct Ensembl ID lookup).

---

## Why Results May Differ from the Previous Analysis

### 1. Variant Calling Improvements
**Read orientation bias correction:** This run uses `LearnReadOrientationModel` + `FilterMutectCalls` with the orientation model artifact file. The previous analysis likely omitted this step. OxoG artifacts (C→A changes in oxidative damage) and FFPE strand artifacts can inflate false positive calls, especially in older snap-frozen or FFPE samples. This run filters these systematically.

**Panel of Normals (PoN):** Built using two normals from this cohort (423_D0_old + 451_D0_old) with matched reference. If the previous run used a different PoN or no PoN, germline and recurrent artifact calls may differ.

**Contamination estimation:** `CalculateContamination` was run per sample. Samples with contamination > 5% would have additional false-positive filtering applied.

### 2. Allele and Database Differences
**SnpEff database:** This run uses `GRCm38.86` (Ensembl 86-based SnpEff database). If the previous run used a different Ensembl release (e.g., 75, 95, 102), functional effect classifications (synonymous vs. missense, impact scores) may differ for splice-region variants, UTR changes, and multi-allelic sites.

**Reference genome:** GRCm38/mm10 throughout. Confirmed consistent between WES and RNA-seq.

### 3. Peptide Extraction Method
**Protein FASTA version:** This run uses Ensembl release 102 protein sequences (`Mus_musculus.GRCm38.pep.all.fa.gz`). Differences in release may change the WT protein sequence used for k-mer extraction, altering which peptides span the mutation.

**Frameshift handling:** The previous analysis likely produced placeholder sequences or skipped frameshifts entirely. This run uses Biopython + Ensembl cDNA FASTA to translate actual novel sequences — identifying peptides like AAPQVVLGL (Hmgcs2) that a placeholder approach would miss.

**pyensembl compatibility:** pyensembl 2.3.13 is incompatible with polars 0.19.19 on Python 3.13, causing a Rust panic that prevents the full sliding-window approach. This run fell back to protein FASTA lookups + targeted k-mer extraction around each mutation, producing peptides only at the mutation site rather than all k-mers of the full protein.

### 4. MHC Binding Prediction
**Local vs. API:** The previous run likely used locally installed netMHCpan (binary). This run uses the IEDB REST API (`tools-cluster-interface.iedb.org`) with netMHCpan_el method. Small numerical differences in rank_EL% are expected between API versions.

**Mouse H-2 allele support:** MHCflurry 2.1.5 (installed but not used) does not support mouse H-2 alleles. The IEDB API was used for H-2-Kb and H-2-Db predictions — consistent with the prior approach if that also used IEDB.

### 5. Clonality (PyClone-VI vs. PyClone)
**Algorithm change:** PyClone-VI uses variational inference (faster) vs. the original PyClone which uses MCMC. Both use a beta-binomial model but PyClone-VI converges differently and may assign mutations to different clusters.

**Copy number input:** This run uses cnvkit 0.9.13 per-mutation copy number from WES tumor-normal analysis. 5 tumor BAMs processed (34_D52_old not on S3 → diploid fallback). Highly aneuploid tumors with 69–126 segments per sample; major/minor CN assigned per mutation locus from `.call.cns` segments.

**Multi-sample mode (fixed):** The initial run used multi-sample mode which discarded 1,327/1,330 mutations not present in all 6 samples, leaving only 3 mutations. The corrected approach uses per-sample mode independently, recovering 1,330 unique mutations with CCF assigned.

**Purity correction:** Initial purity=1.0 for all samples → 0 clonal mutations. Per-sample VAF-based purity (0.18–1.0) → 121 clonal mutations. VAF-based purity + cnvkit CN → **138 clonal mutations** (final).

### 6. Expression Data
**Strandedness error (fixed):** The previous featureCounts run used `-s 2` (reverse-stranded) on a forward-stranded library, assigning ~6% of reads. This silently produced near-zero counts and triggered fallback to STAR GeneCounts (unstranded column). The current run confirmed forward-strandedness from STAR GeneCounts column ratios and re-ran featureCounts with `-s 1`, increasing assignment by ~2–2.5× vs. the unstranded fallback.

**CPM gene symbol mapping:** featureCounts outputs Ensembl gene IDs (ENSMUSG…). Gene symbol lookup requires mapping via GTF — not done in the previous run, meaning expression was likely `NaN` for all peptides in the prior composite score. This run resolves symbols via `Mus_musculus.GRCm38.102.gtf.gz` (55,487 genes), achieving 7,865/7,892 peptides with expression data.

**Missing sample (42_D122_old):** Previously omitted from the merged matrix. Corrected in this run.

---

## Software Environment Summary

| Software | Version | Source | Purpose |
|---------|---------|--------|---------|
| GATK4 | 4.x | Broad Institute | WES preprocessing, variant calling, filtering |
| samtools | 1.x | htslib | BAM sorting, indexing |
| bcftools | 1.x | htslib | VCF manipulation |
| SnpEff | 5.x | Cingolani Lab | Functional variant annotation |
| STAR | 2.7.x | Dobin Lab | RNA-seq alignment (2-pass) |
| Arriba | 2.5.1 | Uhlmann Lab | Fusion gene detection |
| featureCounts | 2.x | Subread | Gene-level count quantification |
| PyClone-VI | latest | Roth Lab | Clonal evolution inference |
| Biopython | 1.86 | Biopython consortium | Frameshift novel sequence translation |
| IEDB API | netMHCpan_el | La Jolla Institute | MHC-I binding prediction (H-2-Db primary) |
| STAR-Fusion | 1.13.0 | Haas Lab | Complementary fusion caller (installed, planned) |
| Python | 3.13 | PSF | Peptide extraction, ranking |
| Ensembl | Release 102 | EMBL-EBI | GRCm38 protein + cDNA FASTA |

---

---

## Software — Additional Tools (Step 13+)

| Software | Version | Source | Purpose |
|---------|---------|--------|---------|
| cnvkit | 0.9.13 | Talevich Lab | WES tumor purity + copy number estimation |
| DNAcopy | Bioconductor | CBS | cnvkit segmentation (CBS algorithm) |
| Arriba | 2.5.1 | Uhlmann Lab | Fusion detection (D88_old, D99_old, D109_new re-run) |

---

## Key Files — S3 Final Results

| Path | Description |
|------|-------------|
| `final-results/ranked/neoantigens_final.tsv` | **6,793 peptides — final ranked table** (purity-corrected CCF, cnvkit CN, RNA validation) |
| `final-results/ranked/strong_binders_final.tsv` | **597 binders (SB+WB+PB, H-2-Db)** sorted by composite score |
| `final-results/ranked/neoantigens_full_annotated.tsv` | Pre-rebuild table (purity+RNA validation integrated, pre-cnvkit) |
| `final-results/peptides/all_peptides_novel.tsv` | 6,700 somatic peptides after normal proteome filter (13.9% removed) |
| `final-results/netmhcpan/all_predictions_H2Db.tsv` | H-2-Db predictions only |
| `final-results/netmhcpan/all_predictions.tsv` | H-2-Kb + H-2-Db merged |
| `final-results/peptides/fusions_all.tsv` | 08b Arriba fusion peptides (germline-filtered) |
| `final-results/peptides/fusions_dot_translated.tsv` | 08d dot-fusion ORF translations (12 peptides) |
| `final-results/rnaseq/raw_counts_s1.txt` | featureCounts output (-s 1, forward-stranded) |
| `final-results/pyclone/loci_cnvkit.tsv` | Final per-sample PyClone-VI CCF (VAF purity + cnvkit CN) |
| `final-results/pyclone/clonality_cnvkit.tsv` | Final clonality summary (138 clonal, 884 subclonal, 308 minor) |
| `final-results/pyclone/loci_per_sample_purity.tsv` | VAF-purity-only PyClone-VI CCF (intermediate) |
| `final-results/docs/pipeline_validation_report.md` | This document |

*Document updated: 2026-03-16*
*Pipeline code: `/home/ec2-user/NeoAntigen2026-aws-rerun/code/`*
*Results: `s3://neoantigen2026-rerun/results/res_20260311_225555/`*
*Final results: `s3://neoantigen2026-rerun/final-results/`*
