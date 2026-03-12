#!/usr/bin/env Rscript
# =============================================================================
# R04_expression.R — DESeq2 normalization + stage-aware DEA + limma batch correction
#
# Sample stages (defined in samples_rnaseq.tsv):
#   normal        — D0 (reference baseline)
#   early         — D20, D21
#   transitioning — D52
#   mid           — D88, D99
#   late          — D109, D122
#
# Statistical design:
#   DEA: DESeq2 model ~ batch + stage
#        Batch (old/new) is a covariate — corrects for batch in the GLM.
#        DO NOT remove batch from counts before DEA (inflates false positives).
#   Expression matrix: limma::removeBatchEffect applied to VST values ONLY
#        for visualisation and step-12 ranking. Never used as DEA input.
#
# DEA contrasts run:
#   early         vs normal  (D20/D21 vs D0)
#   transitioning vs normal  (D52 vs D0)
#   mid           vs normal  (D88/D99 vs D0)
#   late          vs normal  (D109/D122 vs D0)
#   late          vs early   (disease trajectory)
#   mid           vs early   (disease trajectory)
#
# Inputs:
#   ${RNASEQ_RESULTS}/counts/raw_counts.txt  (featureCounts or STAR merge)
#   ${RNASEQ_META}  (samples_rnaseq.tsv — 7 cols incl. stage)
#
# Outputs:
#   expression_matrix.tsv          — batch-corrected VST; WES sample IDs
#   expression_matrix_full.tsv     — all RNASeq samples
#   expression_by_sample.tsv       — long format
#   dea_results.tsv                — all contrasts combined (stage, gene, log2FC, padj)
#   rnaseq_qc.pdf                  — PCA (coloured by stage), volcanos, size factors
#
# All outputs uploaded to ${S3_RNASEQ_OUT}/
# =============================================================================

local_lib <- path.expand("~/R/library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

# ── Install/load packages ─────────────────────────────────────────────────────
bioc_pkgs <- c("DESeq2", "limma", "biomaRt", "vsn")
cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "patchwork",
               "RColorBrewer", "ggrepel", "stringr", "matrixStats")

for (pkg in cran_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org",
                         lib = local_lib, quiet = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org",
                     lib = local_lib, quiet = TRUE)

for (pkg in bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
        BiocManager::install(pkg, lib = local_lib, ask = FALSE, update = FALSE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── Paths from environment ────────────────────────────────────────────────────
results_dir  <- Sys.getenv("RESULTS_DIR",  "/home/ec2-user/results")
rnaseq_res   <- Sys.getenv("RNASEQ_RESULTS",
                            file.path(results_dir, "rnaseq"))
rnaseq_meta  <- Sys.getenv("RNASEQ_META",
                            "/home/ec2-user/NeoAntigen2026-aws-rerun/code/rnaseq/samples_rnaseq.tsv")
s3_rnaseq    <- Sys.getenv("S3_RNASEQ_OUT", "")   # set by config_rnaseq.sh
s3_counts    <- file.path(s3_rnaseq, "counts/raw_counts.txt")

counts_file  <- file.path(rnaseq_res, "counts", "raw_counts.txt")
out_dir      <- rnaseq_res
dir.create(file.path(out_dir, "counts"), recursive = TRUE, showWarnings = FALSE)

message("RNASEQ_RESULTS : ", rnaseq_res)
message("RNASEQ_META    : ", rnaseq_meta)

# ── Helper: fetch from S3 ─────────────────────────────────────────────────────
s3_fetch <- function(s3_path, local_path) {
    if (!file.exists(local_path) && nchar(s3_path) > 5) {
        message("Fetching from S3: ", s3_path)
        ret <- system(paste0("aws s3 cp \"", s3_path, "\" \"", local_path, "\""))
        if (ret != 0) message("[WARN] aws s3 cp returned non-zero for: ", s3_path)
    }
    file.exists(local_path)
}

# ── Load counts ───────────────────────────────────────────────────────────────
if (!s3_fetch(s3_counts, counts_file)) {
    stop("raw_counts.txt not found locally or on S3. Run R03_quantify.sh first.")
}

message("Loading counts: ", counts_file)
raw <- read_tsv(counts_file, comment = "#", show_col_types = FALSE)

# featureCounts has columns: Geneid Chr Start End Strand Length <bam1> <bam2> ...
# STAR merge has columns: Geneid <sample1> <sample2> ...
# Detect which format we have
gene_id_col <- "Geneid"
if (!gene_id_col %in% names(raw)) {
    # Try first column
    gene_id_col <- names(raw)[1]
}

# Drop featureCounts annotation columns (Chr, Start, End, Strand, Length)
fc_annot_cols <- c("Chr", "Start", "End", "Strand", "Length")
count_cols <- setdiff(names(raw), c(gene_id_col, fc_annot_cols))

# Extract gene IDs and count matrix
gene_ids <- raw[[gene_id_col]]
count_mat <- as.matrix(raw[, count_cols])
rownames(count_mat) <- gene_ids

# Clean up BAM path prefixes from featureCounts column names
# e.g. "/home/ec2-user/tmp/rnaseq/bams/443_D21_new.Aligned.sortedByCoord.out.bam"
# → "443_D21_new"
colnames(count_mat) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "",
                            basename(colnames(count_mat)))

message(sprintf("Raw counts: %d genes × %d samples", nrow(count_mat), ncol(count_mat)))

# ── Load sample metadata ──────────────────────────────────────────────────────
meta <- read_tsv(rnaseq_meta, comment = "#", show_col_types = FALSE,
                 col_names = c("sample_id","wes_id","timepoint","batch","role","stage","fastq_pattern"))
meta <- meta[!grepl("^\\s*#", meta$sample_id), ]  # extra safety

message(sprintf("Sample manifest: %d samples", nrow(meta)))

# Keep only samples present in count matrix
meta <- meta[meta$sample_id %in% colnames(count_mat), ]
count_mat <- count_mat[, meta$sample_id, drop = FALSE]
message(sprintf("After manifest filtering: %d genes × %d samples",
                nrow(count_mat), ncol(count_mat)))

if (ncol(count_mat) == 0) {
    stop("No samples match between counts file and manifest. Check sample IDs.")
}

# ── Filter low-count genes ────────────────────────────────────────────────────
keep <- rowSums(count_mat) >= 10
count_mat <- count_mat[keep, , drop = FALSE]
message(sprintf("After low-count filter (rowSums>=10): %d genes remain", nrow(count_mat)))

# ── DESeq2: normalization + stage-aware DEA ───────────────────────────────────
# Design: ~ batch + stage
#   batch (old/new) is a covariate — absorbs sequencing batch variance in GLM.
#   stage (normal/early/transitioning/mid/late) is the factor of interest.
#
# IMPORTANT: batch correction is done in the model here.
#   limma::removeBatchEffect is applied LATER only to the VST matrix for
#   visualisation and step-12 expression scoring. It is NEVER used as DEA input.

STAGE_LEVELS <- c("normal", "early", "transitioning", "mid", "late")
meta$stage <- factor(meta$stage, levels = STAGE_LEVELS)

message("Running DESeq2 (design: ~ batch + stage)...")

col_data <- data.frame(
    sample_id = meta$sample_id,
    batch     = factor(meta$batch),
    stage     = meta$stage,
    row.names = meta$sample_id
)

dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = col_data,
    design    = ~ batch + stage
)

dds <- DESeq(dds, parallel = FALSE, quiet = FALSE)
norm_counts <- counts(dds, normalized = TRUE)
message("DESeq2 complete.")

# ── DEA contrasts ─────────────────────────────────────────────────────────────
# All contrasts vs baseline (normal) plus trajectory comparisons.
dea_contrasts <- list(
    early_vs_normal         = c("stage", "early",         "normal"),
    transitioning_vs_normal = c("stage", "transitioning", "normal"),
    mid_vs_normal           = c("stage", "mid",           "normal"),
    late_vs_normal          = c("stage", "late",          "normal"),
    mid_vs_early            = c("stage", "mid",           "early"),
    late_vs_early           = c("stage", "late",          "early")
)

dea_all <- lapply(names(dea_contrasts), function(cname) {
    contrast <- dea_contrasts[[cname]]
    # Check both levels are present in the data
    if (!all(contrast[2:3] %in% levels(dds$stage))) {
        message(sprintf("  [SKIP] %s — stage level not in data", cname))
        return(NULL)
    }
    tryCatch({
        res <- results(dds,
            contrast             = contrast,
            alpha                = 0.05,
            independentFiltering = TRUE
        )
        df <- as.data.frame(res) %>%
            tibble::rownames_to_column("gene_id") %>%
            mutate(
                contrast  = cname,
                numerator = contrast[2],
                ref_stage = contrast[3],
                sig       = !is.na(padj) & padj < 0.05,
                direction = case_when(
                    sig & log2FoldChange >  1 ~ "UP",
                    sig & log2FoldChange < -1 ~ "DOWN",
                    TRUE                      ~ "NS"
                )
            )
        message(sprintf("  DEA %-30s  UP=%d  DOWN=%d",
                        cname,
                        sum(df$direction == "UP"),
                        sum(df$direction == "DOWN")))
        df
    }, error = function(e) {
        message(sprintf("  [WARN] DEA failed for %s: %s", cname, conditionMessage(e)))
        NULL
    })
})

dea_df <- bind_rows(Filter(Negate(is.null), dea_all)) %>%
    arrange(contrast, padj, desc(abs(log2FoldChange)))
message(sprintf("DEA complete: %d rows across %d contrasts",
                nrow(dea_df), n_distinct(dea_df$contrast)))

# ── VST transformation for visualisation (NOT for DEA input) ──────────────────
# blind=FALSE: use fitted model dispersions (recommended when design is known).
message("Running VST transformation (for visualisation only)...")
vst_mat <- assay(vst(dds, blind = FALSE))

message(sprintf("VST matrix: %d genes × %d samples", nrow(vst_mat), ncol(vst_mat)))

# ── PCA before batch correction ───────────────────────────────────────────────
pca_before <- prcomp(t(vst_mat), scale. = FALSE)
pca_before_df <- as.data.frame(pca_before$x[, 1:2])
pca_before_df$sample_id <- rownames(pca_before_df)
pca_before_df <- left_join(pca_before_df, meta, by = "sample_id")
var_before <- round(100 * summary(pca_before)$importance[2, 1:2], 1)

stage_colours <- c("normal" = "#4DAF4A", "early" = "#FF7F00",
                   "transitioning" = "#984EA3", "mid" = "#E41A1C",
                   "late" = "#377EB8")

fig_pca_before <- ggplot(pca_before_df,
       aes(x = PC1, y = PC2, colour = stage, shape = batch, label = sample_id)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20) +
    scale_colour_manual(values = stage_colours) +
    labs(title = "PCA before batch correction (VST)",
         x = sprintf("PC1 (%.1f%%)", var_before[1]),
         y = sprintf("PC2 (%.1f%%)", var_before[2])) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

# ── limma batch correction ────────────────────────────────────────────────────
# Matches the approach used in step 11 (11_rnaseq_expression.R)
message("Applying limma::removeBatchEffect (batch = old/new)...")
batch_vec <- meta$batch[match(colnames(vst_mat), meta$sample_id)]

corrected_mat <- limma::removeBatchEffect(vst_mat, batch = batch_vec)
message("Batch correction complete.")

# ── PCA after batch correction ────────────────────────────────────────────────
pca_after <- prcomp(t(corrected_mat), scale. = FALSE)
pca_after_df <- as.data.frame(pca_after$x[, 1:2])
pca_after_df$sample_id <- rownames(pca_after_df)
pca_after_df <- left_join(pca_after_df, meta, by = "sample_id")
var_after <- round(100 * summary(pca_after)$importance[2, 1:2], 1)

fig_pca_after <- ggplot(pca_after_df,
       aes(x = PC1, y = PC2, colour = stage, shape = batch, label = sample_id)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20) +
    scale_colour_manual(values = stage_colours) +
    labs(title = "PCA after batch correction (VST + limma)",
         x = sprintf("PC1 (%.1f%%)", var_after[1]),
         y = sprintf("PC2 (%.1f%%)", var_after[2])) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

# ── Size factor plot ──────────────────────────────────────────────────────────
sf_df <- data.frame(
    sample_id    = names(sizeFactors(dds)),
    size_factor  = sizeFactors(dds)
)
sf_df <- left_join(sf_df, meta, by = "sample_id")

fig_sf <- ggplot(sf_df, aes(x = reorder(sample_id, size_factor),
                              y = size_factor, fill = batch)) +
    geom_col(colour = "white") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c("old" = "#E41A1C", "new" = "#377EB8")) +
    labs(title = "DESeq2 Size Factors",
         x = "Sample", y = "Size Factor") +
    theme_classic(base_size = 11) +
    coord_flip() +
    theme(plot.title = element_text(face = "bold"))

# ── Map Ensembl gene IDs to symbols ──────────────────────────────────────────
message("Mapping Ensembl IDs to gene symbols...")
gene_ids_to_map <- rownames(corrected_mat)

# Check if IDs look like Ensembl IDs
is_ensembl <- grepl("^ENSMUSG", gene_ids_to_map[1])

gene_map <- NULL
if (is_ensembl) {
    # Try biomaRt first
    tryCatch({
        mart <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                                  host = "https://nov2020.archive.ensembl.org")  # GRCm38.86
        bm <- biomaRt::getBM(
            attributes = c("ensembl_gene_id", "mgi_symbol"),
            filters    = "ensembl_gene_id",
            values     = gene_ids_to_map,
            mart       = mart
        )
        gene_map <- bm[bm$mgi_symbol != "", ]
        gene_map <- gene_map[!duplicated(gene_map$ensembl_gene_id), ]
        message(sprintf("biomaRt: mapped %d / %d Ensembl IDs to gene symbols",
                        nrow(gene_map), length(gene_ids_to_map)))
    }, error = function(e) {
        message("[WARN] biomaRt query failed: ", conditionMessage(e))
        message("[WARN] Using Ensembl IDs as gene symbols.")
    })
}

# Build a gene_id → symbol vector (fall back to gene_id if no symbol)
if (!is.null(gene_map) && nrow(gene_map) > 0) {
    sym_vec <- setNames(gene_map$mgi_symbol, gene_map$ensembl_gene_id)
    gene_symbols <- ifelse(gene_ids_to_map %in% names(sym_vec),
                           sym_vec[gene_ids_to_map],
                           gene_ids_to_map)
} else {
    # Strip version suffix if present (e.g. ENSMUSG00000000001.3 → ENSMUSG00000000001)
    gene_symbols <- sub("\\.[0-9]+$", "", gene_ids_to_map)
}

# Handle duplicate symbols (keep highest-variance)
if (anyDuplicated(gene_symbols)) {
    message("Resolving duplicate gene symbols (keeping highest variance)...")
    row_vars <- matrixStats::rowVars(corrected_mat)
    dup_df <- data.frame(
        symbol   = gene_symbols,
        row_var  = row_vars,
        orig_idx = seq_along(gene_symbols)
    )
    dup_df <- dup_df[order(dup_df$symbol, -dup_df$row_var), ]
    keep_idx <- !duplicated(dup_df$symbol)
    keep_rows <- sort(dup_df$orig_idx[keep_idx])
    corrected_mat <- corrected_mat[keep_rows, , drop = FALSE]
    gene_symbols   <- gene_symbols[keep_rows]
    message(sprintf("After deduplication: %d unique gene symbols", nrow(corrected_mat)))
}

rownames(corrected_mat) <- gene_symbols

# ── Build output matrices ─────────────────────────────────────────────────────
# Determine which samples to include in the "primary" output
# (WES-linked tumor samples + normal; matching step 10 expectations)
wes_linked <- meta$sample_id[meta$wes_id != "NA"]
normal_sids <- meta$sample_id[meta$role == "normal"]
primary_sids <- union(wes_linked, normal_sids)
primary_sids <- intersect(primary_sids, colnames(corrected_mat))

message(sprintf("Primary output samples (%d): %s",
                length(primary_sids), paste(primary_sids, collapse = ", ")))

# Map column names: use wes_id where available, else sample_id
sid_to_wes <- setNames(meta$wes_id, meta$sample_id)
primary_wes_cols <- ifelse(sid_to_wes[primary_sids] != "NA",
                            sid_to_wes[primary_sids],
                            primary_sids)

expr_primary <- corrected_mat[, primary_sids, drop = FALSE]
colnames(expr_primary) <- primary_wes_cols

expr_full <- corrected_mat  # all RNASeq samples

# ── Write DEA outputs ─────────────────────────────────────────────────────────
# Add gene symbols to DEA results (reuse the gene_symbols vector built above)
if (exists("gene_symbols") && length(gene_symbols) == nrow(corrected_mat)) {
    sym_lookup <- setNames(gene_symbols, rownames(corrected_mat))
    dea_df <- dea_df %>%
        mutate(gene_symbol = ifelse(gene_id %in% names(sym_lookup),
                                    sym_lookup[gene_id], gene_id))
}

out_dea <- file.path(out_dir, "dea_results.tsv")
write_tsv(dea_df, out_dea)
message("Wrote: ", out_dea)

# ── Write TSV outputs ─────────────────────────────────────────────────────────
# expression_matrix.tsv (primary — consumed by step 12)
expr_primary_df <- as.data.frame(expr_primary)
expr_primary_df <- tibble::rownames_to_column(expr_primary_df, var = "gene_id")
out_primary <- file.path(out_dir, "expression_matrix.tsv")
write_tsv(expr_primary_df, out_primary)
message("Wrote: ", out_primary)

# expression_matrix_full.tsv (all samples)
expr_full_df <- as.data.frame(expr_full)
expr_full_df <- tibble::rownames_to_column(expr_full_df, var = "gene_id")
out_full <- file.path(out_dir, "expression_matrix_full.tsv")
write_tsv(expr_full_df, out_full)
message("Wrote: ", out_full)

# expression_by_sample.tsv (long format)
expr_long <- expr_primary_df %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample_id", values_to = "expr_log2")
out_long <- file.path(out_dir, "expression_by_sample.tsv")
write_tsv(expr_long, out_long)
message("Wrote: ", out_long)

# ── QC PDF (PCA + size factors + per-contrast volcanos) ───────────────────────

# Faceted volcano plot: one panel per contrast
vol_df <- dea_df %>%
    mutate(
        gene_label = ifelse(direction != "NS" & abs(log2FoldChange) > 2 &
                            !is.na(padj) & padj < 0.01,
                            coalesce(gene_symbol, gene_id), NA_character_),
        contrast = factor(contrast, levels = names(dea_contrasts))
    )

n_up   <- sum(dea_df$direction == "UP",   na.rm = TRUE)
n_down <- sum(dea_df$direction == "DOWN", na.rm = TRUE)

fig_volcano <- ggplot(vol_df,
       aes(x = log2FoldChange, y = -log10(pvalue + 1e-300), colour = direction)) +
    geom_point(size = 0.7, alpha = 0.5) +
    ggrepel::geom_text_repel(aes(label = gene_label), size = 2.2,
                              na.rm = TRUE, max.overlaps = 10) +
    scale_colour_manual(values = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NS" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    facet_wrap(~ contrast, nrow = 2, scales = "free_y") +
    labs(title = "Volcano plots per DEA contrast",
         x = "log2 Fold Change", y = "-log10(p-value)",
         subtitle = sprintf("Total across all contrasts: UP=%d  DOWN=%d  (padj<0.05, |log2FC|>1)",
                            n_up, n_down)) +
    theme_classic(base_size = 9) +
    theme(plot.title   = element_text(face = "bold"),
          strip.text   = element_text(size = 7, face = "bold"),
          legend.position = "bottom")

# Page 1: PCA + size factors
page1 <- (fig_pca_before | fig_pca_after) / (fig_sf | patchwork::plot_spacer()) +
    patchwork::plot_annotation(
        title    = "RNASeq QC Report — PCA & Size Factors",
        subtitle = sprintf("DESeq2 + VST + limma batch correction | %s",
                           format(Sys.time(), "%Y-%m-%d %H:%M")),
        theme = theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )

# Page 2: per-contrast volcano plots
page2 <- fig_volcano +
    patchwork::plot_annotation(
        title    = "RNASeq DEA — Stage Contrasts",
        subtitle = sprintf("DESeq2 design: ~ batch + stage | %s",
                           format(Sys.time(), "%Y-%m-%d %H:%M")),
        theme = theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )

out_qc <- file.path(out_dir, "rnaseq_qc.pdf")
pdf(out_qc, width = 16, height = 12)
print(page1)
print(page2)
dev.off()
message("Wrote: ", out_qc)

# ── Upload to S3 ──────────────────────────────────────────────────────────────
if (nchar(s3_rnaseq) > 5) {
    message("Uploading to S3: ", s3_rnaseq)
    system(paste("aws s3 cp", shQuote(out_primary), paste0(s3_rnaseq, "/expression_matrix.tsv")))
    system(paste("aws s3 cp", shQuote(out_full),    paste0(s3_rnaseq, "/expression_matrix_full.tsv")))
    system(paste("aws s3 cp", shQuote(out_long),    paste0(s3_rnaseq, "/expression_by_sample.tsv")))
    system(paste("aws s3 cp", shQuote(out_qc),  paste0(s3_rnaseq, "/rnaseq_qc.pdf")))
    system(paste("aws s3 cp", shQuote(out_dea), paste0(s3_rnaseq, "/dea_results.tsv")))
    message("S3 upload complete.")
} else {
    message("[WARN] S3_RNASEQ_OUT not set — skipping S3 upload")
}

# ── Summary ───────────────────────────────────────────────────────────────────
message("══════════════════════════════════════════════════════")
message(sprintf("  Genes in expression matrix : %d", nrow(expr_primary_df)))
message(sprintf("  Primary output samples     : %d (%s)",
                length(primary_wes_cols),
                paste(primary_wes_cols, collapse = ", ")))
message(sprintf("  Full matrix samples        : %d", ncol(expr_full)))
message(sprintf("  Output directory           : %s", out_dir))
message("══════════════════════════════════════════════════════")
message("R04_expression.R complete.")
