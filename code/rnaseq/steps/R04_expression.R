#!/usr/bin/env Rscript
# =============================================================================
# R04_expression.R — DESeq2 normalization + DEA + limma batch correction
#
# Inputs:
#   - ${RNASEQ_RESULTS}/counts/raw_counts.txt  (featureCounts output or STAR merge)
#   - ${RNASEQ_META}  (samples_rnaseq.tsv — for batch/role metadata)
#
# Outputs:
#   - expression_matrix.tsv           — batch-corrected log2; rows=gene symbols,
#                                        cols=WES sample IDs (tumor + normal only)
#   - expression_matrix_full.tsv      — all RNASeq samples including extras
#   - expression_by_sample.tsv        — long format (gene, sample, expr_log2)
#   - dea_results.tsv                 — DEA: tumor-vs-D0 log2FC + padj per gene
#   - dea_results_per_sample.tsv      — per-sample DEA results (each tumor vs D0)
#   - rnaseq_qc.pdf                   — PCA, size factors, volcano plots
#
# All outputs are uploaded to ${S3_RNASEQ_OUT}/
#
# Design:
#   1. Load counts; filter low-expression genes (rowSums >= 10)
#   2. DESeq2: full model (~ batch + condition); DEA tumors vs D0 normal
#   3. VST transformation
#   4. limma::removeBatchEffect (batch = old / new)
#   5. Map Ensembl IDs to gene symbols via biomaRt (or local cache)
#   6. Write outputs and upload to S3
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
                 col_names = c("sample_id","wes_id","timepoint","batch","role","fastq_pattern"))
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

# ── DESeq2: normalization + DEA ───────────────────────────────────────────────
# condition = "normal" (D0) vs "tumor" (all other roles).
# DEA is tumor-vs-normal; positive log2FC = upregulated in tumors.
message("Running DESeq2 (normalization + DEA)...")

col_data <- data.frame(
    sample_id = meta$sample_id,
    batch     = factor(meta$batch),
    condition = factor(ifelse(meta$role == "normal", "normal", "tumor"),
                       levels = c("normal", "tumor")),
    row.names = meta$sample_id
)

dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData   = col_data,
    design    = ~ batch + condition
)

# Full DESeq2 (estimates size factors, dispersions, and GLM)
dds <- DESeq(dds, parallel = FALSE, quiet = FALSE)
norm_counts <- counts(dds, normalized = TRUE)

# ── DEA: all tumors vs D0 normal (pooled) ─────────────────────────────────────
message("Extracting DEA results (tumor vs normal)...")
dea_res <- results(dds,
    contrast       = c("condition", "tumor", "normal"),
    alpha          = 0.05,
    independentFiltering = TRUE
)
dea_df <- as.data.frame(dea_res) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    mutate(
        sig       = !is.na(padj) & padj < 0.05,
        direction = case_when(
            sig & log2FoldChange >  1  ~ "UP",
            sig & log2FoldChange < -1  ~ "DOWN",
            TRUE                       ~ "NS"
        )
    )

message(sprintf("DEA (tumor vs normal): %d UP, %d DOWN (padj<0.05, |log2FC|>1)",
                sum(dea_df$direction == "UP"),
                sum(dea_df$direction == "DOWN")))

# ── Per-timepoint DEA: each tumor sample group vs normal ─────────────────────
# Uses Wald test on the full DESeq object via LRT or individual contrasts.
# For longitudinal data we run each timepoint separately.
tumor_sids <- meta$sample_id[meta$role == "tumor"]
normal_sids_meta <- meta$sample_id[meta$role == "normal"]
message(sprintf("Per-sample DEA: %d tumor samples vs %d normal",
                length(tumor_sids), length(normal_sids_meta)))

per_sample_dea <- list()
for (sid in tumor_sids) {
    tryCatch({
        sel <- c(normal_sids_meta, sid)
        sel <- sel[sel %in% colnames(count_mat)]
        sub_counts <- count_mat[, sel, drop = FALSE]
        sub_meta   <- meta[meta$sample_id %in% sel, ]
        sub_col    <- data.frame(
            batch     = factor(sub_meta$batch),
            condition = factor(ifelse(sub_meta$role == "normal", "normal", "tumor"),
                               levels = c("normal", "tumor")),
            row.names = sub_meta$sample_id
        )
        # If only one batch level, drop it from design
        des <- if (length(unique(sub_meta$batch)) > 1) ~ batch + condition else ~ condition
        dds_s <- DESeqDataSetFromMatrix(sub_counts, sub_col, design = des)
        dds_s <- DESeq(dds_s, quiet = TRUE)
        res_s <- results(dds_s, contrast = c("condition","tumor","normal"), alpha = 0.05)
        df_s  <- as.data.frame(res_s) %>%
            tibble::rownames_to_column("gene_id") %>%
            mutate(sample_id = sid,
                   sig = !is.na(padj) & padj < 0.05,
                   direction = case_when(
                       sig & log2FoldChange >  1 ~ "UP",
                       sig & log2FoldChange < -1 ~ "DOWN",
                       TRUE ~ "NS"))
        per_sample_dea[[sid]] <- df_s
        message(sprintf("  %s: UP=%d  DOWN=%d", sid,
                        sum(df_s$direction == "UP"),
                        sum(df_s$direction == "DOWN")))
    }, error = function(e) {
        message(sprintf("  [WARN] DEA failed for %s: %s", sid, conditionMessage(e)))
    })
}
dea_per_sample_df <- bind_rows(per_sample_dea)

# VST transformation (variance stabilizing — better for PCA/batch correction)
# blind=FALSE: use the fitted model (more accurate when design is known)
message("Running VST transformation...")
vst_mat <- assay(vst(dds, blind = FALSE))

message(sprintf("VST matrix: %d genes × %d samples", nrow(vst_mat), ncol(vst_mat)))

# ── PCA before batch correction ───────────────────────────────────────────────
pca_before <- prcomp(t(vst_mat), scale. = FALSE)
pca_before_df <- as.data.frame(pca_before$x[, 1:2])
pca_before_df$sample_id <- rownames(pca_before_df)
pca_before_df <- left_join(pca_before_df, meta, by = "sample_id")
var_before <- round(100 * summary(pca_before)$importance[2, 1:2], 1)

fig_pca_before <- ggplot(pca_before_df,
       aes(x = PC1, y = PC2, colour = batch, shape = role, label = sample_id)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20) +
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
       aes(x = PC1, y = PC2, colour = batch, shape = role, label = sample_id)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 2.8, max.overlaps = 20) +
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
    if (nrow(dea_per_sample_df) > 0) {
        dea_per_sample_df <- dea_per_sample_df %>%
            mutate(gene_symbol = ifelse(gene_id %in% names(sym_lookup),
                                        sym_lookup[gene_id], gene_id))
    }
}

out_dea        <- file.path(out_dir, "dea_results.tsv")
out_dea_per    <- file.path(out_dir, "dea_results_per_sample.tsv")
write_tsv(dea_df,             out_dea)
write_tsv(dea_per_sample_df,  out_dea_per)
message("Wrote: ", out_dea)
message("Wrote: ", out_dea_per)

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

# ── QC PDF (PCA + size factors + volcano) ─────────────────────────────────────
# Volcano plot: tumor vs normal DEA
vol_df <- dea_df %>%
    mutate(
        gene_label = ifelse(direction != "NS" & abs(log2FoldChange) > 2 &
                            !is.na(padj) & padj < 0.01,
                            coalesce(gene_symbol, gene_id), NA_character_)
    )

fig_volcano <- ggplot(vol_df,
       aes(x = log2FoldChange, y = -log10(pvalue + 1e-300), colour = direction)) +
    geom_point(size = 0.8, alpha = 0.6) +
    ggrepel::geom_text_repel(aes(label = gene_label), size = 2.5,
                              na.rm = TRUE, max.overlaps = 20) +
    scale_colour_manual(values = c("UP" = "#E41A1C", "DOWN" = "#377EB8", "NS" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
    labs(title = "Volcano: Tumor vs D0 Normal (all tumors pooled)",
         x = "log2 Fold Change", y = "-log10(p-value)",
         subtitle = sprintf("UP=%d  DOWN=%d  (padj<0.05, |log2FC|>1)",
                            sum(dea_df$direction == "UP"),
                            sum(dea_df$direction == "DOWN"))) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

combined_qc <- (fig_pca_before | fig_pca_after) / (fig_sf | fig_volcano) +
    patchwork::plot_annotation(
        title    = "RNASeq QC & DEA Report",
        subtitle = sprintf("DESeq2 + VST + limma batch correction | %s",
                           format(Sys.time(), "%Y-%m-%d %H:%M")),
        theme = theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )

out_qc <- file.path(out_dir, "rnaseq_qc.pdf")
ggsave(out_qc, combined_qc, width = 16, height = 14, device = "pdf")
message("Wrote: ", out_qc)

# ── Upload to S3 ──────────────────────────────────────────────────────────────
if (nchar(s3_rnaseq) > 5) {
    message("Uploading to S3: ", s3_rnaseq)
    system(paste("aws s3 cp", shQuote(out_primary), paste0(s3_rnaseq, "/expression_matrix.tsv")))
    system(paste("aws s3 cp", shQuote(out_full),    paste0(s3_rnaseq, "/expression_matrix_full.tsv")))
    system(paste("aws s3 cp", shQuote(out_long),    paste0(s3_rnaseq, "/expression_by_sample.tsv")))
    system(paste("aws s3 cp", shQuote(out_qc),      paste0(s3_rnaseq, "/rnaseq_qc.pdf")))
    system(paste("aws s3 cp", shQuote(out_dea),     paste0(s3_rnaseq, "/dea_results.tsv")))
    system(paste("aws s3 cp", shQuote(out_dea_per), paste0(s3_rnaseq, "/dea_results_per_sample.tsv")))
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
