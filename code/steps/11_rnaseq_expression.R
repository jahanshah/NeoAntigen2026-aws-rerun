#!/usr/bin/env Rscript
# =============================================================================
# Step 11 — RNASeq Expression Integration
# Downloads normalised RNA-seq counts, removes old/new batch effects,
# and produces a per-gene per-sample expression matrix used for
# neoantigen prioritisation in Step 10.
#
# Batch correction: limma::removeBatchEffect (old vs new sequencing run)
# Outputs:
#   results/rnaseq/expression_matrix.tsv      — batch-corrected log2 expression
#   results/rnaseq/expression_by_sample.tsv   — long format: gene × sample × expr
#   results/rnaseq/rnaseq_report.pdf/png      — QC figures (PCA before/after)
# =============================================================================

local_lib <- path.expand("~/R/library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

for (pkg in c("ggplot2","dplyr","tidyr","readr","patchwork",
              "RColorBrewer","scales","ggrepel","limma")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org",
                         lib = local_lib, quiet = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
OUT_DIR  <- "/home/ec2-user/results/rnaseq"
TMP_DIR  <- "/tmp/rnaseq_counts"
S3_COUNTS <- "s3://bam-wes/NeoAntigen-aws/data/RNASeq/count_files"
S3_OUT   <- "s3://bam-wes/NeoAntigen-aws/results/rnaseq"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Download counts ───────────────────────────────────────────────────────────
norm_local <- file.path(TMP_DIR, "normalizedCounts.txt")
raw_local  <- file.path(TMP_DIR, "rawCounts.txt")

if (!file.exists(norm_local)) {
    system(paste0("aws s3 cp \"", S3_COUNTS,
                  "/Longitudinal GEM-Mouse Project (RNA)_normalizedCounts.txt\" \"",
                  norm_local, "\""))
}
if (!file.exists(raw_local)) {
    system(paste0("aws s3 cp \"", S3_COUNTS,
                  "/Longitudinal GEM-Mouse Project (RNA)_rawCounts.txt\" \"",
                  raw_local, "\""))
}

if (!file.exists(norm_local)) stop("Normalised counts not found.")

counts_norm <- read_tsv(norm_local, show_col_types = FALSE)
message(sprintf("Normalised counts: %d genes × %d samples",
                nrow(counts_norm), ncol(counts_norm) - 1))

# ── Parse sample metadata ─────────────────────────────────────────────────────
# Column name pattern: "DAY SAMPLE_ID_Day_NN_Old/New"
# Assign: sample_id (matching WES), timepoint (numeric), batch (old/new)
col_map <- data.frame(
    col_name = names(counts_norm)[-1],
    stringsAsFactors = FALSE
)

col_map <- col_map %>%
    mutate(
        sample_id = case_when(
            grepl("S428|Day_20_New",    col_name, ignore.case = TRUE) ~ "428_D20_new",
            grepl("S443|Day_21_New",    col_name, ignore.case = TRUE) ~ "443_D21_new",
            grepl("S34|Day_52_Old",     col_name, ignore.case = TRUE) ~ "34_D52_old",
            grepl("S80|Day_88_Old",     col_name, ignore.case = TRUE) ~ "D88_old",
            grepl("S32.*Day_99.*Old|Day_99.*Old.*1", col_name, ignore.case = TRUE) ~ "D99_old",
            grepl("S36|Day_99.*New.*2", col_name, ignore.case = TRUE) ~ "36_D99_new",
            grepl("S38|Day_99.*New.*3", col_name, ignore.case = TRUE) ~ "38_D99_new",
            grepl("S2661|Day_109_New",  col_name, ignore.case = TRUE) ~ "D109_new",
            grepl("S42|Day_122_Old",    col_name, ignore.case = TRUE) ~ "42_D122_old",
            grepl("Day_0_Old",          col_name, ignore.case = TRUE) ~ "423_D0_old",
            TRUE ~ col_name
        ),
        timepoint = as.integer(
            gsub(".*?(\\d+).*", "\\1",
                 gsub("Day_", "", regmatches(col_name,
                        regexpr("Day_\\d+", col_name))))
        ),
        batch = ifelse(grepl("old", col_name, ignore.case = TRUE), "old", "new")
    )

message("Sample map:")
print(col_map)

# ── Build expression matrix ───────────────────────────────────────────────────
gene_ids  <- counts_norm$Symbol
expr_mat  <- as.matrix(counts_norm[, -1])
rownames(expr_mat) <- gene_ids
colnames(expr_mat) <- col_map$sample_id

# Log2 transform (normalised values are already normalised counts)
# Add small pseudocount
expr_log2 <- log2(expr_mat + 0.5)

# ── Batch effect removal (limma) ──────────────────────────────────────────────
batch <- col_map$batch
message(sprintf("Removing batch effects (old=%d  new=%d samples)",
                sum(batch == "old"), sum(batch == "new")))

expr_corrected <- limma::removeBatchEffect(expr_log2, batch = batch)

# ── PCA before / after batch correction ──────────────────────────────────────
pca_plot <- function(mat, meta, title) {
    pca <- prcomp(t(mat), scale. = FALSE)
    pct  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    df   <- data.frame(
        PC1 = pca$x[,1], PC2 = pca$x[,2],
        sample_id = meta$sample_id,
        batch     = meta$batch,
        timepoint = meta$timepoint
    )
    ggplot(df, aes(PC1, PC2, colour = batch, label = sample_id)) +
        geom_point(size = 3, aes(shape = batch)) +
        geom_text_repel(size = 3, max.overlaps = 20) +
        scale_colour_manual(values = c("old" = "#E41A1C", "new" = "#377EB8")) +
        labs(title = title,
             x = sprintf("PC1 (%.1f%%)", pct[1]),
             y = sprintf("PC2 (%.1f%%)", pct[2])) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold"))
}

fig1 <- pca_plot(expr_log2,    col_map, "PCA — Before Batch Correction")
fig2 <- pca_plot(expr_corrected, col_map, "PCA — After Batch Correction (limma)")

# Density of expression before / after
df_dens <- rbind(
    data.frame(expr = as.vector(expr_log2),     type = "Before"),
    data.frame(expr = as.vector(expr_corrected), type = "After")
)
fig3 <- ggplot(df_dens, aes(x = expr, colour = type)) +
    geom_density(linewidth = 0.9) +
    scale_colour_manual(values = c("Before" = "#E41A1C", "After" = "#377EB8")) +
    labs(title = "Expression Density Before/After Batch Correction",
         x = "log2(normalised count + 0.5)", y = "Density", colour = "") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

# Top variable genes heatmap (top 50)
var_genes  <- order(apply(expr_corrected, 1, var), decreasing = TRUE)[1:50]
heatmap_df <- expr_corrected[var_genes, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "expr") %>%
    left_join(col_map %>% select(sample_id, batch, timepoint), by = "sample_id") %>%
    mutate(sample_label = paste0(sample_id, "\n(", batch, ")"))

fig4 <- ggplot(heatmap_df,
               aes(x = reorder(sample_label, timepoint),
                   y = reorder(gene, expr),
                   fill = expr)) +
    geom_tile() +
    scale_fill_gradient2(low = "#313695", mid = "white", high = "#a50026",
                         midpoint = median(heatmap_df$expr),
                         name = "log2\nexpr") +
    labs(title = "Top 50 Variable Genes (after batch correction)",
         x = "Sample", y = "Gene") +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
          axis.text.y = element_text(size = 6),
          plot.title  = element_text(face = "bold"))

combined <- (fig1 | fig2) / (fig3 | fig4) +
    plot_annotation(
        title    = "RNASeq Expression Analysis",
        subtitle = sprintf("Batch correction: limma::removeBatchEffect  |  %d genes × %d samples  |  %s",
                           nrow(expr_corrected), ncol(expr_corrected),
                           format(Sys.time(), "%Y-%m-%d %H:%M")),
        theme = theme(plot.title    = element_text(size = 14, face = "bold"),
                      plot.subtitle = element_text(size = 8,  colour = "grey40"))
    )

ggsave(file.path(OUT_DIR, "rnaseq_report.pdf"), combined, width = 16, height = 14)
ggsave(file.path(OUT_DIR, "rnaseq_report.png"), combined, width = 16, height = 14, dpi = 150)
message("Saved: rnaseq_report.pdf/png")

# ── Write outputs ─────────────────────────────────────────────────────────────
# Wide matrix
expr_wide <- as.data.frame(expr_corrected)
expr_wide <- tibble::rownames_to_column(expr_wide, "gene_id")
write_tsv(expr_wide, file.path(OUT_DIR, "expression_matrix.tsv"))
message(sprintf("Saved: expression_matrix.tsv (%d genes)", nrow(expr_wide)))

# Long format for joining to neoantigen candidates
expr_long <- expr_wide %>%
    pivot_longer(-gene_id, names_to = "sample_id", values_to = "expr_log2") %>%
    left_join(col_map %>% select(sample_id, batch, timepoint), by = "sample_id")
write_tsv(expr_long, file.path(OUT_DIR, "expression_by_sample.tsv"))
message(sprintf("Saved: expression_by_sample.tsv (%d rows)", nrow(expr_long)))

# Upload
system(paste("aws s3 sync", OUT_DIR, S3_OUT))
message("Uploaded to S3: ", S3_OUT)
message("Step 11 complete.")
