#!/usr/bin/env Rscript
# =============================================================================
# summarize_variants.R
# Reads variant calling TSV from variant_calling.sh
# Outputs: summary table (TSV + HTML) and QC figures (PDF + PNG)
# Usage: Rscript summarize_variants.R <results.tsv> <outdir>
# =============================================================================

local_lib <- path.expand("~/R/library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "readr", "patchwork",
                   "RColorBrewer", "scales")
for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org",
                         lib = local_lib, quiet = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# --- Args --------------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
tsv     <- ifelse(length(args) >= 1, args[1],
                  "/home/ec2-user/results/variant_calling/variant_calling_results.tsv")
outdir  <- ifelse(length(args) >= 2, args[2],
                  "/home/ec2-user/results/variant_calling")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- Load data ---------------------------------------------------------------
message("Reading: ", tsv)
df <- read_tsv(tsv, col_types = cols(.default = "c")) %>%
    mutate(
        total_variants = as.numeric(total_variants),
        pass_variants  = as.numeric(pass_variants),
        snv_pass       = as.numeric(snv_pass),
        indel_pass     = as.numeric(indel_pass),
        # Extract timepoint from sample name
        timepoint = as.numeric(gsub(".*_(D)(\\d+)_.*", "\\2", sample)),
        cohort    = ifelse(grepl("new", sample), "new", "old")
    )

message("\n===== Variant Calling Summary =====")
print(df %>% select(sample, pass_variants, snv_pass, indel_pass, status))
message("====================================\n")

# Save summary table
write_tsv(df, file.path(outdir, "variant_summary_table.tsv"))

# Colour palette
status_col <- c("PASS" = "#2ca25f", "FAIL_MUTECT2" = "#d7301f",
                "FAIL_FILTER" = "#f0a500")
cohort_col <- c("old" = "#4393c3", "new" = "#d6604d")

# --- Figure 1: Total vs PASS variants per sample ----------------------------
fig1 <- df %>%
    filter(!is.na(pass_variants)) %>%
    select(sample, total_variants, pass_variants) %>%
    pivot_longer(-sample, names_to = "type", values_to = "count") %>%
    mutate(type = recode(type,
                         total_variants = "Total (raw)",
                         pass_variants  = "PASS (filtered)"),
           sample = factor(sample, levels = df$sample[order(df$timepoint)])) %>%
    ggplot(aes(x = sample, y = count, fill = type)) +
    geom_col(position = "dodge", colour = "white") +
    scale_fill_manual(values = c("Total (raw)" = "#aec7e8",
                                 "PASS (filtered)" = "#1f77b4"),
                      name = NULL) +
    scale_y_continuous(labels = comma) +
    coord_flip() +
    labs(title = "Somatic Variants per Sample",
         subtitle = "Raw Mutect2 calls vs PASS after FilterMutectCalls",
         x = NULL, y = "Number of Variants") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

# --- Figure 2: SNV vs Indel breakdown (PASS only) ---------------------------
fig2 <- df %>%
    filter(!is.na(snv_pass)) %>%
    select(sample, snv_pass, indel_pass) %>%
    pivot_longer(-sample, names_to = "var_type", values_to = "count") %>%
    mutate(var_type = recode(var_type, snv_pass = "SNV", indel_pass = "Indel"),
           sample   = factor(sample, levels = df$sample[order(df$timepoint)])) %>%
    ggplot(aes(x = sample, y = count, fill = var_type)) +
    geom_col(colour = "white") +
    scale_fill_manual(values = c("SNV" = "#e6550d", "Indel" = "#3182bd"),
                      name = "Variant Type") +
    scale_y_continuous(labels = comma) +
    coord_flip() +
    labs(title = "SNV vs Indel Breakdown (PASS)",
         x = NULL, y = "Number of PASS Variants") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "top")

# --- Figure 3: Tumour mutational burden over time ----------------------------
fig3 <- df %>%
    filter(!is.na(timepoint), !is.na(pass_variants)) %>%
    ggplot(aes(x = timepoint, y = pass_variants, colour = cohort,
               label = sample)) +
    geom_line(aes(group = cohort), linewidth = 0.8, linetype = "dashed",
              alpha = 0.6) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    scale_colour_manual(values = cohort_col, name = "Cohort") +
    scale_y_continuous(labels = comma) +
    scale_x_continuous(breaks = c(0, 20, 52, 99, 122),
                       labels = paste0("D", c(0, 20, 52, 99, 122))) +
    labs(title = "Somatic Mutational Burden Over Time",
         subtitle = "PASS variants per timepoint",
         x = "Timepoint", y = "PASS Variants") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

# --- Figure 4: Filter pass rate ---------------------------------------------
fig4 <- df %>%
    filter(!is.na(total_variants), total_variants > 0) %>%
    mutate(pass_rate = 100 * pass_variants / total_variants,
           sample = factor(sample, levels = sample[order(timepoint)])) %>%
    ggplot(aes(x = sample, y = pass_rate, fill = cohort)) +
    geom_col(colour = "white") +
    scale_fill_manual(values = cohort_col, name = "Cohort") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 105)) +
    coord_flip() +
    labs(title = "Variant Filter Pass Rate",
         subtitle = "PASS / Total Mutect2 calls",
         x = NULL, y = "Pass Rate (%)") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

# Install ggrepel if needed
if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org",
                     lib = local_lib, quiet = TRUE)
if (requireNamespace("ggrepel", quietly = TRUE))
    library(ggrepel)

# --- Combine figures ---------------------------------------------------------
combined <- (fig1 | fig2) / (fig3 | fig4) +
    plot_annotation(
        title    = "Somatic Variant Calling QC Report — GATK Mutect2",
        subtitle = sprintf("Reference: mm10/GRCm38  |  Normal: 423_D0_old  |  Generated: %s",
                           format(Sys.time(), "%Y-%m-%d %H:%M")),
        theme = theme(plot.title    = element_text(size = 15, face = "bold"),
                      plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )

pdf_path <- file.path(outdir, "variant_calling_report.pdf")
png_path <- file.path(outdir, "variant_calling_report.png")
ggsave(pdf_path, combined, width = 16, height = 14, device = "pdf")
ggsave(png_path, combined, width = 16, height = 14, dpi = 150, device = "png")
message("PDF: ", pdf_path)
message("PNG: ", png_path)

# --- Plain HTML summary table ------------------------------------------------
html_path <- file.path(outdir, "variant_summary_table.html")
tbl <- df %>%
    transmute(
        Sample        = sample,
        Normal        = normal,
        `Total Calls` = formatC(total_variants, format = "d", big.mark = ","),
        `PASS`        = formatC(pass_variants,  format = "d", big.mark = ","),
        `SNV (PASS)`  = formatC(snv_pass,       format = "d", big.mark = ","),
        `Indel (PASS)`= formatC(indel_pass,      format = "d", big.mark = ","),
        `Pass Rate`   = ifelse(!is.na(total_variants) & total_variants > 0,
                               sprintf("%.1f%%", 100*pass_variants/total_variants), "NA"),
        Status        = status
    )
bg <- ifelse(tbl$Status == "PASS", "#d4edda", "#f8d7da")
rows <- apply(tbl, 1, function(r) {
    paste0("<tr style='background:", bg[match(r["Status"], tbl$Status)[1]], "'>",
           paste0("<td>", r, "</td>", collapse = ""), "</tr>")
})
html <- paste0(
    "<!DOCTYPE html><html><head><style>",
    "table{border-collapse:collapse;font-family:sans-serif;font-size:13px}",
    "th,td{border:1px solid #ccc;padding:6px 12px}",
    "th{background:#4a4a4a;color:white}</style></head><body>",
    "<h2>Somatic Variant Calling Summary</h2>",
    "<p>Tool: GATK4 Mutect2 | Reference: mm10/GRCm38 | Normal: 423_D0_old (D0_old)</p>",
    "<table><tr>", paste0("<th>", colnames(tbl), "</th>", collapse = ""), "</tr>",
    paste(rows, collapse = "\n"), "</table></body></html>"
)
writeLines(html, html_path)
message("HTML: ", html_path)
message("\nDone.")
