#!/usr/bin/env Rscript
# =============================================================================
# summarize_bam_qc.R
# Reads BAM validation TSV produced by validate_bams.sh
# Outputs: formatted summary table (TSV + HTML) and QC figures (PDF + PNG)
# Usage: Rscript summarize_bam_qc.R <results.tsv> <outdir>
# =============================================================================

# --- Dependencies ------------------------------------------------------------
local_lib <- path.expand("~/R/library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "scales",
                   "readr", "patchwork", "RColorBrewer")
optional_pkgs <- c("kableExtra")  # needs system libs; graceful fallback

for (pkg in c(required_pkgs, optional_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing: ", pkg)
    tryCatch(
      install.packages(pkg, repos = "https://cloud.r-project.org",
                       lib = local_lib, quiet = TRUE),
      error = function(e) message("  [WARN] Could not install ", pkg, ": ", e$message)
    )
  }
  if (requireNamespace(pkg, quietly = TRUE)) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  } else if (pkg %in% required_pkgs) {
    stop("Required package not available: ", pkg)
  } else {
    message("[WARN] Optional package not available: ", pkg, " — skipping HTML table")
  }
}
# kableExtra's save_kable requires pandoc; use plain HTML if unavailable
USE_KABLE <- requireNamespace("kableExtra", quietly = TRUE) &&
             isTRUE(rmarkdown::pandoc_available())

# --- Args --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
results_tsv <- ifelse(length(args) >= 1, args[1],
                      "/home/ec2-user/results/bam_validation/bam_validation_results.tsv")
outdir      <- ifelse(length(args) >= 2, args[2],
                      "/home/ec2-user/results/bam_validation")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# --- Load data ---------------------------------------------------------------
message("Reading: ", results_tsv)
df <- read_tsv(results_tsv, col_types = cols(.default = "c")) %>%
  mutate(
    size_gb       = as.numeric(size_gb),
    total_reads   = as.numeric(total_reads),
    mapped_reads  = as.numeric(mapped_reads),
    mapped_pct    = as.numeric(mapped_pct),
    duplication_pct = as.numeric(duplication_pct)
  )

# QC status colour palette
qc_colours <- c("PASS" = "#2ca25f", "WARN" = "#f0a500", "FAIL" = "#d7301f")

# --- Console summary ---------------------------------------------------------
message("\n===== BAM Validation Summary =====")
message(sprintf("Total BAM files examined : %d", nrow(df)))
message(sprintf("PASS                     : %d", sum(df$passed_qc == "PASS", na.rm = TRUE)))
message(sprintf("WARN                     : %d", sum(df$passed_qc == "WARN", na.rm = TRUE)))
message(sprintf("FAIL                     : %d", sum(df$passed_qc == "FAIL", na.rm = TRUE)))
message(sprintf("Index (.bai) present     : %d / %d",
                sum(df$index_present == "YES", na.rm = TRUE), nrow(df)))
message(sprintf("Coordinate-sorted        : %d / %d",
                sum(df$sort_order == "coordinate", na.rm = TRUE), nrow(df)))
message("===================================\n")

# --- Clean summary table (wide, formatted) -----------------------------------
summary_table <- df %>%
  select(
    Sample        = sample,
    `Size (GB)`   = size_gb,
    Quickcheck    = quickcheck,
    `Sort Order`  = sort_order,
    Reference     = reference_build,
    `Read Groups` = read_group,
    `Index (.bai)`= index_present,
    `Total Reads` = total_reads,
    `Mapped (%)`  = mapped_pct,
    `Dup (%)`     = duplication_pct,
    QC            = passed_qc,
    Notes         = notes
  ) %>%
  mutate(
    `Total Reads` = ifelse(is.na(`Total Reads`), "NA",
                           formatC(`Total Reads`, format = "d", big.mark = ",")),
    `Size (GB)`   = sprintf("%.2f", `Size (GB)`),
    `Mapped (%)`  = ifelse(is.na(`Mapped (%)`), "NA", sprintf("%.1f%%", `Mapped (%)`)),
    `Dup (%)`     = ifelse(is.na(`Dup (%)`),     "NA", sprintf("%.1f%%", `Dup (%)`))
  )

write_tsv(summary_table, file.path(outdir, "bam_qc_summary_table.tsv"))
message("Summary table saved: ", file.path(outdir, "bam_qc_summary_table.tsv"))

# --- Figure 1: QC Status Overview (bar) --------------------------------------
fig1 <- df %>%
  count(passed_qc) %>%
  ggplot(aes(x = passed_qc, y = n, fill = passed_qc)) +
  geom_col(width = 0.55, colour = "white") +
  geom_text(aes(label = n), vjust = -0.4, size = 5, fontface = "bold") +
  scale_fill_manual(values = qc_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "BAM QC Status Overview",
       x = "QC Status", y = "Number of BAM files") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"))

# --- Figure 2: File size per sample ------------------------------------------
fig2 <- df %>%
  mutate(sample = factor(sample, levels = sample[order(size_gb)])) %>%
  ggplot(aes(x = sample, y = size_gb, fill = passed_qc)) +
  geom_col(colour = "white") +
  scale_fill_manual(values = qc_colours, name = "QC") +
  coord_flip() +
  labs(title = "BAM File Sizes", x = NULL, y = "Size (GB)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# --- Figure 3: Mapping rate per sample ---------------------------------------
fig3 <- df %>%
  filter(!is.na(mapped_pct)) %>%
  mutate(sample = factor(sample, levels = sample[order(mapped_pct)])) %>%
  ggplot(aes(x = sample, y = mapped_pct, fill = passed_qc)) +
  geom_col(colour = "white") +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "#d7301f", linewidth = 0.8) +
  annotate("text", x = 0.6, y = 81, label = "80% threshold",
           colour = "#d7301f", size = 3.5, hjust = 0) +
  scale_fill_manual(values = qc_colours, name = "QC") +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(title = "Mapping Rate per BAM", x = NULL, y = "Mapped Reads (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# --- Figure 4: Duplication rate per sample -----------------------------------
fig4 <- df %>%
  filter(!is.na(duplication_pct)) %>%
  mutate(sample = factor(sample, levels = sample[order(duplication_pct)])) %>%
  ggplot(aes(x = sample, y = duplication_pct, fill = passed_qc)) +
  geom_col(colour = "white") +
  scale_fill_manual(values = qc_colours, name = "QC") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  coord_flip() +
  labs(title = "PCR Duplication Rate per BAM", x = NULL, y = "Duplicate Reads (%)") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# --- Figure 5: QC flags heatmap (checklist) ----------------------------------
flag_df <- df %>%
  transmute(
    sample,
    `Quickcheck OK`    = ifelse(quickcheck == "PASS", 1, 0),
    `Coord Sorted`     = ifelse(sort_order == "coordinate", 1, 0),
    `Has Read Groups`  = ifelse(grepl("^YES", read_group), 1, 0),
    `Index Present`    = ifelse(index_present == "YES", 1, 0),
    `Mapping >80%`     = ifelse(!is.na(mapped_pct) & mapped_pct >= 80, 1, 0)
  ) %>%
  pivot_longer(-sample, names_to = "Check", values_to = "Status")

fig5 <- flag_df %>%
  mutate(Label = ifelse(Status == 1, "PASS", "FAIL")) %>%
  ggplot(aes(x = Check, y = sample, fill = Label)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = Label), size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("PASS" = "#2ca25f", "FAIL" = "#d7301f")) +
  labs(title = "BAM QC Checklist", x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        plot.title = element_text(face = "bold"),
        legend.position = "none",
        panel.grid = element_blank())

# --- Combine and save figures ------------------------------------------------
combined <- (fig1 | fig2) / (fig3 | fig4) / fig5 +
  plot_annotation(
    title    = "BAM File Validation QC Report",
    subtitle = sprintf("S3: s3://bam-wes/NeoAntigen-aws/data/bam/  |  Generated: %s",
                       format(Sys.time(), "%Y-%m-%d %H:%M")),
    theme = theme(plot.title    = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

pdf_path <- file.path(outdir, "bam_qc_report.pdf")
png_path <- file.path(outdir, "bam_qc_report.png")

ggsave(pdf_path, combined, width = 16, height = 18, device = "pdf")
ggsave(png_path, combined, width = 16, height = 18, dpi = 150, device = "png")

message("QC report (PDF): ", pdf_path)
message("QC report (PNG): ", png_path)

# --- HTML summary table ------------------------------------------------------
html_path <- file.path(outdir, "bam_qc_summary_table.html")

if (USE_KABLE) {
  summary_table %>%
    kbl(format = "html", caption = "BAM Validation Summary") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                  full_width = FALSE) %>%
    column_spec(11,
                background = case_when(
                  summary_table$QC == "PASS" ~ "#d4edda",
                  summary_table$QC == "WARN" ~ "#fff3cd",
                  summary_table$QC == "FAIL" ~ "#f8d7da",
                  TRUE                       ~ "white"
                ),
                bold = TRUE) %>%
    save_kable(html_path)
  message("HTML table (kableExtra): ", html_path)
} else {
  # Plain HTML fallback
  qc_bg <- ifelse(summary_table$QC == "PASS", "#d4edda",
           ifelse(summary_table$QC == "WARN", "#fff3cd", "#f8d7da"))
  html_rows <- apply(summary_table, 1, function(row) {
    bg <- qc_bg[match(row["QC"], summary_table$QC)[1]]
    cells <- paste0("<td>", row, "</td>", collapse = "")
    paste0("<tr style='background:", bg, "'>", cells, "</tr>")
  })
  html_content <- paste0(
    "<!DOCTYPE html><html><head><style>",
    "table{border-collapse:collapse;font-family:sans-serif;font-size:13px}",
    "th,td{border:1px solid #ccc;padding:6px 10px}",
    "th{background:#4a4a4a;color:white}</style></head><body>",
    "<h2>BAM Validation Summary</h2><table>",
    "<tr>", paste0("<th>", colnames(summary_table), "</th>", collapse = ""), "</tr>",
    paste(html_rows, collapse = "\n"),
    "</table></body></html>"
  )
  writeLines(html_content, html_path)
  message("HTML table (plain): ", html_path)
}

message("\nDone.")
