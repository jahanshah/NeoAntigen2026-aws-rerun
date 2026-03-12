#!/usr/bin/env Rscript
# =============================================================================
# Step 10 — Neoantigen Candidate Filtering, Ranking & Final Report
#
# Merges:
#   - netMHCpan predictions  (WES missense/frameshift + fusion junction peptides)
#   - PyClone clonal data    (cluster assignments, cellular prevalence)
#   - RNASeq expression      (batch-corrected log2 expression, Step 11)
#   - Peptide metadata       (mutation type: missense, frameshift, fusion)
#
# Priority score (0–6):
#   Binding      : SB=3  WB=1
#   Mutation type: missense=2  frameshift/fusion=1  (novel epitopes preferred)
#   Clonality    : Clonal=1
#   Expression   : expressed (expr_log2 > 1 in that sample) = 1  (adj. tumour)
#
# Outputs:
#   results/final/neoantigen_candidates.tsv   — all SB+WB with scores
#   results/final/neoantigen_top_ranked.tsv   — SB only, ranked, de-duplicated
#   results/final/binding_summary.tsv         — counts by allele/length/class
#   results/final/neoantigen_report.pdf/png   — 6-panel QC/summary figure
# =============================================================================

local_lib <- path.expand("~/R/library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

for (pkg in c("ggplot2","dplyr","tidyr","readr","patchwork",
              "RColorBrewer","scales","ggrepel","stringr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos = "https://cloud.r-project.org",
                         lib = local_lib, quiet = TRUE)
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ── Paths ─────────────────────────────────────────────────────────────────────
# Use RESULTS_DIR env var (set by config.sh) if available; fall back to default.
BASE_DIR    <- Sys.getenv("RESULTS_DIR",
                           file.path(Sys.getenv("HOME", "/home/ec2-user"),
                                     "results",
                                     Sys.getenv("RUN_ID", "")))
if (!nzchar(BASE_DIR) || !file.exists(dirname(BASE_DIR))) {
    BASE_DIR <- "/home/ec2-user/results"
}
MHC_DIR     <- file.path(BASE_DIR, "netmhcpan")
PEP_DIR     <- file.path(BASE_DIR, "peptides")
PYCLONE_DIR <- file.path(BASE_DIR, "pyclone/tables")
RNASEQ_DIR  <- file.path(BASE_DIR, "rnaseq")
OUT_DIR     <- file.path(BASE_DIR, "final")

# S3 paths — use S3_RESULTS env var (set by config.sh) if available;
# otherwise construct from RUN_ID or fall back to new bucket root.
run_id_env  <- Sys.getenv("RUN_ID", "")
s3_results_env <- Sys.getenv("S3_RESULTS", "")
S3_BASE     <- if (nzchar(s3_results_env)) {
    s3_results_env
} else if (nzchar(run_id_env)) {
    paste0("s3://neoantigen2026-rerun/results/", run_id_env)
} else {
    "s3://neoantigen2026-rerun/results"
}
S3_OUT      <- file.path(S3_BASE, "final")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

s3_fetch <- function(s3_path, local_path) {
    if (!file.exists(local_path))
        system(paste0("aws s3 cp \"", s3_path, "\" \"", local_path, "\""))
    file.exists(local_path)
}

# ── Load netMHCpan predictions ────────────────────────────────────────────────
mhc_file <- file.path(MHC_DIR, "all_predictions.tsv")
s3_fetch(file.path(S3_BASE, "netmhcpan/all_predictions.tsv"), mhc_file)
if (!file.exists(mhc_file)) stop("netMHCpan predictions not found. Run Step 9 first.")

mhc <- read_tsv(mhc_file, show_col_types = FALSE) %>%
    mutate(binder_class = case_when(
        el_rank < 0.5 ~ "SB",
        el_rank < 2.0 ~ "WB",
        TRUE           ~ "NB"
    )) %>%
    filter(binder_class != "NB")
message(sprintf("netMHCpan: %d binder rows (SB+WB)", nrow(mhc)))

# ── Load peptide metadata (WES + fusions combined) ────────────────────────────
load_pep_meta <- function() {
    wes <- file.path(PEP_DIR, "all_peptides.tsv")
    fus <- file.path(PEP_DIR, "fusions_all.tsv")
    s3_fetch(file.path(S3_BASE, "peptides/all_peptides.tsv"), wes)
    s3_fetch(file.path(S3_BASE, "peptides/fusions_all.tsv"), fus)
    dfs <- list()
    if (file.exists(wes)) dfs[["wes"]] <- read_tsv(wes, show_col_types = FALSE)
    if (file.exists(fus)) dfs[["fus"]] <- read_tsv(fus, show_col_types = FALSE)
    if (length(dfs) == 0) return(NULL)
    bind_rows(dfs) %>%
        distinct(peptide, .keep_all = TRUE) %>%
        mutate(
            mut_type = case_when(
                grepl("fusion",       effect, ignore.case = TRUE) ~ "fusion",
                grepl("frameshift",   effect, ignore.case = TRUE) ~ "frameshift",
                grepl("missense",     effect, ignore.case = TRUE) ~ "missense",
                grepl("stop|start",   effect, ignore.case = TRUE) ~ "stop/start",
                grepl("inframe",      effect, ignore.case = TRUE) ~ "inframe_indel",
                TRUE ~ "other"
            )
        )
}
pep_meta <- load_pep_meta()
if (!is.null(pep_meta)) {
    message(sprintf("Peptide metadata: %d unique peptides", nrow(pep_meta)))
    mhc <- mhc %>%
        left_join(
            pep_meta %>% select(peptide, sample, mut_id, gene, transcript,
                                effect, hgvsp, is_frameshift, mut_type) %>%
                         distinct(peptide, .keep_all = TRUE),
            by = "peptide"
        )
} else {
    message("[WARN] No peptide metadata — skipping join")
    mhc <- mhc %>% mutate(sample = NA, mut_id = NA, gene = NA,
                           effect = NA, hgvsp = NA, mut_type = NA,
                           is_frameshift = NA)
}

# ── Load PyClone clonality ────────────────────────────────────────────────────
loci_f <- file.path(PYCLONE_DIR, "loci_clusters.tsv")
cp_f   <- file.path(PYCLONE_DIR, "cluster_prevalence.tsv")
s3_fetch(file.path(S3_BASE, "pyclone/tables/loci_clusters.tsv"),     loci_f)
s3_fetch(file.path(S3_BASE, "pyclone/tables/cluster_prevalence.tsv"), cp_f)

if (file.exists(loci_f) && file.exists(cp_f)) {
    loci   <- read_tsv(loci_f, show_col_types = FALSE)
    cp     <- read_tsv(cp_f,   show_col_types = FALSE)
    max_cp <- cp %>%
        group_by(cluster_id) %>%
        summarise(max_cp = max(mean_cp, na.rm = TRUE), .groups = "drop") %>%
        mutate(clonality = ifelse(max_cp >= 0.5, "Clonal", "Subclonal"))
    mhc <- mhc %>%
        left_join(loci %>% select(mutation_id, cluster_id) %>% distinct(),
                  by = c("mut_id" = "mutation_id")) %>%
        left_join(max_cp, by = "cluster_id")
    message(sprintf("PyClone: %d clusters", n_distinct(mhc$cluster_id, na.rm = TRUE)))
} else {
    message("[WARN] PyClone data unavailable — clonality set to NA")
    mhc <- mhc %>% mutate(cluster_id = NA, max_cp = NA, clonality = NA)
}

# ── Load RNASeq expression ────────────────────────────────────────────────────
expr_f <- file.path(RNASEQ_DIR, "expression_matrix.tsv")
s3_fetch(file.path(S3_BASE, "rnaseq/expression_matrix.tsv"), expr_f)

if (file.exists(expr_f)) {
    expr_wide <- read_tsv(expr_f, show_col_types = FALSE)

    # Median expression per gene across tumour samples
    tumour_cols <- setdiff(names(expr_wide), c("gene_id",
                           "423_D0_old", "443_D21_new", "D88_old",
                           "D99_old", "D109_new"))
    expr_median <- expr_wide %>%
        mutate(expr_median = apply(
            select(., all_of(intersect(tumour_cols, names(.)))),
            1, median, na.rm = TRUE)) %>%
        select(gene_id, expr_median)

    # Join to MHC results by gene_id (ENSMUSG ID or gene symbol)
    mhc <- mhc %>%
        left_join(expr_median, by = c("gene" = "gene_id")) %>%
        mutate(expressed = !is.na(expr_median) & expr_median > 1.0)
    message(sprintf("Expression: %.1f%% of candidates have expression data",
                    100 * mean(!is.na(mhc$expr_median))))
} else {
    message("[WARN] Expression data unavailable — run Step 11 first")
    mhc <- mhc %>% mutate(expr_median = NA, expressed = NA)
}

# ── Priority scoring ──────────────────────────────────────────────────────────
# Binding strength: SB=3, WB=1
# Mutation type   : missense=2, frameshift/fusion=1  (lower EL rank for novel)
# Clonality       : Clonal +1
# Expression      : expressed +1
mhc <- mhc %>%
    mutate(
        binding_score   = ifelse(binder_class == "SB", 3L, 1L),
        mutation_score  = case_when(
            mut_type %in% c("missense")                            ~ 2L,
            mut_type %in% c("frameshift","fusion","stop/start")    ~ 1L,
            TRUE                                                   ~ 1L
        ),
        clonality_score = ifelse(!is.na(clonality) & clonality == "Clonal", 1L, 0L),
        expr_score      = ifelse(!is.na(expressed)  & expressed,             1L, 0L),
        priority_score  = binding_score + mutation_score + clonality_score + expr_score
    ) %>%
    arrange(desc(priority_score), el_rank)

# ── Write outputs ─────────────────────────────────────────────────────────────
write_tsv(mhc, file.path(OUT_DIR, "neoantigen_candidates.tsv"))
message(sprintf("Candidates: %d rows", nrow(mhc)))

top_sb <- mhc %>%
    filter(binder_class == "SB") %>%
    distinct(peptide, allele, .keep_all = TRUE)
write_tsv(top_sb, file.path(OUT_DIR, "neoantigen_top_ranked.tsv"))
message(sprintf("Top SB candidates: %d unique (peptide × allele)", nrow(top_sb)))

write_tsv(
    mhc %>% count(allele, plen, binder_class, mut_type, name = "n"),
    file.path(OUT_DIR, "binding_summary.tsv")
)

# ── Figures ───────────────────────────────────────────────────────────────────
n_sb  <- sum(mhc$binder_class == "SB")
n_wb  <- sum(mhc$binder_class == "WB")
pal_type <- c("missense"    = "#E41A1C",
              "frameshift"  = "#FF7F00",
              "fusion"      = "#984EA3",
              "stop/start"  = "#377EB8",
              "inframe_indel" = "#4DAF4A",
              "other"       = "grey60")

# Fig 1: Binders by allele, coloured by mutation type
fig1 <- mhc %>%
    mutate(mut_type = replace_na(mut_type, "other"),
           binder_class = factor(binder_class, c("SB","WB"))) %>%
    ggplot(aes(x = allele, fill = mut_type)) +
    geom_bar(colour = "white") +
    facet_wrap(~ binder_class, scales = "free_y") +
    scale_fill_manual(values = pal_type, name = "Mutation type") +
    labs(title = "Binder Counts by Allele and Mutation Type",
         x = "MHC Allele", y = "Count") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1))

# Fig 2: EL rank distribution by mutation type
fig2 <- mhc %>%
    mutate(mut_type = replace_na(mut_type, "other")) %>%
    ggplot(aes(x = el_rank, fill = mut_type)) +
    geom_histogram(bins = 50, colour = "white", alpha = 0.85) +
    facet_wrap(~ mut_type, scales = "free_y") +
    scale_fill_manual(values = pal_type) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
    labs(title = "EL Rank by Mutation Type (dashed = SB threshold)",
         x = "EL Rank (%)", y = "Count") +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold"))

# Fig 3: Clonality × expression × binding heatmap
fig3_data <- mhc %>%
    filter(!is.na(clonality), binder_class == "SB") %>%
    count(clonality, expressed = replace_na(expressed, FALSE),
          mut_type = replace_na(mut_type, "other"), name = "n")
if (nrow(fig3_data) > 0) {
    fig3 <- ggplot(fig3_data, aes(x = clonality, y = mut_type, fill = n)) +
        geom_tile(colour = "white") +
        facet_wrap(~ expressed,
                   labeller = labeller(expressed = c("TRUE" = "Expressed",
                                                      "FALSE" = "Not expressed"))) +
        scale_fill_gradient(low = "white", high = "#E41A1C", name = "Count") +
        labs(title = "Strong Binders: Clonality × Expression × Type",
             x = "Clonality", y = "Mutation type") +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold"))
} else {
    fig3 <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Awaiting PyClone data", size = 5) +
        theme_void() + labs(title = "Clonality × Expression (pending PyClone)")
}

# Fig 4: Expression vs EL rank scatter (top 200)
if (any(!is.na(mhc$expr_median))) {
    fig4 <- mhc %>%
        filter(!is.na(expr_median)) %>%
        slice_head(n = 200) %>%
        mutate(mut_type = replace_na(mut_type, "other")) %>%
        ggplot(aes(x = expr_median, y = -log10(el_rank + 0.01),
                   colour = mut_type, shape = binder_class)) +
        geom_point(alpha = 0.7, size = 2) +
        geom_hline(yintercept = -log10(0.5), linetype = "dashed", colour = "grey40") +
        scale_colour_manual(values = pal_type) +
        labs(title = "Expression vs Binding Affinity (top 200)",
             x = "Median tumour expression (log2)", y = "-log10(EL Rank)",
             colour = "Type", shape = "Binder") +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold"))
} else {
    fig4 <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Expression data pending (Step 11)", size = 5) +
        theme_void() + labs(title = "Expression vs Binding (pending)")
}

# Fig 5: Top 20 priority peptides
top20 <- top_sb %>%
    arrange(desc(priority_score), el_rank) %>%
    slice_head(n = 20) %>%
    mutate(
        label      = paste0(replace_na(gene,"?"), "\n", peptide),
        mut_type   = replace_na(mut_type, "other"),
        clonality  = replace_na(as.character(clonality), "unknown")
    )
fig5 <- ggplot(top20,
               aes(x = allele,
                   y = reorder(label, priority_score),
                   colour = mut_type, size = priority_score)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(2, 7)) +
    scale_colour_manual(values = pal_type) +
    labs(title = "Top 20 Strong Binder Candidates",
         x = "Allele", y = "Peptide",
         colour = "Type", size = "Priority") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 7))

# Fig 6: Priority score distribution
fig6 <- mhc %>%
    mutate(binder_class = factor(binder_class, c("SB","WB"))) %>%
    ggplot(aes(x = factor(priority_score), fill = binder_class)) +
    geom_bar(colour = "white") +
    scale_fill_manual(values = c("SB" = "#E41A1C", "WB" = "#FF7F00")) +
    labs(title = "Priority Score Distribution",
         x = "Priority Score", y = "Count", fill = "Binder") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

combined <- (fig1 | fig2) / (fig3 | fig4) / (fig5 | fig6) +
    plot_annotation(
        title    = "Neoantigen Candidate Report",
        subtitle = sprintf(
            "SB=%d  WB=%d  |  Missense + Frameshift + Fusions  |  %s",
            n_sb, n_wb, format(Sys.time(), "%Y-%m-%d %H:%M")),
        theme = theme(
            plot.title    = element_text(size = 14, face = "bold"),
            plot.subtitle = element_text(size = 9,  colour = "grey40"))
    )

ggsave(file.path(OUT_DIR, "neoantigen_report.pdf"),
       combined, width = 18, height = 20, device = "pdf")
ggsave(file.path(OUT_DIR, "neoantigen_report.png"),
       combined, width = 18, height = 20, dpi = 150)
message("Saved: neoantigen_report.pdf/png")

# ── Upload ─────────────────────────────────────────────────────────────────────
system(paste("aws s3 sync", OUT_DIR, S3_OUT))
message("Uploaded to S3: ", S3_OUT)

# ── Console summary ────────────────────────────────────────────────────────────
message("══════════════════════════════════════════════════════")
message(sprintf("  Strong Binders (SB, EL<0.5%%):   %d", n_sb))
message(sprintf("  Weak Binders   (WB, EL<2.0%%):   %d", n_wb))
if (!all(is.na(mhc$mut_type))) {
    for (mt in c("missense","frameshift","fusion","stop/start","inframe_indel")) {
        n_mt <- sum(mhc$mut_type == mt & mhc$binder_class == "SB", na.rm = TRUE)
        if (n_mt > 0) message(sprintf("    SB %-18s %d", paste0(mt,":"), n_mt))
    }
}
if (!all(is.na(mhc$clonality))) {
    message(sprintf("  Clonal SB binders:             %d",
                    sum(mhc$binder_class == "SB" & mhc$clonality == "Clonal", na.rm = TRUE)))
}
if (!all(is.na(mhc$expressed))) {
    message(sprintf("  Expressed SB binders:          %d",
                    sum(mhc$binder_class == "SB" & mhc$expressed, na.rm = TRUE)))
}
message(sprintf("  Output: %s", OUT_DIR))
message("══════════════════════════════════════════════════════")
message("Step 10 complete.")
