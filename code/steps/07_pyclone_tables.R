#!/usr/bin/env Rscript
# =============================================================================
# Step 7 — PyClone Cluster and Loci Tables + Figures
# Reads PyClone trace output and builds:
#   - cluster summary table (cellular prevalence per timepoint)
#   - loci assignment table (mutation → cluster)
#   - clonal evolution figures
# =============================================================================

local_lib <- path.expand("~/R/library")
dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

for (pkg in c("ggplot2","dplyr","tidyr","readr","patchwork","RColorBrewer",
              "scales","ggrepel")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg, repos="https://cloud.r-project.org",
                         lib=local_lib, quiet=TRUE)
    suppressPackageStartupMessages(library(pkg, character.only=TRUE))
}

PYCLONE_DIR <- "/home/ec2-user/results/pyclone"
OUT_DIR     <- "/home/ec2-user/results/pyclone/tables"
S3_OUT      <- "s3://bam-wes/NeoAntigen-aws/results/pyclone/tables"
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

TIMEPOINTS <- c("428_D20_new"=20, "34_D52_old"=52,
                "36_D99_new"=99,  "38_D99_new"=99, "42_D122_old"=122)

# ── Load PyClone output ───────────────────────────────────────────────────────
loci_f    <- file.path(PYCLONE_DIR, "output/tables/loci.tsv")
cluster_f <- file.path(PYCLONE_DIR, "output/tables/cluster.tsv")

if (!file.exists(loci_f)) {
    # Try S3
    system(paste("aws s3 cp s3://bam-wes/NeoAntigen-aws/results/pyclone/output/tables/loci.tsv",    loci_f))
    system(paste("aws s3 cp s3://bam-wes/NeoAntigen-aws/results/pyclone/output/tables/cluster.tsv", cluster_f))
}

if (!file.exists(loci_f)) {
    stop("PyClone loci table not found. Run Step 6 first.")
}

loci    <- read_tsv(loci_f,    show_col_types=FALSE)
cluster <- read_tsv(cluster_f, show_col_types=FALSE)

message(sprintf("Loci table:    %d rows, %d cols", nrow(loci), ncol(loci)))
message(sprintf("Cluster table: %d rows, %d cols", nrow(cluster), ncol(cluster)))

# ── Cluster summary table ─────────────────────────────────────────────────────
clust_summary <- cluster %>%
    group_by(cluster_id) %>%
    summarise(
        n_mutations = n_distinct(mutation_id),
        .groups = "drop"
    )

# Cellular prevalence per sample per cluster (mean of posterior)
cp_cols <- grep("cellular_prevalence$", names(loci), value=TRUE)
if (length(cp_cols) > 0) {
    cp_long <- loci %>%
        select(mutation_id, cluster_id, all_of(cp_cols)) %>%
        pivot_longer(all_of(cp_cols), names_to="sample", values_to="cp") %>%
        mutate(sample = gsub("_cellular_prevalence", "", sample),
               timepoint = TIMEPOINTS[sample]) %>%
        group_by(cluster_id, sample, timepoint) %>%
        summarise(mean_cp = mean(cp, na.rm=TRUE),
                  sd_cp   = sd(cp,   na.rm=TRUE), .groups="drop")

    write_tsv(cp_long, file.path(OUT_DIR, "cluster_prevalence.tsv"))
    message("Saved: cluster_prevalence.tsv")
} else {
    cp_long <- NULL
    message("[WARN] No cellular prevalence columns found in loci table")
}

# Loci assignment table
loci_out <- loci %>%
    select(mutation_id, cluster_id, any_of(c("gene","effect","aa_change",cp_cols)))
write_tsv(loci_out, file.path(OUT_DIR, "loci_clusters.tsv"))

write_tsv(clust_summary, file.path(OUT_DIR, "cluster_summary.tsv"))
message("Saved: cluster_summary.tsv, loci_clusters.tsv")

# ── Figures ───────────────────────────────────────────────────────────────────
n_clusters <- n_distinct(loci$cluster_id)
pal <- colorRampPalette(brewer.pal(min(9, n_clusters), "Set1"))(n_clusters)
names(pal) <- sort(unique(loci$cluster_id))

# Fig 1: Cluster sizes (mutations per cluster)
fig1 <- clust_summary %>%
    mutate(cluster_id = factor(cluster_id)) %>%
    ggplot(aes(x=reorder(cluster_id,-n_mutations), y=n_mutations,
               fill=cluster_id)) +
    geom_col(colour="white") +
    scale_fill_manual(values=pal) +
    labs(title="Mutations per Cluster", x="Cluster", y="Number of Mutations") +
    theme_classic(base_size=12) +
    theme(legend.position="none", plot.title=element_text(face="bold"))

# Fig 2: Clonal evolution — cellular prevalence over time
if (!is.null(cp_long)) {
    fig2 <- cp_long %>%
        filter(!is.na(timepoint)) %>%
        mutate(cluster_id=factor(cluster_id)) %>%
        ggplot(aes(x=timepoint, y=mean_cp, colour=cluster_id, group=cluster_id)) +
        geom_ribbon(aes(ymin=pmax(0,mean_cp-sd_cp),
                        ymax=pmin(1,mean_cp+sd_cp),
                        fill=cluster_id), alpha=0.15, colour=NA) +
        geom_line(linewidth=1.2) +
        geom_point(size=3) +
        scale_colour_manual(values=pal, name="Cluster") +
        scale_fill_manual(values=pal,   name="Cluster") +
        scale_x_continuous(breaks=c(20,52,99,122),
                           labels=paste0("D",c(20,52,99,122))) +
        scale_y_continuous(limits=c(0,1), labels=percent) +
        labs(title="Clonal Evolution Over Time",
             subtitle="Mean posterior cellular prevalence ± SD",
             x="Timepoint", y="Cellular Prevalence") +
        theme_classic(base_size=12) +
        theme(plot.title=element_text(face="bold"))
} else {
    fig2 <- ggplot() + annotate("text",x=0.5,y=0.5,label="No CP data") +
        theme_void() + labs(title="Cellular Prevalence (unavailable)")
}

# Fig 3: VAF distribution per cluster
vaf_cols <- grep("variant_allele_frequency|vaf", names(loci), ignore.case=TRUE, value=TRUE)
if (length(vaf_cols) > 0) {
    fig3 <- loci %>%
        select(mutation_id, cluster_id, all_of(vaf_cols[1])) %>%
        rename(vaf=all_of(vaf_cols[1])) %>%
        mutate(cluster_id=factor(cluster_id)) %>%
        ggplot(aes(x=vaf, fill=cluster_id)) +
        geom_histogram(bins=40, colour="white", alpha=0.8) +
        facet_wrap(~cluster_id, scales="free_y", ncol=3) +
        scale_fill_manual(values=pal) +
        scale_x_continuous(labels=percent) +
        labs(title="VAF Distribution per Cluster", x="VAF", y="Count") +
        theme_classic(base_size=11) + theme(legend.position="none",
                                            plot.title=element_text(face="bold"))
} else {
    fig3 <- ggplot() + annotate("text",x=0.5,y=0.5,label="No VAF data") +
        theme_void()
}

combined <- (fig1 | fig2) / fig3 +
    plot_annotation(
        title    = "PyClone Clonal Analysis",
        subtitle = sprintf("Multi-sample beta-binomial MCMC | %d clusters | Generated: %s",
                           n_clusters, format(Sys.time(),"%Y-%m-%d %H:%M")),
        theme=theme(plot.title=element_text(size=15,face="bold"),
                    plot.subtitle=element_text(size=9,colour="grey40"))
    )

ggsave(file.path(OUT_DIR,"pyclone_report.pdf"), combined, width=16, height=13, device="pdf")
ggsave(file.path(OUT_DIR,"pyclone_report.png"), combined, width=16, height=13, dpi=150)
message("Saved: pyclone_report.pdf/png")

system(paste("aws s3 sync", OUT_DIR, S3_OUT))
message("Step 7 complete.")
