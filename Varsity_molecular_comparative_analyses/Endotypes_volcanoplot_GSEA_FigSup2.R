# Endotypes_volcanoplot_GSEA_FigSup2.R 
# Objective: visualize trajectory DEGs (NR vs responder classes) with volcano plots
# and Hallmark fgsea enrichment used in the Endotypes manuscript

# This code is corresponding to Extended Figure 2

# ------------------------------------------------------------------
# Libraries (kept aligned with the Rmd for reproducibility)
# ------------------------------------------------------------------
library(dplyr)
library(reshape)
library(tidyr)
library(stringr)
library(devtools)
library(DESeq2)
library(NMF)
library(matrixStats)
library(ggfortify)
library(EnhancedVolcano)
library(gridExtra)
library(sva)
library(readxl)
library(edgeR)
library(EnsDb.Hsapiens.v79)
library(ggplot2)
library(ggvenn)
library(ComplexHeatmap)
library(biomaRt)
library(GSVA)
library(pheatmap)
library(RColorBrewer)
library(circlize)
library(readr)
library(GGally)
library(gprofiler2)
scale_min_max <- function(x) { (x - 0) / (max(x) - 0) }
z_scale <- function(x) { (x - mean(x)) / sd(x) }
library(msigdbr)
library(fgsea)
library(glue)
library(Cairo)

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------
# Helper functions 
# ------------------------------------------------------------------

# Build the DEG file path (adjust deg_dir if files move)
deg_filepath <- function(treatment, timepoint, control, deg_dir = "DGE_results/Trajectory_DEGs") {
  glue("{deg_dir}/DEG_all_{treatment}_{timepoint}_non_respondervs{control}.txt")
}

# Load a DEG table, de-duplicate by GeneName, and set rownames
load_deg_table <- function(treatment, timepoint, control, deg_dir = "DGE_results/Trajectory_DEGs") {
  path <- deg_filepath(treatment, timepoint, control, deg_dir)
  if (!file.exists(path)) stop("DEG file not found: ", path)
  read.delim(path, sep = ";", check.names = FALSE) %>%
    dplyr::filter(!duplicated(GeneName)) %>%
    tibble::column_to_rownames("GeneName")
}

# Count up/down genes using raw p-values and adjusted p-values
count_degs <- function(deg_df, p_cutoff = 0.05, padj_cutoff = 0.05) {
  list(
    up_raw   = nrow(deg_df %>% dplyr::filter(log2FoldChange > 0, pvalue < p_cutoff)),
    down_raw = nrow(deg_df %>% dplyr::filter(log2FoldChange < 0, pvalue < p_cutoff)),
    up_adj   = nrow(deg_df %>% dplyr::filter(log2FoldChange > 0, padj   < padj_cutoff)),
    down_adj = nrow(deg_df %>% dplyr::filter(log2FoldChange < 0, padj   < padj_cutoff))
  )
}

# Build the volcano plot; parameters mirror the Rmd chunk
make_volcano_plot <- function(deg_df,
                              label_genes,
                              title_text,
                              fc_cutoff = 1,
                              p_cutoff = 0.05,
                              x_limits = c(-5, 7.5),
                              y_limits = c(0, 8),
                              counts = NULL) {
  plot_obj <- EnhancedVolcano(
    deg_df,
    lab = rownames(deg_df),
    x = "log2FoldChange",
    y = "pvalue",
    ylim = y_limits,
    xlim = x_limits,
    selectLab = label_genes,
    xlab = bquote(~Log[2]~ "fold change"),
    ylab = bquote(~-Log[10]~ "pvalue"),
    pCutoff = p_cutoff,
    FCcutoff = fc_cutoff,
    pointSize = 2.0,
    labSize = 3,
    title = title_text,
    subtitle = "",
    caption = "",
    boxedLabels = FALSE,
    colAlpha = 3/5,
    legendPosition = "None",
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = "black",
    max.overlaps = 30
  )
  if (!is.null(counts)) {
    plot_obj <- plot_obj +
      annotate(
        "text", x = Inf, y = Inf,
        label = sprintf(
          "Up (raw): %d\nDown (raw): %d\nUp (adj): %d\nDown (adj): %d",
          counts$up_raw, counts$down_raw, counts$up_adj, counts$down_adj
        ),
        hjust = 1.1, vjust = 1.1, size = 3.5, color = "black"
      )
  }
  plot_obj
}

# ------------------------------------------------------------------
# Volcano plot for a single comparison 
# ------------------------------------------------------------------

treatment <- "Adalimumab"
timepoint <- "Visit 9"
control <- "partial_responder"

DEG_all <- load_deg_table(treatment, timepoint, control)

pp <- read.csv("./files/hallmark_pathways_genes_all.csv") %>%
  dplyr::filter(Hallmark_Name == "ANGIOGENESIS")

angiogenesis_all <- lapply(pp$Gene_List, function(txt) {
  strsplit(txt, ",")[[1]] %>% trimws()
})[[1]]

agniogeensis_filtered <- c(
  "SPP1", "STC1", "TIMP1", "CXCL6", "OLR1", "THBD", "PGLYRP1", "PRG2",
  "FSTL1", "KCNJ8", "FGFR1", "PF4", "CCND2", "MSX1"
)

counts_single <- count_degs(DEG_all, p_cutoff = 0.05, padj_cutoff = 0.05)

volcano_plot <- make_volcano_plot(
  DEG_all,
  label_genes = c(agniogeensis_filtered, "ANGPT2"),
  title_text = glue("NR vs {control} for {treatment} at {timepoint}"),
  fc_cutoff = 1,
  p_cutoff = 0.05,
  x_limits = c(-5, 7.5),
  y_limits = c(0, 8),
  counts = counts_single
)

print(volcano_plot)

# Export (same code as Rmd; leave commented to keep behavior unchanged)
# CairoPDF(file = glue("Figs/Manuscript_EPS/FIG2A_{treatment}_{timepoint}_{control}_Volcano_allsamples_W14.pdf"), width = 10, height = 10)
# plot(volcano_plot)
# dev.off()
# cairo_ps(glue("Figs/Manuscript_EPS/FIG2A_{treatment}_{timepoint}_{control}_Volcano_allsamples_W14.eps"), width = 10, height = 10)
# print(volcano_plot)
# dev.off()

# ------------------------------------------------------------------
# Iterative volcano plotting over treatment/timepoint/control grids
# ------------------------------------------------------------------

pp <- read.csv("./files/hallmark_pathways_genes_all.csv") %>%
  dplyr::filter(Hallmark_Name == "ANGIOGENESIS")

angiogenesis_all <- lapply(pp$Gene_List, function(txt) {
  strsplit(txt, ",")[[1]] %>% trimws()
})[[1]]

agniogeensis_filtered <- c(
  "SPP1", "STC1", "TIMP1", "CXCL6", "OLR1", "THBD", "PGLYRP1", "PRG2",
  "FSTL1", "KCNJ8", "FGFR1", "PF4", "CCND2", "MSX1"
)

treatments <- c("Adalimumab", "Vedolizumab")
timepoints <- c("Visit 1", "Visit 9")
controls <- c("partial_responder", "responder", "super_responder")

for (treatment in treatments) {
  for (timepoint in timepoints) {
    for (control in controls) {
      fn <- deg_filepath(treatment, timepoint, control)
      if (!file.exists(fn)) {
        warning("File not found: ", fn)
        next
      }
      DEG_all <- load_deg_table(treatment, timepoint, control)
      counts_loop <- count_degs(DEG_all, p_cutoff = 0.05, padj_cutoff = 0.05)

      volcano_plot <- make_volcano_plot(
        DEG_all,
        label_genes = c(agniogeensis_filtered, "ANGPT2"),
        title_text = glue("NR vs {control} for {treatment} at {timepoint}"),
        fc_cutoff = 1,
        p_cutoff = 0.05,
        x_limits = c(-5, 7.5),
        y_limits = c(0, 8),
        counts = counts_loop
      )

      base <- glue("{treatment}_{timepoint}_{control}")

      # Uncomment to export
      # out_pdf <- glue("Figs/Manuscript_EPS/Fig2_new/FIG2A_{base}_Volcano_allsamples_W14.pdf")
      # CairoPDF(file = out_pdf, width = 10, height = 10)
      # print(volcano_plot)
      # dev.off()
      #
      # out_eps <- glue("Figs/Manuscript_EPS/Fig2_new/FIG2A_{base}_Volcano_allsamples_W14.eps")
      # cairo_ps(file = out_eps, width = 10, height = 10)
      # print(volcano_plot)
      # dev.off()

      message("Done: ", base)
    }
  }
}

# ------------------------------------------------------------------
# Functional analysis using MSigDB - Set the parameters accordingly (e.g. Treatment, timepoint and control)
# ------------------------------------------------------------------

treatment <- "Adalimumab" #"Vedolizumab"
timepoint <- "Visit 9" #"Visit 1"; "Visit 9";"Visit 28";
control <- "responder"

DEG_results <- read.delim(
  deg_filepath(treatment, timepoint, control),
  sep = ";"
)

DE <- DEG_results %>%
  dplyr::arrange(desc(log2FoldChange))

gene_list <- DE$log2FoldChange
names(gene_list) <- DE$GeneName
gene_list <- sort(gene_list, decreasing = TRUE)

msigdbr_species <- "Homo sapiens"
msigdbr_category <- "H"  # Hallmark gene sets

msigdbr_gene_sets <- msigdbr(species = msigdbr_species, category = msigdbr_category)
msigdbr_list <- split(msigdbr_gene_sets$gene_symbol, msigdbr_gene_sets$gs_name)

fgsea_results <- fgsea(pathways = msigdbr_list, stats = gene_list, nperm = 1000)
fgsea_results$leadingEdge <- NULL
fgsea_results$pathway <- gsub("^HALLMARK_", "", fgsea_results$pathway)

write.table(
  fgsea_results,
  glue("GSEA_results/GSEA_tables/VARSITY_GSEA_Final/{treatment}_{timepoint}_NR_vs_{control}_allsamples_W14.csv"),
  sep = ",",
  row.names = FALSE,
  col.names = TRUE
)

fgsea_df <- fgsea_results %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.1) %>%
  dplyr::mutate(pathway = reorder(pathway, NES))

dotplot <- ggplot(fgsea_df, aes(x = NES, y = pathway)) +
  geom_point(aes(size = -log10(padj), color = padj)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(
    title = "GSEA Enrichment Analysis",
    x = "Normalized Enrichment Score (NES)",
    y = "Pathway",
    size = "-log10(Adjusted P-value)",
    color = "Adjusted P-value"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(dotplot)

# Uncomment to export plot
# png(glue("~/Varsity/Trajectory_analysis/Figs/Volcano_heatmap_allsamplesW14/{treatment}_{timepoint}_{control}_dotplot_allsamples_W14.png"), width = 1200, height = 1000, res = 150)
# print(dotplot)
# dev.off()
