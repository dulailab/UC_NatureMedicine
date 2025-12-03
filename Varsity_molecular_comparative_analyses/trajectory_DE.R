# Trajectory analysis (R script version of trajectory_DE.Rmd)
# Objective: perform PCA/DESeq2-based trajectory DEs, volcano summaries,
# cohorts; logic/order matches the Rmd with configurable paths.
# Code order mirrors the Rmd; logic is unchanged. Paths for inputs/outputs are
# defined up front for easier reuse when sharing.

# ------------------------------------------------------------------
# Directory and file paths (adjust as needed)
# ------------------------------------------------------------------
project_dir <- getwd()
files_dir <- file.path(project_dir, "files")
trajectory_deg_dir <- file.path(project_dir, "Extended Data Fig2", "Trajectory_DEGs")
de_results_dir <- file.path(project_dir, "DE_results")

trajectory_labels_file <- file.path(files_dir, "Jie's labels.csv")

sample_metadata_filtered_file <- file.path(files_dir, "RNAseq_data/sample_metadata_filtered.txt")
varsity_fpkm_file <- file.path(files_dir, "RNAseq_data/varsity_FPKM.txt")
varsity_count_file <- file.path(files_dir, "RNAseq_data/varsity_count.txt")
varsity_geneannot_file <- file.path(files_dir, "RNAseq_data/gene_annotation.txt")

# dir.create(de_results_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(fe_results_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(manuscript_eps_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create(fig2_new_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------
# Libraries (same set as Rmd)
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
library(org.Hs.eg.db)
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
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------
# Loading expression, annotation and meta data
# ------------------------------------------------------------------
sample.metadata.filtered <- read.delim(file = sample_metadata_filtered_file, sep = "\t", header = TRUE)
varsity.count <- read.delim(file = varsity_count_file, sep = "\t", header = TRUE)
varsity.FPKM <- read.delim(file = varsity_fpkm_file, sep = "\t", header = TRUE)
gene.annotation <- read.delim(file = varsity_geneannot_file, sep = "\t", header = TRUE)

varsity.count$GeneID <- rownames(varsity.count)
varsity.count <- varsity.count %>%
  dplyr::left_join(gene.annotation[, c("GeneID", "GeneName")], by = "GeneID") %>%
  dplyr::filter(!duplicated(GeneName))
rownames(varsity.count) <- varsity.count$GeneName
varsity.count$GeneID <- NULL
varsity.count$GeneName <- NULL

# removing the duplicates
sample.metadata.filtered <- sample.metadata.filtered %>% dplyr::filter(!duplicated(Originating_ID))
rownames(sample.metadata.filtered) <- sample.metadata.filtered$Originating_ID

# ------------------------------------------------------------------
# Base filtering by treatment/timepoint and trajectory labels
# ------------------------------------------------------------------
treatment <- "Vedolizumab"
visit <- "Visit 1"

sample.metadata.vedo.BL <- sample.metadata.filtered %>% dplyr::filter(Treatment == treatment, Visit_Type == visit)
varsity.count.vedo.BL <- varsity.count[, sample.metadata.vedo.BL$Originating_ID]
varsity.FPKM.vedo.BL <- varsity.FPKM[, sample.metadata.vedo.BL$Originating_ID]
all(colnames(varsity.count.vedo.BL) == sample.metadata.vedo.BL$Originating_ID)

vedo_trajectory_jie <- read.delim(trajectory_labels_file, sep = ",", header = TRUE)

sample.metadata.vedo.BL <- sample.metadata.vedo.BL %>% dplyr::left_join(vedo_trajectory_jie, by = c("SubjectID" = "SUBJID"))
table(sample.metadata.vedo.BL$trajectory)


# ------------------------------------------------------------------
# DEG analysis
# ------------------------------------------------------------------
str(sample.metadata.vedo.BL)

comparison <- c("non_responder", "super_responder")

t_meta <- sample.metadata.vedo.BL %>% dplyr::filter(trajectory %in% comparison)
pseudoCount <- log2(varsity.FPKM.vedo.BL + 1)
t_pseudoCount <- pseudoCount[, as.character(t_meta$Sample)]
t_count <- varsity.count.vedo.BL[, as.character(t_meta$Sample)]

dds <- DESeqDataSetFromMatrix(
  countData = round(t_count),
  colData = t_meta,
  design = ~trajectory
)

Deseq_dd <- dds
Deseq_dd$label <- factor(Deseq_dd$trajectory, levels = comparison)

gene.annotation <- gene.annotation %>% dplyr::filter(!duplicated(GeneName))
rownames(gene.annotation) <- gene.annotation$GeneName
expression.count.annot <- gene.annotation[rownames(t_count), ]

mcols(Deseq_dd) <- expression.count.annot
Deseq_dd <- DESeq(Deseq_dd)

DEG.comparison <- results(Deseq_dd, contrast = c("trajectory", comparison), alpha = 0.05)
summary(DEG.comparison)
DEG_all <- as.data.frame(DEG.comparison)
DEG.comparison_filtered <- DEG_all %>% dplyr::filter(padj < 0.05)

pos <- nrow(DEG.comparison_filtered[DEG.comparison_filtered$log2FoldChange > 0, ])
neg <- nrow(DEG.comparison_filtered[DEG.comparison_filtered$log2FoldChange < 0, ])
rest <- nrow(t_count) - (pos + neg)

DEG.df <- data.frame(
  category = c("Up Regulated", "Down Regulated", "Not Significant"),
  count = c(pos, neg, rest)
)
DEG.df$fraction <- DEG.df$count / sum(DEG.df$count)
DEG.df$ymax <- cumsum(DEG.df$fraction)
DEG.df$ymin <- c(0, head(DEG.df$ymax, n = -1))
DEG.df$labelPosition <- (DEG.df$ymax + DEG.df$ymin) / 2
DEG.df$label <- paste0(DEG.df$category, "\n %", round(DEG.df$fraction * 100, 2), "\n value: ", DEG.df$count)

mycols <- c("#0073C2FF", "grey68", "#CD534CFF")

ggplot(DEG.df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = mycols)) +
  geom_rect() +
  geom_label(x = 3.5, aes(y = labelPosition, label = label), size = 3.5) +
  scale_fill_brewer(palette = 4) +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#BE2A3E", "grey68", "skyblue3"))

DEG_all$GeneName <- rownames(DEG_all)
DEG_all <- DEG_all %>% left_join(gene.annotation, by = "GeneName")
DEG_all$GeneName <- toupper(DEG_all$GeneName)
DEG_all <- DEG_all %>% arrange(padj)

DEG.table <- as.data.frame(DEG_all)
DEG.table <- DEG.table %>% arrange(padj)

annots <- DEG.table$GeneName
EnhancedVolcano(
  DEG.table,
  lab = annots,
  x = "log2FoldChange",
  y = "padj",
  selectLab = gene.annotation[DEG.table$GeneName[1:10], "GeneName"],
  xlab = bquote(~Log[2]~ "fold change"),
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3,
  title = paste(comparison, collapse = " vs "),
  subtitle = "",
  boxedLabels = FALSE,
  colAlpha = 3 / 5,
  legendPosition = "bottom",
  legendLabSize = 14,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "black",
  max.overlaps = 30
)

# write.table(
#   DEG_all,
#   file.path(de_results_dir, paste0("DEG_all_", treatment, "_", visit, "_", str_replace(comparison[1], "/", "_"), "vs", str_replace(comparison[2], "/", "_"), ".txt")),
#   sep = ";",
#   col.names = TRUE,
#   row.names = FALSE
# )

