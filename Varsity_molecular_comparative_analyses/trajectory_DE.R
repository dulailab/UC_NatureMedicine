# Trajectory analysis (R script version of trajectory_DE.Rmd)
# Objective: perform PCA/DESeq2-based trajectory DEs, volcano/GSEA summaries,
# DEG intersections, LRT clustering, and related plots across Varsity RNA-seq
# cohorts; logic/order matches the Rmd with configurable paths.
# Code order mirrors the Rmd; logic is unchanged. Paths for inputs/outputs are
# defined up front for easier reuse when sharing.

# ------------------------------------------------------------------
# Directory and file paths (adjust as needed)
# ------------------------------------------------------------------
project_dir <- getwd()
files_dir <- file.path(project_dir, "files")
trajectory_deg_dir <- file.path(project_dir, "Extended Data Fig2", "Trajectory_DEGs")
trajectory_analysis_dir <- file.path(project_dir, "Trajectory_analysis")
de_results_dir <- file.path(project_dir, "DE_results")
fe_results_dir <- file.path(project_dir, "FE_results")
figs_dir <- file.path(project_dir, "Figs")
manuscript_eps_dir <- file.path(figs_dir, "Manuscript_EPS")
fig2_new_dir <- file.path(manuscript_eps_dir, "Fig2_new")

trajectory_labels_file <- file.path(trajectory_analysis_dir, "files", "Jie's labels.csv")
trajectory_labels_file_alt <- file.path(trajectory_analysis_dir, "Jie's labels.csv")

sample_metadata_complete_file <- file.path(files_dir, "sample_metadata_complete.txt")
sample_metadata_filtered_file <- file.path(files_dir, "sample_metadata_filtered.txt")
clinical_data_file <- file.path(files_dir, "clinical_data.RDS")
varsity_data_file <- file.path(files_dir, "Varsity_data.RDS")
admayo_file <- file.path(files_dir, "admayo1.csv")
disease_control_file <- file.path(project_dir, "disease control.xlsx")

dir.create(de_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fe_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_eps_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig2_new_dir, recursive = TRUE, showWarnings = FALSE)

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
if (file.exists(source_functions_path_to_use)) source(source_functions_path_to_use)
library(msigdbr)
library(fgsea)
library(Cairo)
options(stringsAsFactors = FALSE)

# ------------------------------------------------------------------
# Loading expression, annotation and meta data
# ------------------------------------------------------------------
sample.metadata <- read.delim(file = sample_metadata_complete_file, sep = ";", header = TRUE)
sample.metadata.filtered <- read.delim(file = sample_metadata_filtered_file, sep = ";", header = TRUE)
clinical_data <- readRDS(clinical_data_file)
data <- readRDS(varsity_data_file)
varsity.count <- data[[1]]
varsity.FPKM <- data[[2]]
gene.annotation <- data[[3]]

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

trajectory_labels_path <- if (file.exists(trajectory_labels_file)) trajectory_labels_file else trajectory_labels_file_alt
vedo_trajectory_jie <- read.delim(trajectory_labels_path, sep = ",", header = TRUE)

sample.metadata.vedo.BL <- sample.metadata.vedo.BL %>% dplyr::left_join(vedo_trajectory_jie, by = c("SubjectID" = "SUBJID"))
table(sample.metadata.vedo.BL$trajectory)

# ------------------------------------------------------------------
# PC analysis and clustering of the groups
# ------------------------------------------------------------------
comparison <- c("non_responder", "super_responder")

pseudoCount <- log2(varsity.FPKM.vedo.BL + 1)
target_samples <- sample.metadata.vedo.BL$Originating_ID
rownames(sample.metadata.vedo.BL) <- sample.metadata.vedo.BL$Originating_ID
sample.metadata.vedo.BL$Sample <- rownames(sample.metadata.vedo.BL)

t_meta <- sample.metadata.vedo.BL %>% dplyr::filter(Originating_ID %in% target_samples, trajectory %in% comparison)
t_pseudoCount <- pseudoCount[, as.character(t_meta$Sample)]

draw_CIM(t_pseudoCount, t_meta, c("trajectory"))

n <- 5000
p <- list()
PCs <- c(1, 2, 3, 4)
ind <- 0
for (i in c(1:(length(PCs) - 1))) {
  for (j in c((i + 1):length(PCs))) {
    ind <- ind + 1
    p[[ind]] <- draw_pca(t_pseudoCount, t_meta, x = i, y = j, n, color = "trajectory", lbl = FALSE, lbl_size = 1)
  }
}

do.call("grid.arrange", c(p[1:4], ncol = 2))
p[[1]]

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

write.table(
  DEG_all,
  file.path(de_results_dir, paste0("DEG_all_", treatment, "_", visit, "_", str_replace(comparison[1], "/", "_"), "vs", str_replace(comparison[2], "/", "_"), ".txt")),
  sep = ";",
  col.names = TRUE,
  row.names = FALSE
)

# ------------------------------------------------------------------
# Functional analysis using fgsea
# ------------------------------------------------------------------
DE_results <- read.delim(file.path(de_results_dir, "DEG_all_Vedolizumab_Visit 1_non_respondervspartial_responder.txt"), sep = ";", header = TRUE)

ensembl <- useEnsembl(biomart = "genes")

organism <- "hsapiens"
mart_dataset <- "hsapiens_gene_ensembl"
query_id <- "ensembl_gene_id"

n <- 2000
cutoff <- 0.05
DE_results <- DE_results %>% dplyr::arrange(padj) %>% head(n)
enrichment_results <- functional_enrichment_dplot(DE_results, organism, mart_dataset, cutoff = 0.05, return_plots = TRUE, activity = "both", n = 40)
enrichment_results[[2]]

# ------------------------------------------------------------------
# Intersection of DEGs
# ------------------------------------------------------------------
Vedo_NR_SR <- read.delim(file.path(de_results_dir, "DEG_all_Adalimumab_Visit 1_non_respondervspartial_responder.txt"), sep = ";", header = TRUE)
Vedo_NR_R <- read.delim(file.path(de_results_dir, "DEG_all_Adalimumab_Visit 1_non_respondervsresponder.txt"), sep = ";", header = TRUE)
Vedo_NR_PR <- read.delim(file.path(de_results_dir, "DEG_all_Adalimumab_Visit 1_non_respondervssuper_responder.txt"), sep = ";", header = TRUE)

Vedo_NR_SR.filter <- Vedo_NR_SR %>% dplyr::filter(padj < 0.05, Source == "protein_coding")
Vedo_NR_R.filter <- Vedo_NR_R %>% dplyr::filter(padj < 0.05, Source == "protein_coding")
Vedo_NR_PR.filter <- Vedo_NR_PR %>% dplyr::filter(padj < 0.05, Source == "protein_coding")

x <- list(
  Ada_NR_SR = Vedo_NR_SR.filter$GeneID,
  Ada_NR_R = Vedo_NR_R.filter$GeneID,
  Ada_NR_PR = Vedo_NR_PR.filter$GeneID
)

ggvenn(
  x,
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
)

length(intersect(Vedo_NR_SR.filter$GeneName, intersect(Vedo_NR_R.filter$GeneName, Vedo_NR_PR.filter$GeneName)))

x.df <- as.data.frame(list_to_matrix(x))
x.df$colsum <- rowSums(x.df)
x.df$GeneID <- rownames(x.df)
x.df <- x.df %>% left_join(gene.annotation, by = "GeneID")

write.table(
  x.df,
  file.path(de_results_dir, paste0("DEG_Intersection_", treatment, "_", visit, ".txt")),
  sep = ";",
  col.names = TRUE,
  row.names = FALSE
)

DE_all <- x.df %>% dplyr::filter(colsum == 3)

DE_all <- DE_all %>%
  dplyr::left_join(Vedo_NR_SR.filter[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "GeneID")], by = "GeneID")

organism <- "hsapiens"
mart_dataset <- "hsapiens_gene_ensembl"
query_id <- "ensembl_gene_id"

n <- 2000
cutoff <- 0.05
DE_results <- DE_all
enrichment_results <- functional_enrichment_dplot(DE_results, organism, mart_dataset, cutoff = 0.05, return_plots = TRUE, activity = "both", n = 40)
enrichment_results[[2]]

# ------------------------------------------------------------------
# Saving results for IPA
# ------------------------------------------------------------------
cols_need <- c("GeneID", "GeneName")
DE_total_complete <- Vedo_NR_PR %>% dplyr::rename(
  Vedo_NR_PR.log2fc = log2FoldChange,
  Vedo_NR_PR.pvalue = pvalue,
  Vedo_NR_PR.padj = padj
)

cols_stats <- c("Vedo_NR_PR.log2fc", "Vedo_NR_PR.pvalue", "Vedo_NR_PR.padj")
DE_total_complete <- DE_total_complete %>% dplyr::select(c(cols_need, cols_stats))

DE_total_complete <- DE_total_complete %>%
  left_join(Vedo_NR_R) %>% dplyr::rename(
    Vedo_NR_R.log2fc = log2FoldChange,
    Vedo_NR_R.pvalue = pvalue,
    Vedo_NR_R.padj = padj
  )

cols4_stats <- c("Vedo_NR_R.log2fc", "Vedo_NR_R.pvalue", "Vedo_NR_R.padj")
DE_total_complete <- DE_total_complete %>% dplyr::select(c(cols_need, cols_stats, cols4_stats))

DE_total_complete <- DE_total_complete %>%
  left_join(Vedo_NR_SR) %>% dplyr::rename(
    Vedo_NR_SR.log2fc = log2FoldChange,
    Vedo_NR_SR.pvalue = pvalue,
    Vedo_NR_SR.padj = padj
  )

cols8_stats <- c("Vedo_NR_SR.log2fc", "Vedo_NR_SR.pvalue", "Vedo_NR_SR.padj")
DE_total_complete <- DE_total_complete %>% dplyr::select(c(cols_need, cols_stats, cols4_stats, cols8_stats))

DE_total_complete[is.na(DE_total_complete$Vedo_NR_PR.log2fc), "Vedo_NR_PR.log2fc"] <- 0
DE_total_complete[is.na(DE_total_complete$Vedo_NR_R.log2fc), "Vedo_NR_R.log2fc"] <- 0
DE_total_complete[is.na(DE_total_complete$Vedo_NR_SR.log2fc), "Vedo_NR_SR.log2fc"] <- 0

DE_total_complete[is.na(DE_total_complete$Vedo_NR_PR.pvalue), "Vedo_NR_PR.pvalue"] <- 1
DE_total_complete[is.na(DE_total_complete$Vedo_NR_R.pvalue), "Vedo_NR_R.pvalue"] <- 1
DE_total_complete[is.na(DE_total_complete$Vedo_NR_SR.pvalue), "Vedo_NR_SR.pvalue"] <- 1

DE_total_complete[is.na(DE_total_complete$Vedo_NR_PR.padj), "Vedo_NR_PR.padj"] <- 1
DE_total_complete[is.na(DE_total_complete$Vedo_NR_R.padj), "Vedo_NR_R.padj"] <- 1
DE_total_complete[is.na(DE_total_complete$Vedo_NR_SR.padj), "Vedo_NR_SR.padj"] <- 1

write.table(
  DE_total_complete,
  file = file.path(de_results_dir, "IPA_DEG_Vedolizumab_BL.txt"),
  sep = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE
)

# ------------------------------------------------------------------
# LRT DE analysis
# ------------------------------------------------------------------
library(DEGreport)

treatment <- "Vedolizumab"
visit <- "Visit 1"

sample.metadata.vedo.BL <- sample.metadata.filtered %>% dplyr::filter(Treatment == treatment, Visit_Type %in% visit)
sample.metadata.vedo.BL <- sample.metadata.vedo.BL %>% dplyr::left_join(vedo_trajectory_jie, by = c("SubjectID" = "SUBJID"))

sample.metadata.vedo.BL$trajectory <- factor(sample.metadata.vedo.BL$trajectory)

varsity.FPKM_filter <- varsity.FPKM[, sample.metadata.vedo.BL$Originating_ID]
varsity.count_filter <- varsity.count[, sample.metadata.vedo.BL$Originating_ID]
pseudoCount_filter <- log2(varsity.FPKM_filter + 1)

dds <- DESeqDataSetFromMatrix(
  countData = round(varsity.count_filter),
  colData = sample.metadata.vedo.BL,
  design = ~trajectory
)

dds3_m <- dds

expression.count.annot <- gene.annotation[rownames(varsity.count_filter), ]
mcols(dds3_m) <- expression.count.annot
dds3_m <- DESeq(dds3_m, test = "LRT", reduced = ~1)

DEGs_df <- results(dds3_m, independentFiltering = FALSE, alpha = 0.05)
summary(DEGs_df)

sig_res_LRT <- DEGs_df %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  dplyr::filter(padj < 0.08) %>% dplyr::arrange(padj)

sig_res_LRT <- sig_res_LRT %>% left_join(gene.annotation, by = c("gene" = "GeneName"))
n <- 2000
n <- ifelse(n > nrow(sig_res_LRT), nrow(sig_res_LRT), n)
clustering_sig_genes <- sig_res_LRT %>%
  dplyr::arrange(padj) %>%
  head(n)

cluster_pseudo <- pseudoCount_filter[clustering_sig_genes$GeneID, ]
metadf <- sample.metadata.vedo.BL
rownames(metadf) <- metadf$Originating_ID

cluster_pseudo <- na.omit(cluster_pseudo)
metadf <- na.omit(metadf)

variance <- apply(cluster_pseudo, 1, var)
cluster_pseudo <- cluster_pseudo[variance > 0, ]

metadf$Visit_Type <- factor(metadf$Visit_Type, levels = c("Visit 1", "Visit 9", "Visit 28"))
clusters <- degPatterns(cluster_pseudo, metadata = metadf, time = "trajectory", plot = TRUE)

head(clusters$df)
clusters$plot

clusters$df <- clusters$df %>% dplyr::left_join(sig_res_LRT, by = c("genes" = "GeneID"))
View(clusters$df)

write.table(
  clusters$df,
  file.path(de_results_dir, paste0("DEG_LRT_", treatment, "_", visit, "_NR-SR_FDR005.txt")),
  sep = ";",
  col.names = TRUE,
  row.names = FALSE
)

DE_results <- clusters$df %>% dplyr::filter(cluster == 2)
DE_results$GeneName <- DE_results$gene
ensembl <- useEnsembl(biomart = "genes")

organism <- "hsapiens"
mart_dataset <- "hsapiens_gene_ensembl"
query_id <- "ensembl_gene_id"

n <- 2000
cutoff <- 0.05
DE_results <- DE_results %>% dplyr::arrange(padj) %>% head(n)
enrichment_results <- functional_enrichment_dplot(DE_results, organism, mart_dataset, cutoff = 0.05, return_plots = TRUE, activity = "both", n = 40)
enrichment_results[[2]]

# ------------------------------------------------------------------
# Functional analysis using MSigDB
# ------------------------------------------------------------------
treatment <- "Vedolizumab"
timepoint <- "Visit 9"
control <- "super_responder"

DEG_results <- read.delim(
  file.path(de_results_dir, glue("{treatment}_{timepoint}_non_responder_vs_{control}_DEG_results_DEGall_New_EndoW14NR.csv")),
  sep = ","
)

DE <- DEG_results %>% dplyr::arrange(desc(log2FoldChange))
gene_list <- DE$log2FoldChange
names(gene_list) <- DE$GeneName
gene_list <- sort(gene_list, decreasing = TRUE)

msigdbr_species <- "Homo sapiens"
msigdbr_category <- "H"

msigdbr_gene_sets <- msigdbr(species = msigdbr_species, category = msigdbr_category)
msigdbr_list <- split(msigdbr_gene_sets$gene_symbol, msigdbr_gene_sets$gs_name)

fgsea_results <- fgsea(pathways = msigdbr_list, stats = gene_list, nperm = 1000)
fgsea_results$leadingEdge <- NULL
fgsea_results$pathway <- gsub("^HALLMARK_", "", fgsea_results$pathway)

fgsea_df <- fgsea_results %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.06) %>%
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

CairoPDF(file = file.path(manuscript_eps_dir, glue("FIG2B_{treatment}_{timepoint}_{control}_Volcano_EndoNR_W14.pdf")), width = 10, height = 10)
plot(dotplot)
dev.off()

cairo_ps(file.path(manuscript_eps_dir, glue("FIG2B_{treatment}_{timepoint}_{control}_Volcano_EndoNR_W14.eps")), width = 10, height = 10)
print(dotplot)
dev.off()

# ------------------------------------------------------------------
# Plotting boxplots
# ------------------------------------------------------------------
meta.t <- sample.metadata.vedo.BL
FPKM.t <- varsity.FPKM.vedo.BL
FPKM.t$GeneID <- rownames(FPKM.t)
FPKM.t$GeneID <- as.character(FPKM.t$GeneID)

FPKM.t <- FPKM.t %>%
  dplyr::left_join(gene.annotation[, c("GeneID", "GeneName")], by = c("GeneID" = "GeneID")) %>%
  dplyr::filter(!duplicated(GeneName), GeneName %in% c("C10orf99", "IL12A", "IL1A", "IL1B", "IL6", "IL23A", "TNF"))

rownames(FPKM.t) <- FPKM.t$GeneName
FPKM.t$GeneID <- FPKM.t$GeneName <- NULL
FPKM.t <- as.data.frame(t(FPKM.t))
FPKM.t$Sample <- rownames(FPKM.t)
FPKM.t <- FPKM.t %>% dplyr::left_join(meta.t[, c("Originating_ID", "trajectory")], by = c("Sample" = "Originating_ID"))
str(FPKM.t)

FPKM.t <- FPKM.t %>%
  mutate(across(IL1A:IL23A, ~log2(as.numeric(trimws(.)) + 1)))

data_long <- FPKM.t %>%
  pivot_longer(cols = IL1A:IL23A, names_to = "gene", values_to = "FPKM")

ggplot(data_long, aes(x = trajectory, y = FPKM, fill = trajectory)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~gene, scales = "free_y") +
  labs(
    title = "Log2(FPKM+1) Values by Trajectory Group for Each Gene - W14 Tissue Samples - Saman's labels",
    x = "Trajectory Group",
    y = "Log2(FPKM+1)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(
    fun.data = function(y) {
      data.frame(y = quantile(y, probs = c(0.25, 0.5, 0.75)), label = c("", "", ""))
    },
    geom = "text", position = position_dodge(width = 0.75), vjust = -0.5
  )


# ------------------------------------------------------------------
# Check RNAseq sample numbers for RNR labels based on clinical, mucosal and disease control
# ------------------------------------------------------------------
View(sample.metadata)
View(sample.metadata.filtered)
View(clinical_data)

control <- read_excel(disease_control_file)
colnames(control) <- c("SUBJID", "diseasecontrol")

tt <- sample.metadata.filtered %>%
  dplyr::left_join(sample.metadata[, c("SubjectID", "Geboes_W52", "Histo_W52", "Endo_W52")], by = "SubjectID") %>%
  distinct()
tt <- tt %>% dplyr::left_join(vedo_trajectory_jie, by = c("SubjectID" = "SUBJID"))
tt <- tt %>% dplyr::left_join(control, by = c("SubjectID" = "SUBJID"))

filtered_tt <- tt %>%
  dplyr::filter(Visit_Type %in% c("Visit 1", "Visit 9"))

complete_subjects <- filtered_tt %>%
  group_by(SubjectID) %>%
  summarize(num_visits = n_distinct(Visit_Type)) %>%
  dplyr::filter(num_visits == 2)

complete_subjects <- complete_subjects %>% dplyr::left_join(sample.metadata.filtered[, c("SubjectID", "Treatment")], by = c("SubjectID"))
complete_subjects <- complete_subjects %>% distinct()
table(complete_subjects$Treatment)

filtered_tt <- tt %>% dplyr::filter(Treatment == "Vedolizumab", has_RNAseq == "Yes")
filtered_tt <- filtered_tt %>% dplyr::select(SubjectID, Remission_W52, trajectory) %>% distinct()
table(filtered_tt$trajectory, filtered_tt$Remission_W52)

admayo <- read.delim(admayo_file, sep = ",")
admayo.filt <- admayo[apply(admayo, 1, function(row) any(grepl("Subscore", row))), ]

admayo.filt <- admayo.filt %>% dplyr::select(c("USUBJID", "PARAM", "AVAL", "AVISIT"))

admayo.wide <- pivot_wider(
  data = admayo.filt,
  id_cols = c(USUBJID, AVISIT),
  names_from = PARAM,
  values_from = AVAL
)
colnames(admayo.wide) <- c("USUBJID", "AVISIT", "PGS", "ES", "RB", "SF")
admayo.wide$USUBJID <- sub("^MLN0002-3026-", "", admayo.wide$USUBJID)
admayo.wide <- admayo.wide %>% dplyr::filter(AVISIT %in% c("Week 14", "Week 52")) %>% dplyr::select(-PGS, -RB, -SF)
str(admayo.wide)

admayo.wide <- admayo.wide %>%
  mutate(
    ES_Status = case_when(
      is.na(ES) ~ NA_character_,
      ES <= 1 ~ "Yes",
      ES > 1 ~ "No"
    )
  )

admayo.wide.final <- admayo.wide %>%
  dplyr::select(USUBJID, AVISIT, ES_Status) %>%
  pivot_wider(
    names_from = AVISIT,
    values_from = ES_Status,
    values_fill = list(ES_Status = NA_character_)
  )

colnames(admayo.wide.final) <- make.names(colnames(admayo.wide.final))
admayo.wide.final <- admayo.wide.final %>% dplyr::rename("EndoResW14" = Week.14, "EndoResW52" = Week.52)

filtered_tt <- filtered_tt %>% dplyr::left_join(admayo.wide.final, by = c("SubjectID" = "USUBJID")) %>% distinct()
table(filtered_tt$trajectory, filtered_tt$EndoResW52)

filtered_tt <- filtered_tt %>% dplyr::left_join(control, by = c("SubjectID" = "SUBJID")) %>% distinct()
table(filtered_tt$trajectory, filtered_tt$diseasecontrol)

# ------------------------------------------------------------------
# Categorizing the patients into different trajectory groups based on their progression
# ------------------------------------------------------------------
clinical_data <- readRDS(clinical_data_file)
vedo_trajectory_jie <- read.delim(file.path(trajectory_analysis_dir, "Jie's labels.csv"), sep = ",", header = TRUE)

admayo <- read.delim(admayo_file, sep = ",")
admayo.filt <- admayo[apply(admayo, 1, function(row) any(grepl("Subscore", row))), ]

admayo.filt <- admayo.filt %>% dplyr::select(c("USUBJID", "PARAM", "AVAL", "AVISIT"))

admayo.wide <- pivot_wider(
  data = admayo.filt,
  id_cols = c(USUBJID, AVISIT),
  names_from = PARAM,
  values_from = AVAL
)
colnames(admayo.wide) <- c("USUBJID", "AVISIT", "PGS", "ES", "RB", "SF")
admayo.wide$USUBJID <- sub("^MLN0002-3026-", "", admayo.wide$USUBJID)
admayo.wide <- admayo.wide %>% dplyr::filter(AVISIT %in% c("Week 14", "Week 52")) %>% dplyr::select(-PGS, -RB, -SF)
str(admayo.wide)

admayo.wide <- admayo.wide %>%
  mutate(
    ES_Status = case_when(
      is.na(ES) ~ NA_character_,
      ES <= 1 ~ "Yes",
      ES > 1 ~ "No"
    )
  )

admayo.wide.final <- admayo.wide %>%
  dplyr::select(USUBJID, AVISIT, ES_Status) %>%
  pivot_wider(
    names_from = AVISIT,
    values_from = ES_Status,
    values_fill = list(ES_Status = NA_character_)
  )

colnames(admayo.wide.final) <- make.names(colnames(admayo.wide.final))
admayo.wide.final <- admayo.wide.final %>% dplyr::rename("EndoResW14" = Week.14, "EndoResW52" = Week.52)

vedo_trajectory_jie <- vedo_trajectory_jie %>% dplyr::left_join(admayo.wide.final, by = c("SUBJID" = "USUBJID"))
vedo_trajectory_jie <- vedo_trajectory_jie %>% dplyr::left_join(clinical_data[, c("SubjectID", "Treatment", "has_RNAseq", "Response_W14", "Remission_W52")], by = c("SUBJID" = "SubjectID"))
vedo_trajectory_jie <- vedo_trajectory_jie %>% distinct()

tt <- vedo_trajectory_jie %>% dplyr::filter(Treatment == "Vedolizumab")
table(tt$trajectory, tt$EndoResW14, useNA = "ifany")
table(tt$EndoResW14, tt$EndoResW52, tt$trajectory, useNA = "ifany")

tt <- vedo_trajectory_jie %>% dplyr::filter(Treatment == "Adalimumab")
table(tt$trajectory, tt$EndoResW14, useNA = "ifany")
table(tt$EndoResW14, tt$EndoResW52, tt$trajectory, useNA = "ifany")

tt <- vedo_trajectory_jie %>% dplyr::filter(Treatment == "Vedolizumab", has_RNAseq == "Yes")
table(tt$trajectory, tt$EndoResW14, useNA = "ifany")
table(tt$EndoResW14, tt$EndoResW52, tt$trajectory, useNA = "ifany")

tt <- vedo_trajectory_jie %>% dplyr::filter(Treatment == "Adalimumab", has_RNAseq == "Yes")
table(tt$trajectory, tt$EndoResW14, useNA = "ifany")
table(tt$EndoResW14, tt$EndoResW52, tt$trajectory, useNA = "ifany")

# ------------------------------------------------------------------
# Perform trajectory DE analysis for the previous results showing Endo response at W14 changed at W52
# ------------------------------------------------------------------
selected_subjects <- vedo_trajectory_jie %>% dplyr::filter(Treatment == "Adalimumab", EndoResW14 == "No", has_RNAseq == "Yes", trajectory %in% c("non_responder", "partial_responder", "responder", "super_responder"))
table(selected_subjects$EndoResW14, selected_subjects$trajectory)

selected_subjects <- selected_subjects %>%
  dplyr::left_join(sample.metadata.filtered[, c("SubjectID", "Originating_ID", "Visit_Type", "Test_Name_or_Sample_Type")], by = c("SUBJID" = "SubjectID")) %>%
  distinct()

sample.metadata <- selected_subjects %>% dplyr::filter(Visit_Type == "Visit 9")
table(sample.metadata$trajectory)

varsity.count.vedo.BL <- varsity.count[, sample.metadata$Originating_ID]
varsity.FPKM.vedo.BL <- varsity.FPKM[, sample.metadata$Originating_ID]
pseudoCount <- log2(varsity.FPKM.vedo.BL + 1)

comparison <- c("non_responder", "super_responder")

t_pseudoCount <- pseudoCount[, as.character(sample.metadata$Originating_ID)]
t_count <- varsity.count.vedo.BL[, as.character(sample.metadata$Originating_ID)]

dds <- DESeqDataSetFromMatrix(
  countData = round(t_count),
  colData = sample.metadata,
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

pos <- nrow(DEG_all[(DEG_all$log2FoldChange > 0 & DEG_all$pvalue < 0.05), ])
neg <- nrow(DEG_all[(DEG_all$log2FoldChange < 0 & DEG_all$pvalue < 0.05), ])
print(paste0("positive adj pvalue=", pos))
print(paste0("positive adj pvalue=", neg))

DEG.comparison_filtered <- DEG_all %>% dplyr::filter(padj < 0.05)

pos <- nrow(DEG.comparison_filtered[DEG.comparison_filtered$log2FoldChange > 0, ])
neg <- nrow(DEG.comparison_filtered[DEG.comparison_filtered$log2FoldChange < 0, ])

print(paste0("positive adj pvalue=", pos))
print(paste0("positive adj pvalue=", neg))
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
  selectLab = gene.annotation[DEG.table$GeneName[1:15], "GeneName"],
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

DE <- DEG_all %>% dplyr::arrange(desc(log2FoldChange))
gene_list <- DE$log2FoldChange
names(gene_list) <- DE$GeneName
gene_list <- sort(gene_list, decreasing = TRUE)
msigdbr_species <- "Homo sapiens"
msigdbr_category <- "H"
msigdbr_gene_sets <- msigdbr(species = msigdbr_species, category = msigdbr_category)
msigdbr_list <- split(msigdbr_gene_sets$gene_symbol, msigdbr_gene_sets$gs_name)
fgsea_results <- fgsea(pathways = msigdbr_list, stats = gene_list, nperm = 1000)

fgsea_df <- fgsea_results %>%
  as.data.frame() %>% dplyr::select(-leadingEdge)

fgsea_df <- fgsea_results %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::mutate(pathway = reorder(pathway, NES))

res <- fgsea_df %>% dplyr::select(-leadingEdge)

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

# ------------------------------------------------------------------
# Running the DE analysis in loop for all trajectory groups vs NR for W14EndoNR
# ------------------------------------------------------------------
visit_t <- "Visit 9"
treatmentt <- "Adalimumab"

selected_subjects <- vedo_trajectory_jie %>% dplyr::filter(Treatment == treatmentt, EndoResW14 == "No", has_RNAseq == "Yes", trajectory %in% c("non_responder", "partial_responder", "responder", "super_responder"))
table(selected_subjects$EndoResW14, selected_subjects$trajectory)

selected_subjects <- selected_subjects %>%
  dplyr::left_join(sample.metadata.filtered[, c("SubjectID", "Originating_ID", "Visit_Type", "Test_Name_or_Sample_Type")], by = c("SUBJID" = "SubjectID")) %>%
  distinct()

comparisons <- list(
  "non_responder_vs_partial_responder" = c("non_responder", "partial_responder"),
  "non_responder_vs_responder" = c("non_responder", "responder"),
  "non_responder_vs_super_responder" = c("non_responder", "super_responder")
)

comp_name <- "non_responder_vs_super_responder"

for (comp_name in names(comparisons)) {
  comparison <- comparisons[[comp_name]]

  sample_metadata <- selected_subjects %>%
    dplyr::filter(trajectory %in% comparison, Visit_Type == visit_t) %>%
    dplyr::distinct()

  t_count <- varsity.count[, as.character(sample_metadata$Originating_ID)]
  t_pseudoCount <- log2(varsity.FPKM[, as.character(sample_metadata$Originating_ID)] + 1)

  dds <- DESeqDataSetFromMatrix(
    countData = round(t_count),
    colData = sample_metadata,
    design = ~trajectory
  )
  dds$label <- factor(dds$trajectory, levels = comparison)

  DESeq_dd <- DESeq(dds)
  DEG.comparison <- results(DESeq_dd, contrast = c("trajectory", comparison), alpha = 0.05)
  DEG_all <- as.data.frame(DEG.comparison)

  DEG_all$GeneID <- rownames(DEG_all)
  DEG_all <- DEG_all %>% left_join(gene.annotation, by = "GeneID")
  DEG_all$GeneName <- toupper(DEG_all$GeneName)
  DEG_all <- DEG_all %>% arrange(padj)

  up_raw <- nrow(DEG_all %>% dplyr::filter(log2FoldChange > 0, pvalue < 0.05))
  down_raw <- nrow(DEG_all %>% dplyr::filter(log2FoldChange < 0, pvalue < 0.05))
  up_adj <- nrow(DEG_all %>% dplyr::filter(log2FoldChange > 0, padj < 0.05))
  down_adj <- nrow(DEG_all %>% dplyr::filter(log2FoldChange < 0, padj < 0.05))

  annots <- DEG_all$GeneName
  volcano_plot <- EnhancedVolcano(
    DEG_all,
    lab = annots,
    x = "log2FoldChange",
    y = "padj",
    selectLab = gene.annotation[DEG_all$GeneName[1:15], "GeneName"],
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
  ) +
    annotate(
      "text", x = Inf, y = Inf,
      label = sprintf("Up (raw): %d\nDown (raw): %d\nUp (adj): %d\nDown (adj): %d", up_raw, down_raw, up_adj, down_adj),
      hjust = 1.1, vjust = 1.1, size = 3.5, color = "black"
    )

  DE <- DEG_all %>% dplyr::arrange(desc(log2FoldChange))
  gene_list <- DE$log2FoldChange
  names(gene_list) <- DE$GeneName
  gene_list <- sort(gene_list, decreasing = TRUE)
  msigdbr_species <- "Homo sapiens"
  msigdbr_category <- "H"
  msigdbr_gene_sets <- msigdbr(species = msigdbr_species, category = msigdbr_category)
  msigdbr_list <- split(msigdbr_gene_sets$gene_symbol, msigdbr_gene_sets$gs_name)
  fgsea_results <- fgsea(pathways = msigdbr_list, stats = gene_list, nperm = 1000)

  fgsea_df <- fgsea_results %>%
    as.data.frame() %>% dplyr::select(-leadingEdge)

  write.table(fgsea_df, file.path(fe_results_dir, paste0(treatmentt, "_", visit_t, "_", comp_name, "_fgsea_results_New.csv")), sep = ",")

  fgsea_df <- fgsea_results %>%
    as.data.frame() %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::mutate(pathway = reorder(pathway, NES))

  res <- fgsea_df %>% dplyr::select(-leadingEdge)

  write.csv(res, file.path(fe_results_dir, paste0(treatmentt, "_", visit_t, "_", comp_name, "_fgsea_results_New_filtered.csv")))

  dotplot <- ggplot(fgsea_df, aes(x = NES, y = pathway)) +
    geom_point(aes(size = -log10(padj), color = padj)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(
      title = paste0(comp_name, "-  FDR < 0.05"),
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway",
      size = "-log10(Adjusted P-value)",
      color = "Adjusted P-value"
    ) +
    theme_minimal()

  ggsave(file.path(figs_dir, paste0(treatmentt, "_", visit_t, "_", comp_name, "_dot_plot.png")), plot = dotplot, width = 9, height = 6)
}

print("Analysis and plotting completed for all comparisons.")
