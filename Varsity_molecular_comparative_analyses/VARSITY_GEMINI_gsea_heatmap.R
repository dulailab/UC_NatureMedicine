# VARSITY_GEMINI_gsea_heatmap.R
# Objective: assemble Hallmark fgsea results across Vedo/Ada timepoints and GEMINI
# contrasts, then render NES heatmaps with significance masking
# with configurable paths at the top.

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
project_dir <- getwd()
fig2_dir <- file.path(project_dir, "Fig2")
fe_results_dir <- file.path(project_dir, "FE_results")
fe_results_varsity_dir <- file.path(fe_results_dir, "VARSITY_GSEA_Final")
fe_results_gemini_dir <- file.path(fe_results_dir, "GEMINI_GSEA_Final")
figs_dir <- file.path(project_dir, "Figs")
manuscript_eps_dir <- file.path(figs_dir, "Manuscript_EPS")

dir.create(fe_results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_eps_dir, recursive = TRUE, showWarnings = FALSE) 

# ------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(purrr)
library(grid)
library(glue)

# ------------------------------------------------------------------
# Helper
# ------------------------------------------------------------------
read_and_process_data <- function(file_path, prefix) {
  data <- read.delim(file_path, sep = ",")
  data <- data %>% dplyr::select(pathway, padj, NES)
  colnames(data)[2:3] <- paste(prefix, colnames(data)[2:3], sep = "_")
  data
}

crossed_pathways <- c(
  "SPERMATOGENESIS",
  "PI3K_AKT_MTOR_SIGNALING",
  "APICAL_JUNCTION",
  "MYC_TARGETS_V2",
  "DNA_REPAIR",
  "NOTCH_SIGNALING",
  "ESTROGEN_RESPONSE_LATE",
  "MITOTIC_SPINDLE",
  "PEROXISOME",
  "HEDGEHOG_SIGNALING",
  "APICAL_SURFACE",
  "P53_PATHWAY",
  "WNT_BETA_CATENIN_SIGNALING",
  "TGF_BETA_SIGNALING"
)

# ------------------------------------------------------------------
# Vedo/Ada x timepoint heatmap
# ------------------------------------------------------------------
vedo_nr_pr_W0 <- read_and_process_data(file.path(fe_results_varsity_dir, "Vedolizumab_Visit 1_non_responder_vs_partial_responder_fgsea_results_New.csv"), "Vedo_NRvPR_W0")
vedo_nr_r_W0 <- read_and_process_data(file.path(fe_results_varsity_dir, "Vedolizumab_Visit 1_non_responder_vs_responder_fgsea_results_New.csv"), "Vedo_NRvR_W0")
vedo_nr_sr_W0 <- read_and_process_data(file.path(fe_results_varsity_dir, "Vedolizumab_Visit 1_non_responder_vs_super_responder_fgsea_results_New.csv"), "Vedo_NRvSR_W0")

vedo_nr_pr_W14 <- read_and_process_data(file.path(fe_results_varsity_dir, "Vedolizumab_Visit 9_non_responder_vs_partial_responder_fgsea_results_New.csv"), "Vedo_NRvPR_W14")
vedo_nr_r_W14 <- read_and_process_data(file.path(fe_results_varsity_dir, "Vedolizumab_Visit 9_non_responder_vs_responder_fgsea_results_New.csv"), "Vedo_NRvR_W14")
vedo_nr_sr_W14 <- read_and_process_data(file.path(fe_results_varsity_dir, "Vedolizumab_Visit 9_non_responder_vs_super_responder_fgsea_results_New.csv"), "Vedo_NRvSR_W14")

ada_nr_pr_W0 <- read_and_process_data(file.path(fe_results_varsity_dir, "Adalimumab_Visit 1_non_responder_vs_partial_responder_fgsea_results_New.csv"), "Ada_NRvPR_W0")
ada_nr_r_W0 <- read_and_process_data(file.path(fe_results_varsity_dir, "Adalimumab_Visit 1_non_responder_vs_responder_fgsea_results_New.csv"), "Ada_NRvR_W0")
ada_nr_sr_W0 <- read_and_process_data(file.path(fe_results_varsity_dir, "Adalimumab_Visit 1_non_responder_vs_super_responder_fgsea_results_New.csv"), "Ada_NRvSR_W0")

ada_nr_pr_W14 <- read_and_process_data(file.path(fe_results_varsity_dir, "Adalimumab_Visit 9_non_responder_vs_partial_responder_fgsea_results_New.csv"), "Ada_NRvPR_W14")
ada_nr_r_W14 <- read_and_process_data(file.path(fe_results_varsity_dir, "Adalimumab_Visit 9_non_responder_vs_responder_fgsea_results_New.csv"), "Ada_NRvR_W14")
ada_nr_sr_W14 <- read_and_process_data(file.path(fe_results_varsity_dir, "Adalimumab_Visit 9_non_responder_vs_super_responder_fgsea_results_New.csv"), "Ada_NRvSR_W14")

pp_df <- purrr::reduce(
  list(
    vedo_nr_pr_W0, vedo_nr_r_W0, vedo_nr_sr_W0,
    vedo_nr_pr_W14, vedo_nr_r_W14, vedo_nr_sr_W14,
    ada_nr_pr_W0, ada_nr_r_W0, ada_nr_sr_W0,
    ada_nr_pr_W14, ada_nr_r_W14, ada_nr_sr_W14
  ),
  full_join, by = "pathway"
)

pp_df$pathway <- gsub("HALLMARK_", "", pp_df$pathway)
pp_df <- pp_df %>% dplyr::filter(!pathway %in% crossed_pathways)

nes_matrix <- as.matrix(pp_df[, grep("NES", names(pp_df))])
padj_matrix <- as.matrix(pp_df[, grep("padj", names(pp_df))])

rownames(nes_matrix) <- pp_df$pathway
colnames(nes_matrix) <- names(pp_df)[grep("NES", names(pp_df))]
rownames(padj_matrix) <- pp_df$pathway
colnames(padj_matrix) <- names(pp_df)[grep("padj", names(pp_df))]

significance_matrix <- apply(padj_matrix, c(1, 2), function(x) ifelse(x < 0.05, "*", ""))

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun_grey <- colorRamp2(c(-2, 0, 2), c("white", "white", "white"))

p2 <- Heatmap(
  nes_matrix,
  name = "NES",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (padj_matrix[i, j] < 0.05) {
      grid.text(sprintf("%.2f%s", nes_matrix[i, j], significance_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
    } else {
      grid.rect(x = x, y = y, width = width, height = height, gp = gpar(fill = "grey93", col = "grey93"))
    }
  },
  heatmap_legend_param = list(title = "NES", at = c(-2, 0, 2), labels = c("-2", "0", "2"))
)

hm_drawn <- draw(p2)
row_ord <- row_order(hm_drawn)
nes_matrix_varsity <- nes_matrix

# ------------------------------------------------------------------
# GEMINI study integrated heatmap
# ------------------------------------------------------------------
GEMINI_NR_SR_W6 <- read.delim(file.path(fe_results_gemini_dir, "fgsea_results_NRvsSR_W6.csv"), sep = ",") %>% dplyr::select(pathway, padj, NES)
GEMINI_NR_R_W6 <- read.delim(file.path(fe_results_gemini_dir, "fgsea_results_NRvsR_W6.txt"), sep = ";") %>% dplyr::select(pathway, padj, NES) %>% dplyr::filter(pathway %in% GEMINI_NR_SR_W6$pathway)
GEMINI_NR_PR_W6 <- read.delim(file.path(fe_results_gemini_dir, "fgsea_results_NRvsPR_W6.txt"), sep = ";") %>% dplyr::select(pathway, padj, NES) %>% dplyr::filter(pathway %in% GEMINI_NR_SR_W6$pathway)
GEMINI_NR_PR_W0 <- read.delim(file.path(fe_results_gemini_dir, "fgsea_results_NRvsPR_W0.txt"), sep = ";") %>% dplyr::select(pathway, padj, NES) %>% dplyr::filter(pathway %in% GEMINI_NR_SR_W6$pathway)
GEMINI_NR_R_W0 <- read.delim(file.path(fe_results_gemini_dir, "fgsea_results_NRvsR_W0.txt"), sep = ";") %>% dplyr::select(pathway, padj, NES) %>% dplyr::filter(pathway %in% GEMINI_NR_SR_W6$pathway)
GEMINI_NR_SR_W0 <- read.delim(file.path(fe_results_gemini_dir, "fgsea_results_NRvsSR_W0.txt"), sep = ";") %>% dplyr::select(pathway, padj, NES) %>% dplyr::filter(pathway %in% GEMINI_NR_SR_W6$pathway)

GEMINI_NR_SR_W6$pathway <- GEMINI_NR_PR_W6$pathway <- GEMINI_NR_PR_W0$pathway <- GEMINI_NR_R_W0$pathway <- GEMINI_NR_SR_W0$pathway <- NULL
pp_df <- cbind(GEMINI_NR_PR_W0, GEMINI_NR_R_W0, GEMINI_NR_SR_W0, GEMINI_NR_PR_W6, GEMINI_NR_R_W6, GEMINI_NR_SR_W6)
pp_df$pathway <- gsub("HALLMARK_", "", pp_df$pathway)

colnames(pp_df) <- make.unique(names(pp_df))
pp_df <- pp_df %>% dplyr::filter(pathway %in% rownames(nes_matrix_varsity))

nes_matrix <- as.matrix(pp_df[, grep("NES", names(pp_df))])
padj_matrix <- as.matrix(pp_df[, grep("padj", names(pp_df))])

nn <- c("GEMINI_NR_PR_W0", "GEMINI_NR_R_W0", "GEMINI_NR_SR_W0", "GEMINI_NR_PR_W6", "GEMINI_NR_R_W6", "GEMINI_NR_SR_W6")
rownames(nes_matrix) <- pp_df$pathway
colnames(nes_matrix) <- nn
rownames(padj_matrix) <- pp_df$pathway
colnames(padj_matrix) <- nn

significance_matrix <- apply(padj_matrix, c(1, 2), function(x) ifelse(x < 0.05, "*", ""))

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

nes_matrix <- nes_matrix[row_ord, ]
significance_matrix <- significance_matrix[row_ord, ]

ttt <- Heatmap(
  nes_matrix,
  name = "NES",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f%s", nes_matrix[i, j], significance_matrix[i, j]), x, y,
      gp = gpar(fontsize = 10)
    )
  },
  heatmap_legend_param = list(
    title = "NES",
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2")
  ),
  show_heatmap_legend = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10)
)

print(ttt)

# CairoPDF(file = file.path(manuscript_eps_dir, "GEMINI_PR_SR_Heatmap_new.pdf"), width = 5, height = 10)
# print(ttt)
# dev.off()
# cairo_ps(file.path(manuscript_eps_dir, "GEMINI_PR_SR_Heatmap_new.eps"), width = 5, height = 10)
# print(ttt)
# dev.off()
