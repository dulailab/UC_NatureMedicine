# Fig2e_analysis.R
# Objective: map UNIFI contrasts to Hallmark/trajectory signatures, run fgsea,
# and combine with Varsity/GEMINI enrichment for dot/heatmap summaries.
# Script version of Fig2e_analysis.Rmd. Logic/order is unchanged; paths are
# configurable at the top to simplify reuse.

# correspondong to Figure 2E

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
project_dir <- if (basename(getwd()) == "Fig2") normalizePath("..", mustWork = FALSE) else getwd()
fig2_dir <- file.path(project_dir, "Fig2")
files_dir <- file.path(project_dir, "files")
unifi_dir <- file.path(files_dir, "UNIFI_data")
fig2_files_dir <- file.path(fig2_dir, "files")
fe_results_dir <- file.path(project_dir, "GSEA_results")

unifi_comparisons_file <- file.path(unifi_dir, "selected_comparisons.xlsx")
unifi_sample_meta_file <- file.path(unifi_dir, "GSE206285_additional_sample_metadata.tsv")
unifi_de_stats_file <- file.path(unifi_dir, "Comparison.FoldChangePValueTable.GSE206285.txt")
unifi_gene_annotation_file <- file.path(unifi_dir, "gene_annotation.xlsx")
hallmark_sigs_file <- file.path(fig2_files_dir, "VDZ_Varsity_trajectory_signatures_all_hallmarks_padj05_V2.xlsx")

unifi_gsea_scores_file <- file.path(files_dir, "RNR_targeted_signatures_gsea", "UNIFI_GSEA_scores.csv")
hallmark_filtered_file <- file.path(files_dir, "RNR_targeted_signatures_gsea", "hallmark_sigs_filtered_new.csv")
fig2e_gsea_out <- file.path(fe_results_dir, "Fig2e_gsea.xlsx")

# ------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------
library(dplyr)
library(readxl)
library(writexl)
library(tidyr)
library(fgsea)
library(stringr)
library(purrr)
library(tibble)

# ------------------------------------------------------------------
# Importing the files
# ------------------------------------------------------------------
comparisons <- read_excel(unifi_comparisons_file) %>%
  dplyr::select(
    c(
      "ComparisonID", "ProjectName", "PlatformName",
      "ComparisonType", "ComparisonCategory", "ComparisonContrast",
      "SampleDataMode", "Case.DiseaseCategory", "Case.DiseaseState"
    )
  )
sample_meta <- read.delim(unifi_sample_meta_file)
DE_stats <- read.delim(unifi_de_stats_file)
gene.annot <- read_excel(unifi_gene_annotation_file)

hallmark_sigs <- read_excel(hallmark_sigs_file) %>%
  dplyr::filter(Name %in% c(
    "ANGIOGENESIS", "EPITHELIAL_MESENCHYMAL_TRANSITION", "IL6_JAK_STAT3_SIGNALING",
    "INFLAMMATORY_RESPONSE", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE",
    "TNFA_SIGNALING_VIA_NFKB", "Neutrophil_degranulation"
  ))

hallmark_sigs_filtered <- hallmark_sigs %>%
  dplyr::select("pathway" = Name, "genes" = NR_up_all)

hallmark_sigs_all <- hallmark_sigs %>%
  dplyr::select("pathway" = Name, "genes" = Gene_List)

# ------------------------------------------------------------------
# Subset DE stats to selected comparisons and relabel columns
# ------------------------------------------------------------------
keep_ids <- comparisons$ComparisonID
keep_cols <- c(
  "GeneID",
  paste0(keep_ids, ".Log2FoldChange"),
  paste0(keep_ids, ".RawPValue"),
  paste0(keep_ids, ".AdjustedPValue")
)

DE_stats_sub <- DE_stats %>%
  dplyr::select(any_of(keep_cols))

contrast_map <- c(
  "ClinicalOutcome:TreatmentGroup => remission at week 8;mucosal healing -> ustekinumab, 130 mg vs placebo" =
    "Remission_w8_MH_Ust130mpk_vs_Pbo",
  "ClinicalOutcome:TreatmentGroup => remission at week 8;mucosal healing -> ustekinumab, 6 mg/kg vs placebo" =
    "Remission_w8_MH_Ust6mpk_vs_Pbo",
  "ClinicalOutcome:TreatmentGroup => no remission at week 8;no mucosal healing -> ustekinumab, 6 mg/kg vs placebo" =
    "NoRemission_w8_NoMH_Ust6mpk_vs_Pbo",
  "ExperimentGroup => no remission at week 8, ustekinumab, 130 mg vs disease control" =
    "NoRemission_w8_Ust130mpk_vs_DC",
  "ExperimentGroup => no remission at week 8, ustekinumab, 6 mg/kg vs disease control" =
    "NoRemission_w8_Ust6mpk_vs_DC",
  "TreatmentGroup:ClinicalOutcome => ustekinumab, 130 mg -> no remission at week 8;no mucosal healing vs remission at week 8;mucosal healing" =
    "NoRemission_w8_NoMH_Ust130mpk_vs_Remission_MH_w8",
  "TreatmentGroup:ClinicalOutcome => ustekinumab, 130 mg -> no remission at week 8 vs remission at week 8" =
    "NoRemission_w8_Ust130mpk_vs_Remission_w8",
  "TreatmentGroup:ClinicalOutcome => ustekinumab, 6 mg/kg -> no remission at week 8;no mucosal healing vs remission at week 8;mucosal healing" =
    "NoRemission_w8_NoMH_Ust6mpk_vs_Remission_MH_w8",
  "ClinicalOutcome:TreatmentGroup => no remission at week 8;no mucosal healing -> ustekinumab, 130 mg vs placebo" =
    "NoRemission_w8_NoMH_Ust130mpk_vs_Pbo"
)

comparisons2 <- comparisons %>%
  mutate(ContrastShort = contrast_map[ComparisonContrast])

id_to_short <- setNames(comparisons2$ContrastShort, comparisons2$ComparisonID)
DE_stats_sub <- DE_stats_sub %>% tibble::column_to_rownames("GeneID")

cols <- colnames(DE_stats_sub)
for (i in seq_along(cols)) {
  for (full_id in names(id_to_short)) {
    pattern <- paste0("^", full_id, "\\.")
    if (grepl(pattern, cols[i])) {
      cols[i] <- sub(pattern, paste0(id_to_short[[full_id]], "."), cols[i])
      break
    }
  }
}
colnames(DE_stats_sub) <- cols

DE_stats_sub <- DE_stats_sub %>%
  tibble::rownames_to_column("GeneID") %>%
  dplyr::left_join(gene.annot[, c("id", "GeneName")], by = c("GeneID" = "id")) %>%
  dplyr::filter(!duplicated(GeneName))

DE_stats_sub <- DE_stats_sub %>%
  tibble::column_to_rownames("GeneName") %>%
  dplyr::select(-GeneID)

str(DE_stats_sub)

# ------------------------------------------------------------------
# Running the GSEA analysis
# ------------------------------------------------------------------
gene_lists <- strsplit(hallmark_sigs_filtered$genes, ",")
pathways <- setNames(gene_lists, hallmark_sigs_filtered$pathway)

cmp_prefixes <- colnames(DE_stats_sub) %>%
  str_remove("\\.(Log2FoldChange|RawPValue|AdjustedPValue)$") %>%
  unique()

de_list <- map(set_names(cmp_prefixes), function(prefix) {
  DE_stats_sub %>%
    rownames_to_column("Gene") %>%
    transmute(
      Gene,
      log2FC = .data[[paste0(prefix, ".Log2FoldChange")]],
      pvalue = .data[[paste0(prefix, ".RawPValue")]],
      padj = .data[[paste0(prefix, ".AdjustedPValue")]]
    )
})

fgsea_res_UNIFI <- imap(de_list, function(df, cmp_name) {
  ranks <- df %>%
    dplyr::filter(!is.na(padj), padj > 0) %>%
    mutate(rank = -log10(padj) * sign(log2FC)) %>%
    dplyr::select(Gene, rank) %>%
    deframe()

  fgsea(pathways = pathways, stats = ranks, nperm = 10000) %>%
    as_tibble() %>%
    mutate(comparison = cmp_name)
}) %>%
  bind_rows()

fgsea_res_UNIFI %>%
  dplyr::filter(comparison == cmp_prefixes[1]) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10)

fgsea_res_UNIFI$leadingEdge <- NULL
fgsea_res_UNIFI <- fgsea_res_UNIFI %>% dplyr::rename("dataset" = comparison)

# write.csv(fgsea_res_UNIFI, file.path(fe_results_dir, "RNR_targeted_signatures_gsea", "UNIFI_GSEA_scores.csv"))

# ------------------------------------------------------------------
# Include enrichment results from other studies including Varsity
# ------------------------------------------------------------------
fgsea_res_UNIFI <- read.table(unifi_gsea_scores_file, sep = ",", header = TRUE, row.names = 1)

fgsea_res_UNIFI <- fgsea_res_UNIFI %>%
  dplyr::filter(dataset %in% c("NoRemission_w8_NoMH_Ust130mpk_vs_Remission_MH_w8"))
fgsea_res_UNIFI$dataset <- gsub("NoRemission_w8_NoMH_Ust130mpk_vs_Remission_MH_w8", "GSE206285_pre_UST", fgsea_res_UNIFI$dataset)

fgsea_res <- read.table(hallmark_filtered_file, sep = ",", header = TRUE)
fgsea_res$dataset <- gsub("GSE73661_pre_VDZ", "GEMINI_pre_VDZ", fgsea_res$dataset)
fgsea_res$dataset <- gsub("Varsity_pre_VDZ", "VARSITY_pre_VDZ", fgsea_res$dataset)
fgsea_res$dataset <- gsub("Varsity_pre_ADA", "VARSITY_pre_ADA", fgsea_res$dataset)

de_names <- c("GSE206285_pre_UST", "GSE12251_pre_IFX", "GSE16879_pre_IFX", "GSE23597_pre_IFX", "GSE73661_pre_IFX", "VARSITY_pre_ADA", "VARSITY_pre_VDZ", "GEMINI_pre_VDZ")

fgsea_res <- rbind(fgsea_res, fgsea_res_UNIFI)
fgsea_res <- fgsea_res %>% dplyr::filter(padj <= 0.05, dataset %in% de_names)
fgsea_res$dataset <- factor(fgsea_res$dataset, levels = de_names)

pp <- ggplot(fgsea_res, aes(x = pathway, y = dataset)) +
  geom_point(aes(
    color = NES,
    size = -log10(padj)
  )) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-3, 3),
    oob = scales::squish,
    name = "NES"
  ) +
  scale_size_continuous(name = "-log10(padj)") +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  ggtitle("GSEA enrichment across IBD RNR studies - NR-up/Filtered")

# CairoPDF(file = file.path(project_dir, "Figs/Manuscript_EPS/Fig2_new/FIG2D_Angio_Validation_filtered_UST.pdf"), width = 9, height = 7)
# plot(pp)
# dev.off()

# ------------------------------------------------------------------
# Create a heatmap instead of the dotplot
# ------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(grid)

fgsea_res_UNIFI <- read.table(
  unifi_gsea_scores_file,
  sep = ",", header = TRUE, row.names = 1
)

fgsea_res_UNIFI <- fgsea_res_UNIFI %>% dplyr::filter(!pathway %in% c("Neutrophil_degranulation"))

fgsea_res_UNIFI <- fgsea_res_UNIFI %>%
  dplyr::filter(dataset %in% c("NoRemission_w8_NoMH_Ust130mpk_vs_Remission_MH_w8")) %>%
  mutate(dataset = gsub(
    "NoRemission_w8_NoMH_Ust130mpk_vs_Remission_MH_w8",
    "GSE206285_pre_UST", dataset
  ))

fgsea_res <- read.table(
  hallmark_filtered_file,
  sep = ",", header = TRUE
)
fgsea_res <- fgsea_res %>% dplyr::filter(!pathway %in% c("Neutrophil_degranulation"))

fgsea_res$dataset <- gsub("GSE73661_pre_VDZ", "GEMINI_pre_VDZ", fgsea_res$dataset)
fgsea_res$dataset <- gsub("Varsity_pre_VDZ", "VARSITY_pre_VDZ", fgsea_res$dataset)
fgsea_res$dataset <- gsub("Varsity_pre_ADA", "VARSITY_pre_ADA", fgsea_res$dataset)

de_names <- c(
  "GEMINI_pre_VDZ", "VARSITY_pre_VDZ", "VARSITY_pre_ADA",
  "GSE12251_pre_IFX", "GSE16879_pre_IFX", "GSE23597_pre_IFX",
  "GSE73661_pre_IFX", "GSE206285_pre_UST"
)

fg_all <- bind_rows(fgsea_res, fgsea_res_UNIFI) %>%
  dplyr::filter(dataset %in% de_names)

sig_paths <- fg_all %>%
  dplyr::filter(!is.na(padj) & padj < 0.05) %>%
  pull(pathway) %>%
  unique()

fg_use <- fg_all %>%
  dplyr::filter(pathway %in% sig_paths)

fg_use$dataset <- factor(fg_use$dataset, levels = de_names)

fg_wide <- fg_use %>%
  dplyr::select(dataset, pathway, NES, padj) %>%
  complete(dataset, pathway) %>%
  arrange(pathway, dataset)

nes_mat <- fg_wide %>%
  dplyr::select(dataset, pathway, NES) %>%
  pivot_wider(names_from = dataset, values_from = NES) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

padj_mat <- fg_wide %>%
  dplyr::select(dataset, pathway, padj) %>%
  pivot_wider(names_from = dataset, values_from = padj) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

nes_mat_capped <- nes_mat
nes_mat_capped[!is.na(nes_mat_capped) & nes_mat_capped > 3] <- 2
nes_mat_capped[!is.na(nes_mat_capped) & nes_mat_capped < -3] <- -2

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

ht <- Heatmap(
  nes_mat_capped,
  name = "NES",
  col = col_fun,
  na_col = "grey90",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 10),
  row_names_gp = gpar(fontsize = 9),
  heatmap_legend_param = list(title = "NES", at = c(-2, 0, 2)),
  use_raster = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    v_raw <- nes_mat[i, j]
    p <- padj_mat[i, j]
    v_cap <- nes_mat_capped[i, j]

    if (!is.na(v_cap) && (is.na(p) || p >= 0.05)) {
      grid.rect(
        x = x, y = y, width = width, height = height,
        gp = gpar(fill = "grey85", col = NA)
      )
    }

    if (!is.na(v_raw) && !is.na(p) && p < 0.05) {
      lab <- sprintf("%.2f*", v_raw)
      grid.text(lab, x = x, y = y, gp = gpar(fontsize = 9, col = "black"))
    }
  }
)

draw(ht)
write_xlsx(fg_all, path = fig2e_gsea_out)
