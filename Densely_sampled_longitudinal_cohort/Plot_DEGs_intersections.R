## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(openxlsx)
library(ComplexHeatmap)
library(UpSetR)


## Loading data ------------------------------------------------------------

tx2gene <- read.table("metadata/tx2gene.txt")
annotation <- tx2gene[, 2:3] %>%
  unique()




## Load DGE analysis results ------------------------

#get files
comparisons <- c("4h", "24h", "72h", "2w", "6w", "14w")
groups <- c("R", "NR")
therapies <- c("IFX", "VDZ")

files <- c()
for(therapy in therapies){
  for(group in groups){
    for(comparison in comparisons){
      
      files[i] <- paste0("results/tp_comparison/0-", comparison, "_comparison/", therapy, "/DESeq2result_0hvs", comparison, "_", group, ".txt")
      names(files)[i] <- paste0(therapy, "_", comparison, "_0h_UC_CD_", group)
      i <- i+1
    }
  }
}

#load result tables
gene_number_data <- data.frame()
listInput <- list()
up_list <- list()
down_list <- list()
for(i in 1:length(files)){
  
  df_results <- read.table(files[i], header = TRUE) %>%
    rownames_to_column("gene_id") %>%
    merge(annotation, by = "gene_id")
  
  df_results <- df_results %>%
    subset(.$padj < 0.05)
  
  up_len <- df_results %>%
    subset(.$log2FoldChange > 0)
  down_len <- df_results %>%
    subset(.$log2FoldChange < 0)
  
  df_results <- df_results[order(df_results$pvalue, decreasing = FALSE), ]
  
  listInput[[i]] <- as.character(df_results$gene_id)
  up_list[[i]] <- as.character(up_len$gene_id)
  down_list[[i]] <- as.character(down_len$gene_id)
  
  gene_number_data <- rbind(gene_number_data, data.frame("Name" = names(files)[i], 
                                                         "Name" = nrow(up_len), 
                                                         "Name" = nrow(down_len)))
}

names(listInput) <- names(files)
names(up_list) <- names(files)
names(down_list) <- names(files)




## Upset plots for up- adn downregulated DEGs at baseline vs later time points for R and NR (IFX and VDZ) ----------------------------------

#upregulated degs (VDZ)
comb_matrix <- make_comb_mat(up_list[13:24], mode = "distinct")
df_upset_hyper <- UpSet(comb_matrix,
                        top_annotation = HeatmapAnnotation(
                          "DEG intersections" = 
                            anno_barplot(comb_size(comb_matrix), 
                                         border = FALSE, 
                                         gp = gpar(col = "#990000", fill = "#990000"), 
                                         height = unit(3, "cm"),
                                         numbers_rot = 90,
                                         add_numbers = TRUE), 
                          annotation_name_side = "left", 
                          gp = gpar(col = "black"),
                          annotation_name_rot = 0),
                        right_annotation = 
                          upset_right_annotation(comb_matrix, 
                                                 gp = gpar(col = "#990000", fill = "#990000"),
                                                 numbers_rot = 0,
                                                 add_numbers = TRUE),
                        comb_col = c("#990000"),
                        comb_order = order(comb_degree(comb_matrix), -comb_size(comb_matrix)))

pdf("results/VDZ_0h_vs_timepoint_upDEG_intersection.pdf", width = 9, height = 5)
df_upset_hyper
dev.off()


#downregulated degs (VDZ)
comb_matrix <- make_comb_mat(down_list[13:24], mode = "distinct")
df_upset_hypo <- UpSet(comb_matrix,
                       top_annotation = HeatmapAnnotation(
                         "DEG intersections" = 
                           anno_barplot(comb_size(comb_matrix), 
                                        border = FALSE, 
                                        gp = gpar(col = "#004C99", fill = "#004C99"), 
                                        height = unit(3, "cm"),
                                        numbers_rot = 90,
                                        add_numbers = TRUE), 
                         annotation_name_side = "left", 
                         gp = gpar(col = "black"),
                         annotation_name_rot = 0),
                       right_annotation = 
                         upset_right_annotation(comb_matrix, 
                                                gp = gpar(col = "#004C99", fill = "#004C99"),
                                                numbers_rot = 0,
                                                add_numbers = TRUE),
                       comb_col = c("#004C99"),
                       comb_order = order(comb_degree(comb_matrix), -comb_size(comb_matrix)))

pdf("results/VDZ_0h_vs_timepoint_downDEG_intersection.pdf", width = 9, height = 5)
df_upset_hypo
dev.off()


#upregulated degs (IFX)
comb_matrix <- make_comb_mat(up_list[1:12], mode = "distinct")
df_upset_hyper <- UpSet(comb_matrix,
                        top_annotation = HeatmapAnnotation(
                          "DEG intersections" = 
                            anno_barplot(comb_size(comb_matrix), 
                                         border = FALSE, 
                                         gp = gpar(col = "#990000", fill = "#990000"), 
                                         height = unit(3, "cm"),
                                         numbers_rot = 90,
                                         add_numbers = TRUE), 
                          annotation_name_side = "left", 
                          gp = gpar(col = "black"),
                          annotation_name_rot = 0),
                        right_annotation = 
                          upset_right_annotation(comb_matrix, 
                                                 gp = gpar(col = "#990000", fill = "#990000"),
                                                 numbers_rot = 0,
                                                 add_numbers = TRUE),
                        comb_col = c("#990000"),
                        comb_order = order(comb_degree(comb_matrix), -comb_size(comb_matrix)))

pdf("results/IFX_0h_vs_timepoint_upDEG_intersection.pdf", width = 9, height = 5)
df_upset_hyper
dev.off()


#downregulated degs (IFX)
comb_matrix <- make_comb_mat(down_list[1:12], mode = "distinct")
df_upset_hypo <- UpSet(comb_matrix,
                       top_annotation = HeatmapAnnotation(
                         "DEG intersections" = 
                           anno_barplot(comb_size(comb_matrix), 
                                        border = FALSE, 
                                        gp = gpar(col = "#004C99", fill = "#004C99"), 
                                        height = unit(3, "cm"),
                                        numbers_rot = 90,
                                        add_numbers = TRUE), 
                         annotation_name_side = "left", 
                         gp = gpar(col = "black"),
                         annotation_name_rot = 0),
                       right_annotation = 
                         upset_right_annotation(comb_matrix, 
                                                gp = gpar(col = "#004C99", fill = "#004C99"),
                                                numbers_rot = 0,
                                                add_numbers = TRUE),
                       comb_col = c("#004C99"),
                       comb_order = order(comb_degree(comb_matrix), -comb_size(comb_matrix)))

pdf("results/IFX_0h_vs_timepoint_downDEG_intersection.pdf", width = 9, height = 5)
df_upset_hypo
dev.off()
