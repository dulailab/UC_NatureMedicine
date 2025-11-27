## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(tximport)
library(DESeq2)
library(genefilter)
library(geneplotter)
library(topGO)
library(DBI)
library(plyr)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(writexl)

source("Functional_enrichment_analysis_helper_functions.R")


## Loading data ------------------------------------------------------------

tx2gene <- read.table("metadata/tx2gene.txt")
annotation <- tx2gene[, 2:3] %>%
  unique()




## ORA for DEGs in R vs NR at all time points (for VDZ and IFX) ----------------------

my_group <- c("VDZ", "IFX")
my_tp <- c("0h", "4h", "24h", "72h", "2w", "6w", "14w")
my_direction <- c("up", "down")

for(d in 1:length(my_direction)){
  for(g in 1:length(my_group)){
    for(t in 1:length(my_tp)){
      
      df_results <- read.table(c(paste0("results/R_NR_comparison/", my_tp[t], "/", my_group[g], "/DESeq2_result_NR_vs_R.txt")), header = TRUE)
      
      df_degs <- df_results %>%
        rownames_to_column("Gene_ID") %>%
        merge(annotation, by.x = "Gene_ID", by.y = "gene_id") %>%
        subset(.$padj < 0.05)
      
      if(my_direction[d] == "up"){
        degs <- df_degs %>%
          subset(.$log2FoldChange > 0) %>%
          .$Gene_ID %>%
          unique()
      }else if(my_direction[d] == "down"){
        degs <- df_degs %>%
          subset(.$log2FoldChange < 0) %>%
          .$Gene_ID %>%
          unique()
      }
      
      if(length(degs) > 2){
        
        overallBaseMean <- as.matrix(df_results[, "baseMean", drop = F]) 
        colnames(overallBaseMean) <- "mean_expression"
        
        topGOResults <- run_GO_ORA(degs, 
                                   overallBaseMean, 
                                   selected_database = "org.Hs.eg.db", 
                                   back_num = 20,
                                   annotation = annotation)
        
        go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
        
        file_filtered <- paste0("results/ORA_GO_R_vs_NR/GOterms_DEGs_RvsNR_", my_group[g], "_", my_tp[t], "_", my_direction[d], ".csv")
        write.csv(go_results_filtered, file_filtered, row.names = FALSE)
      }
    }
  }
}

  



### Dotplot for IFX and VDZ GO terms at baseline ---------------------

my_group <- c("VDZ", "IFX")
my_direction <- c("up", "down")
my_tp <- c("0h")

GO_tables <- list()

for(g in 1:length(my_group)){
  for(d in 1:length(my_direction)){
    
    GO_tables <- append(GO_tables,
                        list(read.csv(paste0("results/ORA_GO_R_vs_NR/GOterms_DEGs_RvsNR_", my_group[g], "_", my_tp, "_", my_direction[d], ".csv"))))
  }
}

names(GO_tables) <- paste0(rep(my_group, each = 2), "_", my_tp, "_", my_direction, "_in_R")

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = names(GO_tables),
                     score_name = score_name,
                     reverse_bool = FALSE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 58)

p <- ggplot(df_plot, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  ylab("GO Term") +
  scale_colour_gradient(high = "#990000", low = "#FF9999")  + 
  theme_bw() + 
  theme(axis.text.y = element_text(hjust = 1, size=13, color = "black"), 
        axis.text.x = element_text(size=13, color = "black", angle = 45, hjust = 1),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 13),
        axis.title = element_blank()) +
  labs(color = expression(paste("-log"[10], "p")), 
       size = paste0("Gene ratio"))

pdf(paste0("results/ORA_GO_R_vs_NR/dotplot_GOterms_DEGs_VDZ_IFX_biopsies_baseline.pdf"), height = 7, width = 7)
print(p)
dev.off()


#get DGE results for genes behind angiogenesis
df_go <- read.csv(paste0("results/ORA_GO_R_vs_NR/GOterms_DEGs_RvsNR_VDZ_0h_down.csv"))
degs <- read.table(c(paste0("results/R_NR_comparison/0h/VDZ/DESeq2_result_NR_vs_R.txt")), header = TRUE) %>%
  rownames_to_column("Gene_ID") %>%
  merge(annotation, by.x = "Gene_ID", by.y = "gene_id") 
  
angio_genes <- df_go %>%
  subset(.$Term == "angiogenesis") %>%
  .$Genes %>%
  str_split(", ") %>%
  unlist()

degs_ango <- degs %>%
  subset(.$Gene_ID %in% angio_genes) %>%
  .[order(.$pvalue, decreasing = FALSE), ]

write.csv(degs_ango, "results/ORA_GO_R_vs_NR/Angiogenesis_DEGs_RvsNR_VDZ_0h_down.csv", row.names = FALSE)





## ORA for DEGs in R and NR at baseline vs later time points (for VDZ and IFX) ----------------------

comparisons <- c("4h", "24h", "72h", "2w", "6w", "14w")
groups <- c("R", "NR")
therapies <- c("IFX", "VDZ")
my_direction <- c("up", "down")

for(therapy in therapies){
  for(group in groups){
    for(comparison in comparisons){
      for(d in 1:length(my_direction)){
        
        df_results <- read.table(paste("results/tp_comparison/", "0-", comparison, "_comparison/", therapy, "/DESeq2result_0hvs", comparison, "_", group, ".txt", sep = ""), header = TRUE)
        
        df_degs <- df_results %>%
          rownames_to_column("Gene_ID") %>%
          merge(annotation, by.x = "Gene_ID", by.y = "gene_id") %>%
          subset(.$padj < 0.05)
        
        if(my_direction[d] == "up"){
          degs <- df_degs %>%
            subset(.$log2FoldChange > 0) %>%
            .$Gene_ID %>%
            unique()
        }else if(my_direction[d] == "down"){
          degs <- df_degs %>%
            subset(.$log2FoldChange < 0) %>%
            .$Gene_ID %>%
            unique()
        }
        
        if(length(degs) > 2){
          
          overallBaseMean <- as.matrix(df_results[, "baseMean", drop = F]) 
          colnames(overallBaseMean) <- "mean_expression"
          
          topGOResults <- run_GO_ORA(degs, 
                                     overallBaseMean, 
                                     selected_database = "org.Hs.eg.db", 
                                     back_num = 20,
                                     annotation = annotation)
          
          go_results_filtered <- filter_go_results(topGOResults, top_num = 600, onto = "BP")
          
          file_filtered <- paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_", group, "_0hvs", comparison, "_", my_direction[d], ".csv")
          write.csv(go_results_filtered, file_filtered, row.names = FALSE)
        }
      }
    }
  }
}




### Dotplot for IFX and VDZ GO terms in baseline vs later time points for R ---------------------

#VDZ - up - R
therapy <- "VDZ"
my_direction <- "up"
group <- "R"
comparisons <- c("4h", "24h", "2w", "6w", "14w")

GO_tables <- list()
for(comparison in comparisons){
  
  GO_tables <- append(GO_tables,
                      list(read.csv(paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_", group, "_0hvs", comparison, "_", my_direction, ".csv"))))
}

names(GO_tables) <- paste0(therapy, "_", group, "_", comparisons, "vs0h_", my_direction)

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = names(GO_tables),
                     score_name = score_name,
                     reverse_bool = FALSE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 58)

p <- ggplot(df_plot, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  ylab("GO Term") +
  scale_colour_gradient(high = "#990000", low = "#FF9999")  + 
  theme_bw() + 
  theme(axis.text.y = element_text(hjust = 1, size=13, color = "black"), 
        axis.text.x = element_text(size=13, color = "black", angle = 45, hjust = 1),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 13),
        axis.title = element_blank()) +
  labs(color = expression(paste("-log"[10], "p")), 
       size = paste0("Gene ratio"))

pdf(paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_biopsies_", group, "_tpvs0h_", my_direction, ".pdf"), height = 8, width = 8)
print(p)
dev.off()




#VDZ - down - R 
therapy <- "VDZ"
my_direction <- "down"
group <- "R"
comparisons <- c("4h", "2w", "6w", "14w")

GO_tables <- list()
for(comparison in comparisons){
  
  GO_tables <- append(GO_tables,
                      list(read.csv(paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_", group, "_0hvs", comparison, "_", my_direction, ".csv"))))
}

names(GO_tables) <- paste0(therapy, "_", group, "_", comparisons, "vs0h_", my_direction)

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = names(GO_tables),
                     score_name = score_name,
                     reverse_bool = FALSE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 58)

p <- ggplot(df_plot, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  ylab("GO Term") +
  scale_colour_gradient(high = "#990000", low = "#FF9999")  + 
  theme_bw() + 
  theme(axis.text.y = element_text(hjust = 1, size=13, color = "black"), 
        axis.text.x = element_text(size=13, color = "black", angle = 45, hjust = 1),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 13),
        axis.title = element_blank()) +
  labs(color = expression(paste("-log"[10], "p")), 
       size = paste0("Gene ratio"))

pdf(paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_biopsies_", group, "_tpvs0h_", my_direction, ".pdf"), height = 8, width = 8)
print(p)
dev.off()




#IFX - up - R 
therapy <- "IFX"
my_direction <- "up"
group <- "R"
comparisons <- c("24h", "72h", "2w", "6w", "14w")

GO_tables <- list()
for(comparison in comparisons){
  
  GO_tables <- append(GO_tables,
                      list(read.csv(paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_", group, "_0hvs", comparison, "_", my_direction, ".csv"))))
}

names(GO_tables) <- paste0(therapy, "_", group, "_", comparisons, "vs0h_", my_direction)

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = names(GO_tables),
                     score_name = score_name,
                     reverse_bool = FALSE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 58)

p <- ggplot(df_plot, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  ylab("GO Term") +
  scale_colour_gradient(high = "#990000", low = "#FF9999")  + 
  theme_bw() + 
  theme(axis.text.y = element_text(hjust = 1, size=13, color = "black"), 
        axis.text.x = element_text(size=13, color = "black", angle = 45, hjust = 1),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 13),
        axis.title = element_blank()) +
  labs(color = expression(paste("-log"[10], "p")), 
       size = paste0("Gene ratio"))

pdf(paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_biopsies_", group, "_tpvs0h_", my_direction, ".pdf"), height = 8, width = 7)
print(p)
dev.off()



#IFX - down - R
therapy <- "IFX"
my_direction <- "down"
group <- "R"
comparisons <- c("24h", "72h", "2w", "6w", "14w")

GO_tables <- list()
for(comparison in comparisons){
  
  GO_tables <- append(GO_tables,
                      list(read.csv(paste0("ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_", group, "_0hvs", comparison, "_", my_direction, ".csv"))))
}

names(GO_tables) <- paste0(therapy, "_", group, "_", comparisons, "vs0h_", my_direction)

score_name <- "Fisher.elim"
df_plot <- get_GO_df(GO_tables = GO_tables,
                     direction = names(GO_tables),
                     score_name = score_name,
                     reverse_bool = FALSE,
                     showTerms = 5,
                     multLines = TRUE,
                     numChar = 58)

p <- ggplot(df_plot, mapping = aes(x = Direction, y = Term, size = Ratio, color = changed_scores)) +
  geom_point() +
  ylab("GO Term") +
  scale_colour_gradient(high = "#990000", low = "#FF9999")  + 
  theme_bw() + 
  theme(axis.text.y = element_text(hjust = 1, size=13, color = "black"), 
        axis.text.x = element_text(size=13, color = "black", angle = 45, hjust = 1),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.placement = "outside", 
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size = 13),
        axis.title = element_blank()) +
  labs(color = expression(paste("-log"[10], "p")), 
       size = paste0("Gene ratio"))

pdf(paste0("results/ORA_GO_baseline_vs_timepoints/GOterms_DEGs_", therapy, "_biopsies_", group, "_tpvs0h_", my_direction, ".pdf"), height = 9, width = 8)
print(p)
dev.off()

