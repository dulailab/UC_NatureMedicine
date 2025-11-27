## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(openxlsx)
library(rstatix)
library(DESeq2)


## Loading data ------------------------------------------------------------

tx2gene <- read.table("metadata/tx2gene.txt")
annotation <- tx2gene[, 2:3] %>%
  unique()

trajectory_genes <- read.xlsx("metadata/Varsity_trajectory_signatureV4.xlsx")

my_colors <- c("R" = "#4393C3", "NR" = "#D6604D")



## Load DGE analysis results ------------

pathways <- trajectory_genes$Name
results <- colnames(trajectory_genes)[2:6]
comparisons <- c("0h", "4h", "24h", "72h", "2w", "6w", "14w")

gene_table <- data.frame()

pathway <- pathways[1]
comparison <- comparisons[1]
result <- results[1]
for(comparison in comparisons){
  
  df_results <- read.table(paste0("results/R_NR_comparison/", comparison, "/VDZ/DESeq2_result_NR_vs_R.txt")) %>%
    rownames_to_column("Gene_ID") %>%
    merge(annotation, by.x = "Gene_ID", by.y = "gene_id")
  
  for(pathway in pathways){
    for(result in results){
      
      df_genes <- trajectory_genes[trajectory_genes$Name == pathway, result] %>%
        str_remove_all(" ") %>%
        str_split(",") %>% 
        unlist()
      
      df_results_temp <- df_results %>%
        subset(.$gene_name %in% df_genes) %>%
        add_column("pathway" = pathway,
                   "comparison" = comparison,
                   "timepoint" = result,
                   "n" = nrow(.))
      
      gene_table <- rbind(gene_table, df_results_temp)
    }
  }
}

gene_table[is.na(gene_table$padj), ]$padj <- 1




## Longitudinal expression of angiogenesis genes-----------------

### CD and UC, VDZ -----------------

#load data
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment == "Vedolizumab")
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point, levels = )
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
count_data <- count_data[, col_data$file_name]

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Active.Inactive)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

normalized_counts <- counts(dds_counts, normalized = TRUE)
vst_counts <- dds_counts %>%
  vst() %>%
  assay() %>%
  as.data.frame()

#select angiogenesis pathway
my_pathway <- "ANGIOGENESIS"

#prepare data frame with counts 
gene_table_temp <- gene_table %>%
  subset(.$timepoint == "NR_up_EndoW14NR") %>%
  .[order(.$padj, decreasing = FALSE), ] %>%
  subset(.$pathway == my_pathway)

my_gene_ids <- gene_table_temp %>%
  .$Gene_ID %>%
  unique()

df_plot <- vst_counts[my_gene_ids, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "file_name") %>%
  merge(col_data, by = "file_name") %>%
  merge(gene_table[, c("Gene_ID", "gene_name")], by = "Gene_ID")
df_plot$Time.Point <- factor(df_plot$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

#prepare data frame with median expression values
df_median <- df_plot %>%
  group_by(Remission.NR.C, gene_name, Time.Point) %>%
  dplyr::summarize("value" = median(value))
df_median$Time.Point <- factor(df_median$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

#prepare data frame with significance stars
df_max <- df_plot %>%
  group_by(gene_name) %>%
  dplyr::summarize("max" = max(value, na.rm = TRUE) + (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))*0.05)

significance <- gene_table_temp %>%
  subset(.$Gene_ID %in% my_gene_ids) %>%
  add_significance(p.col = "padj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  merge(df_max, by = "gene_name")
significance$comparison <- factor(significance$comparison, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

my_width <- ceiling(nrow(df_max)/3)

p <- ggplot(df_plot) +
  geom_line(aes(group = Patient.Nr., x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_point(data = df_median, size = 2, aes(x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_line(data = df_median, aes(group = Remission.NR.C,x = Time.Point, y = value, color = Remission.NR.C), lwd = 2, alpha = 0.5) +
  geom_text(data = significance, aes(x = comparison, label = p.adj_signif, y = max)) +
  scale_color_manual(values = my_colors) +
  facet_wrap(.~gene_name, nrow = 3, scales = "free_y") +
  labs(x = "Timepoint", y = "Variance stabilized normalized expression", title = my_pathway) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pdf(paste0("results/longitudinal_geneExpr_", my_pathway, "_CD_UC_Vedo.pdf"), width = my_width*3.5, height = 8)
plot(p)
dev.off()





### UC, VDZ -----------------

#load data
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment == "Vedolizumab")
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))
col_data <- subset(col_data, col_data$Diagnosis == "CU")

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point, levels = )
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
count_data <- count_data[, col_data$file_name]

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Active.Inactive)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

normalized_counts <- counts(dds_counts, normalized = TRUE)
vst_counts <- dds_counts %>%
  vst() %>%
  assay() %>%
  as.data.frame()

#select angiogenesis pathway
my_pathway <- "ANGIOGENESIS"

#prepare data frame with counts
gene_table_temp <- gene_table %>%
  subset(.$timepoint == "NR_up_EndoW14NR") %>%
  .[order(.$padj, decreasing = FALSE), ] %>%
  subset(.$pathway == my_pathway)

my_gene_ids <- gene_table_temp %>%
  .$Gene_ID %>%
  unique()

#prepare data frame with median expression values
df_plot <- vst_counts[my_gene_ids, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "file_name") %>%
  merge(col_data, by = "file_name") %>%
  merge(gene_table[, c("Gene_ID", "gene_name")], by = "Gene_ID")
df_plot$Time.Point <- factor(df_plot$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

df_median <- df_plot %>%
  group_by(Remission.NR.C, gene_name, Time.Point) %>%
  summarize("value" = median(value))
df_median$Time.Point <- factor(df_median$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

my_width <- ceiling(length(unique(df_plot$gene_name))/3)

p <- ggplot(df_plot) +
  geom_line(aes(group = Patient.Nr., x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_point(data = df_median, size = 2, aes(x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_line(data = df_median, aes(group = Remission.NR.C,x = Time.Point, y = value, color = Remission.NR.C), lwd = 2, alpha = 0.5) +
  scale_color_manual(values = my_colors) +
  facet_wrap(.~gene_name, nrow = 3, scales = "free_y") +
  labs(x = "Timepoint", y = "Variance stabilized normalized expression", title = my_pathway) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pdf(paste0("results/longitudinal_geneExpr_", my_pathway, "_UC_Vedo.pdf"), width = my_width*3.5, height = 8)
plot(p)
dev.off()





### UC, IFX -----------------

#load data
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment %in% c("Infliximab", ""))
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))
col_data <- subset(col_data, col_data$Diagnosis == "CU")

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point, levels = )
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
count_data <- count_data[, col_data$file_name]

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Active.Inactive)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

normalized_counts <- counts(dds_counts, normalized = TRUE)
vst_counts <- dds_counts %>%
  vst() %>%
  assay() %>%
  as.data.frame()

#select angiogenesis pathway
my_pathway <- "ANGIOGENESIS"

#prepare data frame with counts
gene_table_temp <- gene_table %>%
  subset(.$timepoint == "NR_up_EndoW14NR") %>%
  .[order(.$padj, decreasing = FALSE), ] %>%
  subset(.$pathway == my_pathway)

my_gene_ids <- gene_table_temp %>%
  .$Gene_ID %>%
  unique()

df_plot <- vst_counts[my_gene_ids, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "file_name") %>%
  merge(col_data, by = "file_name") %>%
  merge(gene_table[, c("Gene_ID", "gene_name")], by = "Gene_ID")
df_plot$Time.Point <- factor(df_plot$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

#prepare data frame with median expression values
df_median <- df_plot %>%
  group_by(Remission.NR.C, gene_name, Time.Point) %>%
  summarize("value" = median(value))
df_median$Time.Point <- factor(df_median$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

my_width <- ceiling(length(unique(df_plot$gene_name))/3)
p <- ggplot(df_plot) +
  geom_line(aes(group = Patient.Nr., x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_point(data = df_median, size = 2, aes(x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_line(data = df_median, aes(group = Remission.NR.C,x = Time.Point, y = value, color = Remission.NR.C), lwd = 2, alpha = 0.5) +
  scale_color_manual(values = my_colors) +
  facet_wrap(.~gene_name, nrow = 3, scales = "free_y") +
  labs(x = "Timepoint", y = "Variance stabilized normalized expression", title = my_pathway) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pdf(paste0("results/longitudinal_geneExpr_", my_pathway, "_UC_IFX.pdf"), width = my_width*3.5, height = 8)
plot(p)
dev.off()





## Longitudinal expression of IL6-JAK/STAT signaling genes-----------------

### CD and UC, VDZ -----------------

#load data
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment == "Vedolizumab")
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point, levels = )
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
count_data <- count_data[, col_data$file_name]

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Active.Inactive)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

normalized_counts <- counts(dds_counts, normalized = TRUE)
vst_counts <- dds_counts %>%
  vst() %>%
  assay() %>%
  as.data.frame()

#select IL6-JAK/STAT signaling pathway 
my_pathway <- "IL6_JAK_STAT3_SIGNALING"

#prepare data frame with counts 
gene_table_temp <- gene_table %>%
  subset(.$timepoint == "NR_up_EndoW14NR") %>%
  .[order(.$padj, decreasing = FALSE), ] %>%
  subset(.$pathway == my_pathway)

my_gene_ids <- gene_table_temp %>%
  subset(.$padj < 0.05) %>%
  .$Gene_ID %>%
  unique()

df_plot <- vst_counts[my_gene_ids, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "file_name") %>%
  merge(col_data, by = "file_name") %>%
  merge(gene_table[, c("Gene_ID", "gene_name")], by = "Gene_ID")
df_plot$Time.Point <- factor(df_plot$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

#prepare data frame with median expression values
df_median <- df_plot %>%
  group_by(Remission.NR.C, gene_name, Time.Point) %>%
  dplyr::summarize("value" = median(value))
df_median$Time.Point <- factor(df_median$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

#prepare data frame with significance stars
df_max <- df_plot %>%
  group_by(gene_name) %>%
  dplyr::summarize("max" = max(value, na.rm = TRUE) + (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))*0.05)

significance <- gene_table_temp %>%
  subset(.$Gene_ID %in% my_gene_ids) %>%
  add_significance(p.col = "padj", 
                   output.col = "p.adj_signif", 
                   cutpoints = c(0, 0.0001, 0.01, 0.05, 1),
                   symbols = c("***", "**", "*", "ns")) %>%
  merge(df_max, by = "gene_name")
significance$comparison <- factor(significance$comparison, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

my_width <- ceiling(nrow(df_max)/3)

p <- ggplot(df_plot) +
  geom_line(aes(group = Patient.Nr., x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_point(data = df_median, size = 2, aes(x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_line(data = df_median, aes(group = Remission.NR.C,x = Time.Point, y = value, color = Remission.NR.C), lwd = 2, alpha = 0.5) +
  geom_text(data = significance, aes(x = comparison, label = p.adj_signif, y = max)) +
  scale_color_manual(values = my_colors) +
  facet_wrap(.~gene_name, nrow = 3, scales = "free_y") +
  labs(x = "Timepoint", y = "Variance stabilized normalized expression", title = my_pathway) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())


pdf(paste0("results/longitudinal_geneExpr_", my_pathway, "_CD_UC_Vedo.pdf"), width = my_width*3.5, height = 8)
plot(p)
dev.off()




### UC, VDZ -----------------

#load data
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment == "Vedolizumab")
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))
col_data <- subset(col_data, col_data$Diagnosis == "CU")

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point, levels = )
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
count_data <- count_data[, col_data$file_name]

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Active.Inactive)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

normalized_counts <- counts(dds_counts, normalized = TRUE)
vst_counts <- dds_counts %>%
  vst() %>%
  assay() %>%
  as.data.frame()

#select IL6-JAK/STAT signaling pathway 
my_pathway <- "IL6_JAK_STAT3_SIGNALING"

#prepare data frame with counts 
gene_table_temp <- gene_table %>%
  subset(.$timepoint == "NR_up_EndoW14NR") %>%
  .[order(.$padj, decreasing = FALSE), ] %>%
  subset(.$pathway == my_pathway)

my_gene_ids <- gene_table_temp %>%
  subset(.$padj < 0.05) %>%
  .$Gene_ID %>%
  unique()

df_plot <- vst_counts[my_gene_ids, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "file_name") %>%
  merge(col_data, by = "file_name") %>%
  merge(gene_table[, c("Gene_ID", "gene_name")], by = "Gene_ID")
df_plot$Time.Point <- factor(df_plot$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

#prepare data frame with median expression values
df_median <- df_plot %>%
  group_by(Remission.NR.C, gene_name, Time.Point) %>%
  dplyr::summarize("value" = median(value))
df_median$Time.Point <- factor(df_median$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

my_width <- ceiling(length(unique(df_plot$gene_name))/3)

p <- ggplot(df_plot) +
  geom_line(aes(group = Patient.Nr., x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_point(data = df_median, size = 2, aes(x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_line(data = df_median, aes(group = Remission.NR.C,x = Time.Point, y = value, color = Remission.NR.C), lwd = 2, alpha = 0.5) +
  scale_color_manual(values = my_colors) +
  facet_wrap(.~gene_name, nrow = 3, scales = "free_y") +
  labs(x = "Timepoint", y = "Variance stabilized normalized expression", title = my_pathway) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pdf(paste0("results/longitudinal_geneExpr_", my_pathway, "_UC_Vedo.pdf"), width = my_width*3.5, height = 8)
plot(p)
dev.off()



### UC, IFX -----------------

#load data
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment %in% c("Infliximab", ""))
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))
col_data <- subset(col_data, col_data$Diagnosis == "CU")

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point, levels = )
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
count_data <- count_data[, col_data$file_name]

#normalize count data
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Active.Inactive)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

normalized_counts <- counts(dds_counts, normalized = TRUE)
vst_counts <- dds_counts %>%
  vst() %>%
  assay() %>%
  as.data.frame()

#select IL6-JAK/STAT signaling pathway 
my_pathway <- "IL6_JAK_STAT3_SIGNALING"

#prepare data frame with counts 
gene_table_temp <- gene_table %>%
  subset(.$timepoint == "NR_up_EndoW14NR") %>%
  .[order(.$padj, decreasing = FALSE), ] %>%
  subset(.$pathway == my_pathway)

my_gene_ids <- gene_table_temp %>%
  subset(.$padj < 0.05) %>%
  .$Gene_ID %>%
  unique()

df_plot <- vst_counts[my_gene_ids, ] %>%
  rownames_to_column("Gene_ID") %>%
  melt(variable.name = "file_name") %>%
  merge(col_data, by = "file_name") %>%
  merge(gene_table[, c("Gene_ID", "gene_name")], by = "Gene_ID")
df_plot$Time.Point <- factor(df_plot$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

#prepare data frame with median expression values
df_median <- df_plot %>%
  group_by(Remission.NR.C, gene_name, Time.Point) %>%
  dplyr::summarize("value" = median(value))
df_median$Time.Point <- factor(df_median$Time.Point, levels = c("0h", "4h", "24h", "72h", "2w", "6w", "14w"))

my_width <- ceiling(length(unique(df_plot$gene_name))/3)

p <- ggplot(df_plot) +
  geom_line(aes(group = Patient.Nr., x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_point(data = df_median, size = 2, aes(x = Time.Point, y = value, color = Remission.NR.C)) +
  geom_line(data = df_median, aes(group = Remission.NR.C,x = Time.Point, y = value, color = Remission.NR.C), lwd = 2, alpha = 0.5) +
  scale_color_manual(values = my_colors) +
  facet_wrap(.~gene_name, nrow = 3, scales = "free_y") +
  labs(x = "Timepoint", y = "Variance stabilized normalized expression", title = my_pathway) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pdf(paste0("results/longitudinal_geneExpr_", my_pathway, "_UC_IFX.pdf"), width = my_width*3.5, height = 8)
plot(p)
dev.off()









