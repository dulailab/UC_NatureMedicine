## Setup ----------------------------------------
### Bioconductor and CRAN libraries used

library(tidyverse)
library(reshape2)
library(DESeq2)


## DGE analysis in R and NR, VDZ and IFX - baseline vs later time points ------------------------------------------------------------

comparisons <- c("4h", "24h", "72h", "2w", "6w", "14w")
groups <- c("R", "NR")
therapies <- c("Infliximab", "Vedolizumab")

for(therapy in therapies){
  for(group in groups){
    for(comparison in comparisons){
      
      #Initialize parameters
      condition1 <- '0h'
      condition2 <- comparison
      condition_column_name <- 'Time.Point'
      patient_column_name <- 'Patient.Nr.'
      file_column_name <- 'file_name'
      treatment <- c('Vedolizumab', 'Infliximab')
      filtering_options <- data.frame(column_name = c('Treatment', 'Remission.NR.C'), values = c(therapy, group))
      
      #load data
      count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
      col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
      
      analysis_folder <- file.path("results/tp_comparison", paste("0-", condition2, "_comparison", sep = ''))
      output_folder <- file.path(analysis_folder, ifelse(filtering_options$values[1] == '', '',
                                                        (ifelse(filtering_options$values[1] == 'Infliximab', 'IFX/', 'VDZ/'))))
      
      #Get patients for conditions
      pat_use1 <- subset(col_data, col_data[condition_column_name] == condition1)[patient_column_name]
      pat_use2 <- subset(col_data, col_data[condition_column_name] == condition2)[patient_column_name]
      pat_use <- intersect(pat_use1[,1], pat_use2[,1])
      
      #Create sample names
      for(i in 1:length(col_data[,condition_column_name])){
        trt <- ifelse(col_data[i, as.vector(filtering_options$column_name[1])]=='Infliximab', 'IFX', 
                      ifelse(col_data[i, as.vector(filtering_options$column_name[1])]=='Vedolizumab', 'VDZ', 'TZZ'))
        col_data$colnames[i] <- paste(col_data[i, condition_column_name], trt, 
                                      col_data[i, patient_column_name], sep = "_")
      }
      
      #Filter data for selected comparison
      col_data <- subset(col_data, col_data[,patient_column_name] %in% pat_use)
      col_data <- subset(col_data, col_data[condition_column_name] == condition1 | col_data[condition_column_name] == condition2)
      for(i in 1:length(filtering_options$column_name)){
        filter_column <- as.vector(filtering_options$column_name[i])
        value <- filtering_options$values[i]
        if(value == '' & filter_column == 'Treatment'){
          value <- treatment
        }
        col_data <- subset(col_data, col_data[,filter_column] %in% value)
      }
      
      col_data[,patient_column_name] <- factor(as.character(col_data[,patient_column_name]))
      col_data[,condition_column_name] <- factor(col_data[,condition_column_name])
      count_data <- count_data[, as.character(col_data[,file_column_name])]
      
      #run DGE analysis
      dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Patient.Nr. + Time.Point)
      dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
      dds_counts <- estimateSizeFactors(dds_counts)
      
      dds <- DESeq(dds_counts, betaPrior = FALSE)
      res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
      
      #Output results to file
      res_sorted <- res[order(res$padj), ]
      write.table(res_sorted, file = paste(output_folder, "/DESeq2result_", condition1, "vs", condition2, "_", filtering_options$values[2], ".txt", sep = ""), sep = "\t", quote = FALSE)
      
    }
  }
}




## DGE analysis in VDZ - R vs NR at all time points ------------------------------------------------------------

#load data
count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment == "Vedolizumab")
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point)
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

#run pairwise DGE analysis
timepoints <- c("0h", "4h", "24h", "72h", "2w", "6w", "14w")

for (tp in timepoints) {
  
  tp_col_data <- subset(col_data, col_data$Time.Point == tp)
  tp_count_data <- count_data[, as.character(tp_col_data$file_name)]
  
  dds_counts <- DESeqDataSetFromMatrix(countData = tp_count_data, colData = tp_col_data, design = ~ Remission.NR.C)
  dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
  dds_counts <- estimateSizeFactors(dds_counts)
  
  dds <- DESeq(dds_counts, betaPrior = FALSE)
  res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
  
  res_sorted <- res[order(res$padj), ]
  write.table(res_sorted, file = paste("results/R_NR_comparison/", tp, "/VDZ/DESeq2_result_NR_vs_R.txt", sep = ""), sep = "\t", quote = FALSE)
}




## DGE analysis in IFX - R vs NR at all time points ------------------------------------------------------------

#load data
count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
col_data <- read.table("metadata/Sequenced_sample_information.csv", header = TRUE, sep = ',')
col_data <- subset(col_data, col_data$Treatment == "Infliximab")
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point)
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$Active.Inactive <- factor(col_data$Active.Inactive)

#run pairwise DGE analysis
timepoints <- c("0h", "4h", "24h", "72h", "2w", "6w", "14w")

for (tp in timepoints) {
  
  tp_col_data <- subset(col_data, col_data$Time.Point == tp)
  tp_count_data <- count_data[, as.character(tp_col_data$file_name)]

  dds_counts <- DESeqDataSetFromMatrix(countData = tp_count_data, colData = tp_col_data, design = ~ Remission.NR.C)
  dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
  dds_counts <- estimateSizeFactors(dds_counts)
  
  dds <- DESeq(dds_counts, betaPrior = FALSE)
  res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
  
  res_sorted <- res[order(res$padj), ]
  write.table(res_sorted, file = paste("results/R_NR_comparison/", tp, "/IFX/DESeq2_result_NR_vs_R.txt", sep = ""), sep = "\t", quote = FALSE)
}

