library(dplyr)
library(Seurat)
library(ggplot2)
library(CellChat)
library(ComplexHeatmap)
library(SpatialExperiment)
library(nnSVG)

#PBMC datasets were download from the GEO accession at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261334
#the original publication is https://pubmed.ncbi.nlm.nih.gov/39343250/
#set data path
folder_path <- "~/GSE261334/"
# Get all .h5 file names
rds_files <- list.files(path = folder_path, pattern = "\\.h5$", full.names = TRUE)

# Create a vector of object names without file extension
object_names <- gsub("GSM[0-9]+_(.*?)_filtered_counts.h5", "\\1", basename(rds_files))

# read in all the list, the list will be used to generate seurat objects
seurat_list <- lapply(rds_files, Read10X_h5)

# Assign cleaned names to the list
names(seurat_list) <- object_names

# Check
print(object_names)
print(names(seurat_list))

#create seurat objects
# Convert the Read10X_h5 output list into Seurat objects using only the RNA (Gene Expression) assay
seurat_list_obj <- lapply(seurat_list, function(data) {
  expr_data <- data[["Gene Expression"]]
  ADT_data <- data[["Antibody Capture"]]
  
  # Create Seurat object with RNA
  seurat_obj <- CreateSeuratObject(counts = expr_data)
  
  # Create ADT assay and add to the object
  adt_assay <- CreateAssay5Object(counts = ADT_data)
  seurat_obj[["ADT"]] <- adt_assay
  
  return(seurat_obj)
})

#add sample names
for (name in names(seurat_list_obj)) {
  seurat_list_obj[[name]]$sample <- name
}

#remove the CD4 samples, we only need the PBMC datasets
seurat_list_obj <- seurat_list_obj[1:25]


#merge samples 
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_list_obj)

#do a simple QC by filtering the samples with nfeature >= 400
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
merged_seurat <- subset(merged_seurat, subset = percent.mt <= 20)
merged_seurat <- subset(merged_seurat, subset = nFeature_RNA >= 400)



#process the data
merged_seurat[["RNA"]] <- split(merged_seurat[["RNA"]], f = merged_seurat$sample)
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- ScaleData(merged_seurat)
merged_seurat <- RunPCA(merged_seurat)
merged_seurat <- IntegrateLayers(object = merged_seurat, orig.reduction = "pca",method = HarmonyIntegration,
  new.reduction = "harmony", verbose = TRUE)
merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
#need to re-do the clustering and cell typing
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30, reduction = "harmony")
merged_seurat <- FindClusters(merged_seurat, resolution = c(0.1,0.2,0.3,0.4))

#add condition information

#the repsonder non-responder information should be found here: https://github.com/VeroHo/vedo_paper/blob/main/scRNAseq/sample_sheet.txt
merged_seurat$condition <- case_when(merged_seurat$sample %in% unique(merged_seurat$sample)[1:5] ~ "HC",
                                     merged_seurat$sample %in% c("IBD-1_1_PBMC", "IBD-2_1_PBMC", "IBD-3_1_PBMC", "IBD-7_1_PBMC", "IBD-10_1_PBMC") ~ "Post_0week_Res",
                                     merged_seurat$sample %in% c("IBD-1_2_PBMC", "IBD-2_2_PBMC", "IBD-3_2_PBMC", "IBD-7_2_PBMC", "IBD-10_2_PBMC") ~ "Post_6week_Res",
                                     merged_seurat$sample %in% c("IBD-4_1_PBMC", "IBD-5_1_PBMC", "IBD-6_1_PBMC", "IBD-8_1_PBMC", "IBD-9_1_PBMC") ~ "Post_0week_nonRes",
                                     merged_seurat$sample %in% c("IBD-4_2_PBMC", "IBD-5_2_PBMC", "IBD-6_2_PBMC", "IBD-8_2_PBMC", "IBD-9_2_PBMC") ~ "Post_6week_nonRes")


#add the module score of interest


ANGIOGENESIS_NR_up_EndoW14NR <- "STC1,SPP1,FGFR1,TIMP1,CXCL6,THBD,SERPINA5,PGLYRP1,PRG2,OLR1,MSX1"
ANGIOGENESIS_NR_up_all <- "SPP1,STC1,TIMP1,CXCL6,OLR1,THBD,PGLYRP1,PRG2,FSTL1,KCNJ8,FGFR1,PF4,CCND2,MSX1"

ANGIOGENESIS_NR_up_EndoW14NR <- unlist(strsplit(ANGIOGENESIS_NR_up_EndoW14NR, ","))
ANGIOGENESIS_NR_up_all <- unlist(strsplit(ANGIOGENESIS_NR_up_all, ","))

merged_seurat <- AddModuleScore(merged_seurat, features = list(ANGIOGENESIS_NR_up_EndoW14NR,ANGIOGENESIS_NR_up_all),
                                      name = c("ANGIOGENESIS_NR_up_EndoW14NR","ANGIOGENESIS_NR_up_all"))

#annotat the cell types by using resolution = 0.2
Idents(merged_seurat) <- "RNA_snn_res.0.2"
merged_seurat$cell_type_res02 <- case_when(merged_seurat$RNA_snn_res.0.2 == 0 ~ "CD14 Mono",
                                           merged_seurat$RNA_snn_res.0.2 == 1 ~ "CD8 T",
                                           merged_seurat$RNA_snn_res.0.2 == 2 ~ "CD4 T",
                                           merged_seurat$RNA_snn_res.0.2 == 3 ~ "CD4 Memory T",
                                           merged_seurat$RNA_snn_res.0.2 == 4 ~ "B",
                                           merged_seurat$RNA_snn_res.0.2 == 5 ~ "NK",
                                           merged_seurat$RNA_snn_res.0.2 == 6 ~ "naive CD8 T",
                                           merged_seurat$RNA_snn_res.0.2 == 7 ~ "CD16 mono",
                                           merged_seurat$RNA_snn_res.0.2 == 8 ~ "CD14 mono",
                                           merged_seurat$RNA_snn_res.0.2 == 9 ~ "cDCs",
                                           merged_seurat$RNA_snn_res.0.2 == 10 ~ "MAIT",
                                           merged_seurat$RNA_snn_res.0.2 == 11 ~ "B",
                                           merged_seurat$RNA_snn_res.0.2 == 12 ~ "Proliferating T",
                                           merged_seurat$RNA_snn_res.0.2 == 13 ~ "non-classical monocytes",
                                           merged_seurat$RNA_snn_res.0.2 == 14 ~ "pDCs",
                                           merged_seurat$RNA_snn_res.0.2 == 15 ~ "plasma cell precursor",
                                           merged_seurat$RNA_snn_res.0.2 == 16 ~ "megakaryocytes",
                                           merged_seurat$RNA_snn_res.0.2 == 17 ~ "SPP1 MoMac"
                                           )

#subset to monocytes and re-process and re-cluster
Idents(merged_seurat) <- "RNA_snn_res.0.2"
merged_seurat_monocytes <- subset(merged_seurat, idents = c(0,7,8,9,13,17))

merged_seurat_monocytes[["RNA"]] <- split(merged_seurat_monocytes[["RNA"]], f = merged_seurat_monocytes$sample)
merged_seurat_monocytes <- NormalizeData(merged_seurat_monocytes)
merged_seurat_monocytes <- FindVariableFeatures(merged_seurat_monocytes)
merged_seurat_monocytes <- ScaleData(merged_seurat_monocytes)
merged_seurat_monocytes <- RunPCA(merged_seurat_monocytes)
merged_seurat_monocytes <- IntegrateLayers(object = merged_seurat_monocytes, orig.reduction = "pca",method = HarmonyIntegration,
  new.reduction = "harmony", verbose = TRUE)
merged_seurat_monocytes <- RunUMAP(merged_seurat_monocytes, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
#need to re-do the clustering and cell typing
merged_seurat_monocytes <- FindNeighbors(merged_seurat_monocytes, dims = 1:30, reduction = "harmony")
merged_seurat_monocytes <- FindClusters(merged_seurat_monocytes, resolution = c(0.1,0.2,0.3,0.4))


#label the label the monocyte subset by SPP1 expression
merged_seurat_monocytes <- JoinLayers(merged_seurat_monocytes)
#Fetch SPP1 expression as a data frame
df <- FetchData(merged_seurat_monocytes, vars = "SPP1")
df$SPP1_status <- ifelse(df$SPP1 > 0, "SPP1+", "SPP1-")
merged_seurat_monocytes$SPP1_status <- df$SPP1_status

#for figure 2F
#subset to only week 0 (baseline)
Idents(merged_seurat_monocytes) <- "condition"
merged_seurat_monocytes_week0 <- subset(merged_seurat_monocytes, idents = c("Post_0week_Res","Post_0week_nonRes"))

DotPlot(merged_seurat_monocytes_week0, group.by = "SPP1_monocytes", features = c("ANGIOGENESIS_NR_up_all2",ANGIOGENESIS_NR_up_all,"ITGA4","ITGB7", "MMP19", "CCL2"),cols = c("lightgrey","red"))
                        +ggtitle("Monocytes Baseline (week 0)")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
