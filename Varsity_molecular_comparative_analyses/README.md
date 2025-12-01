# Varsity Molecular Comparative Analyses

R scripts that generate volcano plots, fgsea results, heatmaps, and Olink expression trajectories for the VARSITY/GEMINI studies.

## Files
- `Endotypes_volcanoplot_GSEA_FigSup2.R`: Loads trajectory DEG tables, builds volcano plots across treatment/timepoint/control grids, and runs Hallmark fgsea/dotplots for Extended Fig 2. Reads `DGE_results/Trajectory_DEGs` and `files/hallmark_pathways_genes_all.csv`; exports to `Figs/Manuscript_EPS` when uncommented.
- `VARSITY_GEMINI_gsea_heatmap.R`: Assembles Varsity and GEMINI Hallmark fgsea results into NES/padj matrices and renders significance-masked heatmaps; outputs under `Figs/Manuscript_EPS` if export lines are enabled.
- `VDZ_ADA_Vascular_SPP1_lineplots.R`: Processes Olink NPX data and clinical metadata to plot responder vs non-responder trajectories (raw/change/% change) for SPP1 and a vascular inflammation signature. Uses `Fig2/Olink_data` and writes to `Fig2/results/figures`.
- `trajectory_DE.R`: Loads counts/FPKM/metadata from `files/`, runs PCA, DESeq2 contrasts, DEG summaries, GSEA, and figure exports under `Figs/Manuscript_EPS/Fig2_new`.
- `Fig2e_analysis.R`: Maps UNIFI contrasts to Hallmark/trajectory signatures, runs fgsea, and merges Varsity/GEMINI enrichment. Reads `files/UNIFI_data` and `files/RNR_targeted_signatures_gsea`; writes `GSEA_results/Fig2e_gsea.xlsx`.

## Inputs and directories
- DEG tables: `DGE_results/Trajectory_DEGs` and `DGE_results/Trajectory_DEGs_W14EndoNR`.
- fgsea outputs/targets: `GSEA_results/` and `FE_results/` (created on run), `files/RNR_targeted_signatures_gsea/`.
- UNIFI data: `files/UNIFI_data/` (DE stats, sample metadata, gene annotations, selected comparisons).
- Olink data: `Fig2/Olink_data/` (NPX matrix, sample manifest, clinical metadata).
- Varsity expression/metadata: `files/sample_metadata_complete.txt`, `files/sample_metadata_filtered.txt`, `files/clinical_data.RDS`, `files/Varsity_data.RDS`, plus `files/hallmark_pathways_genes_all.csv`.

## Dependencies (common)
`dplyr`, `tidyr`, `readr`, `readxl`, `writexl`, `ggplot2`, `ggpubr`, `ggfortify`, `ggvenn`, `ComplexHeatmap`, `circlize`, `EnhancedVolcano`, `fgsea`, `msigdbr`, `DESeq2`, `edgeR`, `GSVA`, `Cairo`, `glue`, `purrr`, `grid`, `reshape2/reshape`, `NMF`, `gprofiler2`, `DT`, `openxlsx`, `survival/ggsurvfit` (where used). Install via `install.packages(...)` as needed.

