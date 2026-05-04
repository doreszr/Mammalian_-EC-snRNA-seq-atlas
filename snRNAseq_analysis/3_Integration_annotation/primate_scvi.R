########################################################################################################################

### Author: Dorottya Ralbovszki
### Date: 08-24

### Integrating human and renamed baboon (renamed to human orthologues) samples with scVI

### saving Seurat objects for further analysis

########################################################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(reticulate)
library(SeuratWrappers)


#### Load data
human <- readRDS("Data_EC_atlas/Seurat_obj/human/merged.human_annot_diet.rds")

baboon <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamed.merged.baboon_annot_diet.rds")

##### Merging objects #####

# merge for integrative layers
merged.primate <- merge(x = human, y = baboon, project = "PRIMATE_EC")


##### Normalization and scaling #####

DefaultAssay(merged.primate.scvi) <- "RNA"

merged.primate.scvi <- NormalizeData(merged.primate.scvi)
merged.primate.scvi <- FindVariableFeatures(merged.primate.scvi)
merged.primate.scvi <- ScaleData(merged.primate.scvi)
merged.primate.scvi <- RunPCA(merged.primate.scvi)
merged.primate.scvi <- FindNeighbors(merged.primate.scvi, dims = 1:30, reduction = "pca")
merged.primate.scvi <- FindClusters(merged.primate.scvi, resolution = 0.6, cluster.name = "unintegrated_clusters")
merged.primate.scvi <- RunUMAP(merged.primate.scvi, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


##### Integration #####

merged.primate.scvi <- IntegrateLayers(
  object = merged.primate.scvi, method = scVIIntegration, orig.reduction = "pca",
  new.reduction = "integrated.scvi",
  conda_env = "PATH_TO_YOUR_CONDA_ENV", verbose = T)

##### Find neighbors and clusters

merged.primate.scvi <- FindNeighbors(merged.primate.scvi, reduction = "integrated.scvi", dims = 1:30)

# use the same resolution as for other methods for better comparison

merged.primate.scvi <- FindClusters(merged.primate.scvi, resolution = 0.6, cluster.name = "scvi_clusters")

## run UMAP on SCVI integration

merged.primate.scvi <- RunUMAP(merged.primate.scvi, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")

saveRDS(merged.primate.scvi, file = "Data_EC_atlas/Seurat_obj/primate/merged.primate.scvi.rds")
