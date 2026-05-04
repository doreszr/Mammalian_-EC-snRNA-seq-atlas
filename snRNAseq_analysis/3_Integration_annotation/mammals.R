########################################################################################################################

### Author: Dorottya Ralbovszki

### Date: 05-25

### Integration of all mammalian neuronal data sets

### Neuronal cell type annotation

### saving Seurat objects for further analysis

########################################################################################################################


library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(clustree)
library(cowplot)
library(Matrix)
library(sctransform)
library(glmGamPoi)
library(ROCR)
library(parallel)
library(fields)
library(presto)


##### Loading data #####
baboon <- readRDS(file = "Data_EC_atlas/Seurat_obj/primate/rerenamed.baboon.neu_annot_diet.rds")
human <- readRDS(file = "/Data_EC_atlas/Seurat_obj/primate/renamed.human.neu_annot_diet.rds")
mouse <- readRDS(file = "Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot_diet.rds")
bat <- readRDS(file = "Data_EC_atlas/Seurat_obj/fruitbaat/renamed.merged.bat_neu_annot_diet.rds")


##### Merging data sets #####
merged.mammal <- merge(human, y = c(mouse, bat, baboon), project = "Human_Baboon_Mouse_Bat_EC")

##### Rerunning normalization with SCTransform #####
merged.mammal <- SCTransform(merged.mammal, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.

DefaultAssay(merged.mammal) <- "SCT"

##### Dimensional reduction before integration #####
merged.mammal <- RunPCA(merged.mammal)
merged.mammal <- FindNeighbors(merged.mammal, dims = 1:30, reduction = "pca")
merged.mammal <- FindClusters(merged.mammal, resolution = 2, cluster.name = "unintegrated_clusters")
merged.mammal <- RunUMAP(merged.mammal, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")



##### Integration #####

### RPCA ###
merged.mammal <- IntegrateLayers(object = merged.mammal, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                 normalization.method = "SCT", verbose = TRUE)

### Re-join layers - next step
merged.mammal[["RNA"]] <- JoinLayers(merged.mammal[["RNA"]])


##### Find neighbors and clusters
merged.mammal <- FindNeighbors(merged.mammal, reduction = "integrated.rpca", dims = 1:30)

# Get a wide range of resolutions
merged.mammal <- FindClusters(merged.mammal, resolution = 0.2, cluster.name = "rpca_cluster0.2")
merged.mammal <- FindClusters(merged.mammal, resolution = 0.3, cluster.name = "rpca_cluster0.3")
merged.mammal <- FindClusters(merged.mammal, resolution = 0.8, cluster.name = "rpca_cluster0.8")
merged.mammal <- FindClusters(merged.mammal, resolution = 2, cluster.name = "rpca_cluster2")

## run UMAP on RPCA integration
merged.mammal <- RunUMAP(merged.mammal, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")

##### Add metadata #####

# Create combined meta data
merged.mammal@meta.data$evo_neuronal_celltype_species <- paste0(as.character(merged.mammal@meta.data$evolution), "_", as.character(merged.mammal@meta.data$neuronal_celltype_species))
merged.mammal@meta.data$evo_neuronal_celltype <- paste0(as.character(merged.mammal@meta.data$evolution), "_", as.character(merged.mammal@meta.data$neuronal_celltype))
merged.mammal@meta.data$species_neuronal_celltype_species <- paste0(as.character(merged.mammal@meta.data$species), "_", as.character(merged.mammal@meta.data$neuronal_celltype_species))
merged.mammal@meta.data$species_neuronal_celltype <- paste0(as.character(merged.mammal@meta.data$species), "_", as.character(merged.mammal@meta.data$neuronal_celltype))


##### Save final annotated object #####
merged.mammal@meta.data <- merged.mammal@meta.data %>%
  select(orig.ident, percent.mt, nCount_SCT, nFeature_SCT, sample, species, major_cell_types, neuronal_celltype_species, evolution, vo_neuronal_celltype_species, evo_neuronal_celltype, species_neuronal_celltype_species, species_neuronal_celltype)

# subset to IN and EXN
in.mammal <- subset(merged.mammal, idents = c("Pvalb", "Vip_primate", "Sst", "Lamp5", "Sncg_Pax6", "Pvalb_Vipr2", "Sst_Chodl_mouse"))

exn.mammal <- subset(merged.mammal, idents = c("Pvalb", "Vip_primate", "Sst", "Lamp5", "Sncg_Pax6", "Pvalb_Vipr2", "Sst_Chodl_mouse"), invert = T)

# save cleaned up object
saveRDS(merged.mammal, file = "Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

saveRDS(exn.mammal, file = "Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")

saveRDS(in.mammal, file = "Data_EC_atlas/Seurat_obj/mammals/mammalIN_annot.rds")

