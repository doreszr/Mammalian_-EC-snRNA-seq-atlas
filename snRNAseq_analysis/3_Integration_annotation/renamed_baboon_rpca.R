########################################################################################################################

### Author: Dorottya Ralbovszki
### Date: 11-24

### Batch integration and normalization of baboon (renamed to human orthologues) samples using SCTransform and RPCA respectively

### Major cell type annotation

### Sub-setting neuronal cells

### saving Seurat object for further analysis

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
library(DoubletFinder)
library(reticulate)
library(presto)
library(dplyr)



##### Load preprocessed samples #####
baboon1 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamedbaboon1.rds")
baboon2 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamedbaboon2.rds")
baboon3 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamedbaboon3.rds")


##### Add metadata #####
baboon1$sample <- "BABOON1"
baboon2$sample <- "BABOON2"
baboon3$sample <- "BABOON3"


baboon1$species <- "BABOON"
baboon2$species <- "BABOON"
baboon3$species <- "BABOON"

baboon1$evolution <- "PRIMATE"
baboon2$evolution <- "PRIMATE"
baboon3$evolution <- "PRIMATE"


##### Merging objects #####

# merge for integrative layers
merged.baboon <- merge(baboon1, y = c(baboon2, baboon3), add.cell.ids = c("baboon1", "baboon2", "baboon3"), project = "BABOON_EC")

table(merged.baboon$orig.ident)
Idents(merged.baboon)

##### Rerunning normalization with SCTransform #####
merged.baboon <- SCTransform(merged.baboon, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.

DefaultAssay(merged.baboon) <- "SCT"


##### Dimensional reduction before integration #####
merged.baboon <- RunPCA(merged.baboon)
merged.baboon <- FindNeighbors(merged.baboon, dims = 1:30, reduction = "pca")
merged.baboon <- FindClusters(merged.baboon, resolution = 2, cluster.name = "unintegrated_clusters")

# Run UMAP before integration
merged.baboon <- RunUMAP(merged.baboon, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


##### RPCA Integration #####

merged.baboon <- IntegrateLayers(object = merged.baboon, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                 normalization.method = "SCT", verbose = T)

merged.baboon[["RNA"]] <- JoinLayers(merged.baboon[["RNA"]])

##### Find neighbors and clusters

merged.baboon <- FindNeighbors(merged.baboon, reduction = "integrated.rpca", dims = 1:30)

# We can get a wide range of resolutions and they will all be stored in the
# object metadata.

merged.baboon <- FindClusters(merged.baboon,resolution = 0.1, verbose = T, cluster.name = "rpca_cluster0.1")
merged.baboon <- FindClusters(merged.baboon,resolution = 0.2, verbose = T, cluster.name = "rpca_cluster0.2")
merged.baboon <- FindClusters(merged.baboon,resolution = 0.3, cluster.name = "rpca_cluster0.3", verbose = T)


##### Run UMAP on RPCA integration

merged.baboon <- RunUMAP(merged.baboon, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


##### Cell type annotation #####

# make sure you are in the right assay
DefaultAssay(merged.baboon) <- "SCT"

##### Canonical marker genes

#IN
FeaturePlot(merged.baboon, features = "GAD1", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "PVALB", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "SST", reduction = "umap.rpca")



FeaturePlot(merged.baboon, features = "ADARB2", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "VIP", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "LHX6", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "LAMP5", reduction = "umap.rpca")


# EXN
FeaturePlot(merged.baboon, features = "SLC17A7", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "CUX2", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "RORB", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "THEMIS", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "FEZF2", reduction = "umap.rpca")

# ASTR
FeaturePlot(merged.baboon , features = "AQP4", reduction = "umap.rpca")
FeaturePlot(merged.baboon , features = "AQP1", reduction = "umap.rpca")
FeaturePlot(merged.baboon , features = "FGFR3", reduction = "umap.rpca")
FeaturePlot(merged.baboon , features = "SERPINI2", reduction = "umap.rpca")
FeaturePlot(merged.baboon , features = "PLCG1", reduction = "umap.rpca")


# OPC
FeaturePlot(merged.baboon , features = "COL20A1", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "ADARB2", reduction = "umap.rpca")
FeaturePlot(merged.baboon, features = "PDGFRA", reduction = "umap.rpca")

FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "OLIG1")
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "OLIG2")
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "CSPG4")
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "PLP1")
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "C1QB") # microglia
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "CX3CR1") # microglia https://doi.org/10.1038/s43587-023-00424-y
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "VWF") # endothelial cells 
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "FLT1") # endothelial cells doi:10.1523/JNEUROSCI.0237-23.2023
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "MECOM") # endothelial cells doi:10.1523/JNEUROSCI.0237-23.2023
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "CLDN5") # endothelial cells doi:10.1523/JNEUROSCI.0237-23.2023
FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "ABCB1") # endothelial cells doi:10.1523/JNEUROSCI.0237-23.2023


FeaturePlot(merged.baboon, features = "CALB1")
FeaturePlot(merged.baboon, features = "APOE", reduction = "umap.rpca")

FeaturePlot(merged.baboon, reduction = "umap.rpca", features = "OPALIN")


##### Rename clusters to major cell types #####

# make sure you have the right active ident
Idents(merged.baboon) <- merged.baboon$rpca_cluster0.1

## rename clusters
merged.baboon <- RenameIdents(merged.baboon, "0" = "EXN",  "1" = "EXN", "2" = "IN", "3" = "EXN", "4" = "ASTR", "5" = "IN",
                              "6" = "EXN", "7" = "IN", "8" = "EXN", "9" = "EXN", "10" = "OPC", "11" = "Oligo", "12" = "EXN",
                              "13" = "Endothelial", "14" = "IN")


## if all good, save major cell types as new meta data
merged.baboon$major_cell_types <- Idents(merged.baboon)


##### Save Seurat objects for further analyses #####

##### Save your new annotated object
saveRDS(merged.baboon, file = "Data_EC_atlas/Seurat_obj/baboon/renamed.merged.baboon_annot.rds")


### Save stripped down object for further integration

# switch to RNA assay for stripping data set down
DefaultAssay(merged.baboon) <- "RNA"

# Strip down your data set
merged.baboon.diet <- DietSeurat(merged.baboon, assays = "RNA", dimreducs = NULL)
merged.baboon.diet[["RNA"]] <- split(merged.baboon.diet[["RNA"]], f = merged.baboon.diet$orig.ident)


# Save your stripped data set
saveRDS(merged.baboon.diet, file = "Data_EC_atlas/Seurat_obj/baboon/renamed.merged.baboon_annot_diet.rds")


##### Neuronal cells only #####

##### Subset neuronal cell types
neu.merged.baboon <- subset(merged.baboon, idents = c("ASTR", "Oligo", "Endothelial", "OPC"), invert = TRUE)

# switch to RNA assay for stripping data set down
DefaultAssay(neu.merged.baboon) <- "RNA"

# Strip down neuronal subset
merged.baboon.neu <- DietSeurat(neu.merged.baboon, assays = "RNA", dimreducs = NULL)
merged.baboon.neu[["RNA"]] <- split(merged.baboon.neu[["RNA"]], f = merged.baboon.neu$orig.ident)



##### Save your neuronal object

saveRDS(merged.baboon.neu, file = "Data_EC_atlas/Seurat_obj/baboon/renamed.merged.baboon_neu_diet.rds")

