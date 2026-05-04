########################################################################################################################

### Author: Dorottya Ralbovszki
### Date: 11-24

### Normalization and batch integration of published human samples using SCTransform and RPCA respectively

### Major cell type annotation

### Sub-setting neuronal cells

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
library(reticulate)
library(presto)
library(dplyr)


##### Load preprocessed samples #####
human1 <- readRDS("Data_EC_atlas/Seurat_obj/human/human1.rds")
human2 <- readRDS("Data_EC_atlas/Seurat_obj/human/human2.rds")
human3 <- readRDS("Data_EC_atlas/Seurat_obj/human/human3.rds")


##### Add metadata #####
human1$sample <- "HUMAN1"
human2$sample <- "HUMAN2"
human3$sample <- "HUMAN3"

human1$species <- "HUMAN"
human2$species <- "HUMAN"
human3$species <- "HUMAN"

human1$evolution <- "PRIMATE"
human2$evolution <- "PRIMATE"
human3$evolution <- "PRIMATE"



##### Merging objects #####

# merge for integrative layers
merged.human <- merge(human1, y = c(human2, human3), add.cell.ids = c("human1", "human2", "human3"), project = "HUMAN_EC")


##### Rerunning normalization with SCTransform #####
merged.human <- SCTransform(merged.human, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.

DefaultAssay(merged.human) <- "SCT"


##### Dimensional reduction before integration #####
merged.human <- RunPCA(merged.human)
merged.human <- FindNeighbors(merged.human, dims = 1:30, reduction = "pca")
merged.human <- FindClusters(merged.human, resolution = 2, cluster.name = "unintegrated_clusters")

# Run UMAP before integration
merged.human <- RunUMAP(merged.human, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


##### RPCA Integration #####

merged.human <- IntegrateLayers(object = merged.human, method = RPCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.rpca", normalization.method = "SCT", verbose = T)

merged.human[["RNA"]] <- JoinLayers(merged.human[["RNA"]])

##### Find neighbors and clusters

merged.human <- FindNeighbors(merged.human, reduction = "integrated.rpca", dims = 1:30)

# We can get a wide range of resolutions and they will all be stored in the
# object metadata.

merged.human <- FindClusters(merged.human,resolution = 0.1, verbose = T, cluster.name = "rpca_cluster0.1")
merged.human <- FindClusters(merged.human,resolution = 0.2, verbose = T, cluster.name = "rpca_cluster0.2")
merged.human <- FindClusters(merged.human,resolution = 0.3, cluster.name = "rpca_cluster0.3", verbose = T)

##### Run UMAP on RPCA integration

merged.human <- RunUMAP(merged.human, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")



##### Cell type annotation #####

# make sure you are in the right assay
DefaultAssay(merged.human) <- "SCT"

DimPlot(merged.human, reduction = "umap.rpca", label = TRUE, label.size = 10, group.by = "rpca_cluster0.1") + ggtitle("Human RPCA")

##### Canonical marker genes
FeaturePlot(merged.human, reduction = "umap.rpca", features = "GAD1")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "AQP4")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "AQP1")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "OPALIN")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "ADARB2")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "SST")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "SLC17A7")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "PDGFRA")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "OLIG1")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "OLIG2")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "CSPG4")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "PLP1")
FeaturePlot(merged.human, reduction = "umap.rpca", features = "C1QB") # microglia
FeaturePlot(merged.human, reduction = "umap.rpca", features = "CX3CR1") # microglia https://doi.org/10.1038/s43587-023-00424-y
FeaturePlot(merged.human, reduction = "umap.rpca", features = "VWF") # endothelial cells 
FeaturePlot(merged.human, reduction = "umap.rpca", features = "FLT1") # endothelial cells doi:10.1523/JNEUROSCI.0237-23.2023 
FeaturePlot(merged.human, reduction = "umap.rpca", features = "APOE")


##### Rename clusters to major cell types #####

# make sure you have the right active ident
Idents(merged.human) <- merged.human$rpca_cluster0.2

## rename clusters
merged.human <- RenameIdents(merged.human, "0" = "EXN",  "1" = "ASTR", "2" = "OPC", "3" = "EXN", "4" = "IN", "5" = "Oligo",
                             "6" = "IN", "7" = "EXN", "8" = "EXN", "9" = "IN", "10" = "EXN", "11" = "EXN", "12" = "Microglia",
                             "13" = "Oligo/ASTR", "14" = "Endothelial", "15" = "EXN", "16" = "IN", "17" = "Oligo", "18" = "EXN")



## if all good, save major cell types as new meta data
merged.human$major_cell_types <- Idents(merged.human)



##### Save Seurat objects for further analyses #####

### Save your annotated object
saveRDS(merged.human, file = "Data_EC_atlas/Seurat_obj/human/merged.human_annot.rds")


### Save stripped down object for further integration

# switch to RNA assay for stripping data set down
DefaultAssay(merged.human) <- "RNA"

# Strip down your data set
merged.human.diet <- DietSeurat(merged.human, assays = "RNA", dimreducs = NULL)
merged.human.diet[["RNA"]] <- split(merged.human.diet[["RNA"]], f = merged.human.diet$orig.ident)


# Save your stripped data set
saveRDS(merged.human.diet, file = "Data_EC_atlas/Seurat_obj/human/merged.human_annot_diet.rds")



##### Neuronal cell types #####

##### Subset neuronal cell types
Idents(merged.human) <- merged.human$major_cell_types
neu.merged.human <- subset(merged.human, idents = c("ASTR", "Oligo/ASTR", "OPC", "Oligo", "Endothelial", "Microglia"), invert = TRUE)

# check number of cells of each cell type
table(neu.merged.human$major_cell_types)

# switch to RNA assay for stripping data set down
DefaultAssay(neu.merged.human) <- "RNA"

# Strip down neuronal subset
merged.human.neu <- DietSeurat(neu.merged.human, assays = "RNA", dimreducs = NULL)
merged.human.neu[["RNA"]] <- split(merged.human.neu[["RNA"]], f = merged.human.neu$orig.ident)


# Save your stripped neuronal subset
saveRDS(merged.human.neu, file = "Data_EC_atlas/Seurat_obj/human/merged.human_neu_diet.rds")

