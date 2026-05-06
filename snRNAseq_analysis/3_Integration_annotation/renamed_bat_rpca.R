########################################################################################################################

### Author: Dorottya Ralbovszki
### Date: 11-24

### Normalisation and batch integration of fruit bat (renamed to mouse orthologues) samples using SCTransform and RPCA respectively

### Label transfer from Allen EC snRNA-seq reference

### Major cell type annotation by reference (Allen EC snRNA-seq data set) and marker genes

### Sub-setting neuronal cells

### Label transfer from Allen EC snRNA-seq neuronal only reference

### Annotating neuronal subtypes by reference (Allen EC snRNA-seq data set) and marker genes

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



##### Load preprocessed samples #####
bat3 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedbat3.rds")
bat4 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedbat4.rds")
bat5 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedbat5.rds")

##### Add metadata #####
bat3$sample <- "BAT3"
bat4$sample <- "BAT4"
bat5$sample <- "BAT5"

bat3$species <- "BAT"
bat4$species <- "BAT"
bat5$species <- "BAT"

bat3$evolution <- "BAT"
bat4$evolution <- "BAT"
bat5$evolution <- "BAT"


##### Merging objects #####

# merge for integrative layers
merged.bat <- merge(bat3, y = c(bat4, bat5), add.cell.ids = c("bat3", "bat4", "bat5"), project = "BAT_EC")

table(merged.bat$orig.ident)
Idents(merged.bat)

##### Rerunning normalization with SCTransform #####
merged.bat <- SCTransform(merged.bat, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.
DefaultAssay(merged.bat) <- "SCT"


##### Dimensional reduction before integration #####
merged.bat <- RunPCA(merged.bat)
merged.bat <- FindNeighbors(merged.bat, dims = 1:30, reduction = "pca")
merged.bat <- FindClusters(merged.bat, resolution = 2, cluster.name = "unintegrated_clusters")

merged.bat <- RunUMAP(merged.bat, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


##### RPCA Integration #####
merged.bat <- IntegrateLayers(object = merged.bat, method = RPCAIntegration, orig.reduction = "pca",
                              new.reduction = "integrated.rpca", normalization.method = "SCT", verbose = T)

# Joining layers
merged.bat[["RNA"]] <- JoinLayers(merged.bat[["RNA"]])

##### Find neighbors and clusters
merged.bat <- FindNeighbors(merged.bat, reduction = "integrated.rpca", dims = 1:30)

# We can get a wide range of resolutions and they will all be stored in the
# object metadata.

merged.bat <- FindClusters(merged.bat,resolution = 0.1, verbose = T, cluster.name = "rpca_cluster0.1")
merged.bat <- FindClusters(merged.bat,resolution = 0.2, verbose = T, cluster.name = "rpca_cluster0.2")
merged.bat <- FindClusters(merged.bat,resolution = 0.3, verbose = T, cluster.name = "rpca_cluster0.3")

## run UMAP on RPCA integration
merged.bat <- RunUMAP(merged.bat, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


##### Annotate with mouse atlas #####

### Load Allen mouse brain atlas ###

mouse_atlas <- readRDS("Data_EC_atlas/Seurat_obj/mouseAllen/mouse_atlasv5.rds")


# Make sure your reference and query are in the SCT assay
DefaultAssay(merged.bat) <- "SCT"
DefaultAssay(mouse_atlas) <- "SCT"


# Find transfer anchors
re.anchors <- FindTransferAnchors(reference = mouse_atlas, 
                                  query = merged.bat,
                                  normalization.method = 'SCT',
                                  dims = 1:30, reference.reduction = "pca")

# below is the workflow for generating prediction scores for each cell in your query
# object. I have done this 3 times as I want my query object you have 3 sets of cell type ids.

# metadata #1 : class_label
predictions.class <- TransferData(anchorset = re.anchors, 
                                  refdata =  mouse_atlas$class_label,
                                  dims = 1:30)


merged.bat <- AddMetaData(merged.bat, metadata = predictions.class)

merged.bat$class_label <- merged.bat$predicted.id


# metadata #2 : subclass_label
predictions.subclass <- TransferData(anchorset = re.anchors, 
                                     refdata = mouse_atlas$subclass_label,
                                     dims = 1:30)


merged.bat <- AddMetaData(merged.bat, metadata = predictions.subclass)

merged.bat$subclass_label <- merged.bat$predicted.id


# metadata #3 : cluster_label
predictions.cluster <- TransferData(anchorset = re.anchors, 
                                    refdata = mouse_atlas$cluster_label,
                                    dims = 1:30)


merged.bat <- AddMetaData(merged.bat, metadata = predictions.cluster)

merged.bat$cluster_label <- merged.bat$predicted.id

# your query object will now have three new metadata slots containing predicted cell
# types at 3 resolutions.


# Canonical markers

# IN
FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Gad1")
FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Sst")
FeaturePlot(merged.bat, features = "Gad2", reduction = "umap.rpca") # 6, 8, 11
FeaturePlot(merged.bat, features = c("Gad1", "Gad2", "Sst"), reduction = "umap.rpca")

# EXN
FeaturePlot(merged.bat, features = "Slc17a7", reduction = "umap.rpca") # 0, 1, 2, 3, 4, 5, 7, 9, 12, 13, 15
FeaturePlot(merged.bat, features = "Cux2", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = "Rorb", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = "Fezf2", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = c("Slc17a7", "Cux2"), reduction = "umap.rpca")

# Oligo
FeaturePlot(merged.bat, features = "Mbp", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = "Plp1", reduction = "umap.rpca")
FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Opalin") # 14 star part

# OPC
FeaturePlot(merged.bat, features = "Pdgfra", reduction = "umap.rpca")


Idents(merged.bat) <- merged.bat$rpca_cluster0.2
merged.bat <- PrepSCTFindMarkers(merged.bat)

clust10 <- FindMarkers(merged.bat, ident.1 = 10)

# ASTR/VLMC (clust 10)
FeaturePlot(merged.bat, features = "Aqp4", reduction = "umap.rpca")
FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Apoe")
FeaturePlot(merged.bat, features = "Ptgds", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = "Gja1", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = "Selenop", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = "Sparc", reduction = "umap.rpca")
FeaturePlot(merged.bat, features = "Slc15a2", reduction = "umap.rpca") 



##### Rename clusters to major cell types #####

# make sure you have the right active ident
Idents(merged.bat) <- merged.bat$rpca_cluster0.2
merged.bat <- RenameIdents(merged.bat, "0" = "EXN",  "1" = "EXN", "2" = "EXN", "3" = "EXN", "4" = "EXN", "5" = "EXN",
                           "6" = "IN", "7" = "EXN", "8" = "IN", "9" = "EXN", "10" = "VLMC/Astrocyte", "11" = "IN", "12" = "EXN",
                           "13" = "EXN", "14" = "OPC", "15" = "EXN")


merged.bat$major_cell_types <- Idents(merged.bat)


###### Save your annotated clustered object ######
saveRDS(merged.bat, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_annot.rds")




##### Neuronal cell types #####

##### Subset neuronal cell types
Idents(merged.bat) <- merged.bat$major_cell_types

neu.merged.bat <- subset(merged.bat, idents = c("Oligo", "VLMC/Astrocyte/Endothelial", "OPC"), invert = T)


##### Transform object for further annotation #####

# Keep only RNA assay and meta data
# Split the layers again

# Make sure your object is in the RNA assay
DefaultAssay(neu.merged.bat) <- "RNA"


# Keep only the RNA assay with the meta data
merged.bat.neu <- DietSeurat(neu.merged.bat, assays = "RNA", dimreducs = NULL)

# Split  layers for integration
merged.bat.neu[["RNA"]] <- split(merged.bat.neu[["RNA"]], f = merged.bat.neu$orig.ident)


##### Rerunning normalization with SCTransform #####
merged.bat.neu <- SCTransform(merged.bat.neu, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.
DefaultAssay(merged.bat.neu) <- "SCT"


##### Dimensional reduction before integration #####
merged.bat.neu <- RunPCA(merged.bat.neu)
merged.bat.neu <- FindNeighbors(merged.bat.neu, dims = 1:30, reduction = "pca")
merged.bat.neu <- FindClusters(merged.bat.neu, resolution = 2, cluster.name = "unintegrated_clusters")

merged.bat.neu <- RunUMAP(merged.bat.neu, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


##### RPCA Integration #####
merged.bat.neu <- IntegrateLayers(object = merged.bat.neu, method = RPCAIntegration, orig.reduction = "pca",
                                  new.reduction = "integrated.rpca", normalization.method = "SCT", verbose = T)

# Joining layers
merged.bat.neu[["RNA"]] <- JoinLayers(merged.bat.neu[["RNA"]])

##### Find neighbors and clusters
merged.bat.neu <- FindNeighbors(merged.bat.neu, reduction = "integrated.rpca", dims = 1:30)

# We can get a wide range of resolutions and they will all be stored in the
# object metadata.

merged.bat.neu <- FindClusters(merged.bat.neu,resolution = 0.1, verbose = T, cluster.name = "rpca_cluster0.1")
merged.bat.neu <- FindClusters(merged.bat.neu,resolution = 0.2, verbose = T, cluster.name = "rpca_cluster0.2")
merged.bat.neu <- FindClusters(merged.bat.neu,resolution = 0.3, verbose = T, cluster.name = "rpca_cluster0.3")
merged.bat.neu <- FindClusters(merged.bat.neu,resolution = 0.5, cluster.name = "rpca_cluster0.5", verbose = T)
merged.bat.neu <- FindClusters(merged.bat.neu,resolution = 0.8, cluster.name = "rpca_cluster0.8", verbose = T)
merged.bat.neu <- FindClusters(merged.bat.neu,resolution = 1.5, cluster.name = "rpca_cluster1.5", verbose = T)

## run UMAP on RPCA integration
merged.bat.neu <- RunUMAP(merged.bat.neu, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


##### Annotation to Allen mouse brain atlas EC reference neuronal only #####

# import reference.
mouse_atlas.neu <- readRDS("Data_EC_atlas/Seurat_obj/mouseAllen/mouse_atlasv5_neu.rds")

#In this example, the qyuery seurat object is named re
# Make sure your reference and query are in the SCT assay
DefaultAssay(merged.bat.neu) <- "SCT"
DefaultAssay(mouse_atlas.neu) <- "SCT"

# Find transfer anchors
re.anchors <- FindTransferAnchors(reference = mouse_atlas.neu, 
                                  query = merged.bat.neu,
                                  normalization.method = 'SCT',
                                  dims = 1:30)

# below is the workflow for generating prediction scores for each cell in your query
# object. I have done this 3 times as I want my query object you have 3 sets of cell type ids.

# metadata #1 : class_label
predictions.class <- TransferData(anchorset = re.anchors, 
                                  refdata =  mouse_atlas.neu$class_label,
                                  dims = 1:30)

merged.bat.neu <- AddMetaData(merged.bat.neu, metadata = predictions.class)

merged.bat.neu$class_label <- merged.bat.neu$predicted.id

# metadata #2 : subclass_label
predictions.subclass <- TransferData(anchorset = re.anchors, 
                                     refdata = mouse_atlas.neu$subclass_label,
                                     dims = 1:30)

merged.bat.neu <- AddMetaData(merged.bat.neu, metadata = predictions.subclass)

merged.bat.neu$subclass_label <- merged.bat.neu$predicted.id

merged.bat <- merged.bat.neu


##### Transform object for further annotation #####

# Keep only RNA assay and meta data
# Split the layers again

# Make sure your object is in the RNA assay
DefaultAssay(merged.bat.neu) <- "RNA"


# Keep only the RNA assay with the meta data
merged.bat.neu <- DietSeurat(merged.bat.neu, assays = "RNA", dimreducs = NULL)

# Split  layers for integration
merged.bat.neu[["RNA"]] <- split(merged.bat.neu[["RNA"]], f = merged.bat.neu$orig.ident)

# Save your new subset object
saveRDS(merged.bat.neu, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_annot_diet.rds")


##### Subset and save final annotated neuronal data set #####

merged.bat@meta.data <- merged.bat@meta.data %>%
  select(orig.ident, percent.mt, nCount_SCT, nFeature_SCT, sample, species, neuronal_celltype_species, evolution)



in.bat <- subset(merged.bat, idents = c("Pvalb", "Vip/Sncg", "Sst", "Lamp5", "Pvalb_Vipr2", "Lamp5_Lhx6"))

exn.bat <- subset(merged.bat, idents = c("Pvalb", "Vip/Sncg", "Sst", "Lamp5", "Pvalb_Vipr2", "Lamp5_Lhx6"), invert = T)

# Save your new subset object
saveRDS(merged.bat, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_annot.rds")

saveRDS(in.bat, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_IN_annot.rds")

saveRDS(exn.bat, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_EXN_annot.rds")

