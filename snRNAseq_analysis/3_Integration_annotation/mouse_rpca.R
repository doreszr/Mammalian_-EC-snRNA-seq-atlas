########################################################################################################################

### Author: Dorottya Ralbovszki

### Date: 09-24

### Normalization and batch integration of published mouse samples using SCTransform and RPCA respectively

### Label transfer from Allen EC snRNA-seq neuronal only reference

### Neuronal subtype annotation by reference (Allen EC snRNA-seq data set) and marker genes

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
mouse1 <- readRDS("Data_EC_atlas/Seurat_obj/mouse/mouse1.rds")
mouse2 <- readRDS("Data_EC_atlas/Seurat_obj/mouse/mouse2.rds")
mouse3 <- readRDS("Data_EC_atlas/Seurat_obj/mouse/mouse3.rds")

##### Add metadata #####
mouse1$sample <- "MOUSE1"
mouse2$sample <- "MOUSE2"
mouse3$sample <- "MOUSE3"

mouse1$species <- "MOUSE"
mouse2$species <- "MOUSE"
mouse3$species <- "MOUSE"

mouse1$evolution <- "MOUSE"
mouse2$evolution <- "MOUSE"
mouse3$evolution <- "MOUSE"

##### Merging objects #####

# merge for integrative layers
merged.mouse <- merge(mouse1, y = c(mouse2, mouse3), add.cell.ids = c("mouse1", "mouse2", "mouse3"), project = "MOUSE_EC")

table(merged.mouse$orig.ident)
Idents(merged.mouse)

### Rerunning normalization with SCTransform ###
merged.mouse <- SCTransform(merged.mouse, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.

DefaultAssay(merged.mouse) <- "SCT"

### Dimensional reduction before integration ###
merged.mouse <- RunPCA(merged.mouse)
merged.mouse <- FindNeighbors(merged.mouse, dims = 1:30, reduction = "pca")
merged.mouse <- FindClusters(merged.mouse, resolution = 2, cluster.name = "unintegrated_clusters")

merged.mouse <- RunUMAP(merged.mouse, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


##### RPCA Integration #####
merged.mouse <- IntegrateLayers(object = merged.mouse, method = RPCAIntegration, orig.reduction = "pca",
                                new.reduction = "integrated.rpca", normalization.method = "SCT", verbose = T)

# Join layers after integration
merged.mouse[["RNA"]] <- JoinLayers(merged.mouse[["RNA"]])

##### Find neighbors and clusters ###
merged.mouse <- FindNeighbors(merged.mouse, reduction = "integrated.rpca", dims = 1:30)

# We can get a wide range of resolutions and they will all be stored in the
# object metadata.

merged.mouse <- FindClusters(merged.mouse,resolution = 0.1, cluster.name = "rpca_cluster0.1", verbose = T)
merged.mouse <- FindClusters(merged.mouse,resolution = 0.2, cluster.name = "rpca_cluster0.2", verbose = T)
merged.mouse <- FindClusters(merged.mouse,resolution = 0.3, cluster.name = "rpca_cluster0.3", verbose = T)
merged.mouse <- FindClusters(merged.mouse,resolution = 0.5, cluster.name = "rpca_cluster0.5", verbose = T)
merged.mouse <- FindClusters(merged.mouse,resolution = 0.8, cluster.name = "rpca_cluster0.8", verbose = T)
merged.mouse <- FindClusters(merged.mouse,resolution = 1.0, cluster.name = "rpca_cluster1.0", verbose = T)


### run UMAP on RPCA integration ###
merged.mouse <- RunUMAP(merged.mouse, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


##### Annotate with mouse atlas #####

# open mouse atlas neuronal cell types only
mouse_atlas.neu <- readRDS("Data_EC_atlas/Seurat_obj/mouseAllen/mouse_atlasv5_neu.rds")


# Make sure your reference and query are in the SCT assay
DefaultAssay(merged.mouse) <- "SCT"
DefaultAssay(mouse_atlas.neu) <- "SCT"

### Find transfer anchors ###
re.anchors <- FindTransferAnchors(reference = mouse_atlas.neu, 
                                  query = merged.mouse,
                                  normalization.method = 'SCT',
                                  dims = 1:30, reference.reduction = "pca")

# below is the workflow for generating prediction scores for each cell in your query
# object. I have done this 3 times as I want my query object you have 3 sets of cell type ids.

# metadata #1 : class_label
predictions.class <- TransferData(anchorset = re.anchors, 
                                  refdata =  mouse_atlas.neu$class_label,
                                  dims = 1:30)

merged.mouse <- AddMetaData(merged.mouse, metadata = predictions.class)


merged.mouse$class_label <- merged.mouse$predicted.id



# metadata #2 : subclass_label
predictions.subclass <- TransferData(anchorset = re.anchors, 
                                     refdata = mouse_atlas.neu$subclass_label,
                                     dims = 1:30)

merged.mouse <- AddMetaData(merged.mouse, metadata = predictions.subclass)



merged.mouse$subclass_label <- merged.mouse$predicted.id



# metadata #3 : cluster_label
predictions.cluster <- TransferData(anchorset = re.anchors, 
                                    refdata = mouse_atlas.neu$cluster_label,
                                    dims = 1:30)

merged.mouse <- AddMetaData(merged.mouse, metadata = predictions.cluster)


merged.mouse$cluster_label <- merged.mouse$predicted.id


# your query object will now have three new metadata slots containing predicted cell
# types at 3 resolutions.

## Save stripped down object for further integration

# switch to RNA assay for stripping data set down
DefaultAssay(merged.mouse) <- "RNA"

# Strip down your data set
mouse.diet <- DietSeurat(merged.mouse, assays = "RNA", dimreducs = NULL)
mouse.diet[["RNA"]] <- split(mouse.diet[["RNA"]], f = mouse.diet$orig.ident)

# Save your stripped data set
saveRDS(mouse.diet, file = "Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot_diet.rds")



##### Subset and save final annotated neuronal data set #####

merged.mouse@meta.data <- merged.mouse@meta.data %>%
  select(orig.ident, percent.mt, nCount_SCT, nFeature_SCT, sample, species, evolution, neuronal_celltype_species, class_label)


# Subset IN and EXN neurons
in.mouse <- subset(merged.mouse, idents = c("Pvalb", "Vip", "Sst", "Lamp5", "Lamp5_Lhx6", "Sncg", "Pvalb_Vipr2", "Sst_Chodl"))

exn.mouse <- subset(merged.mouse, idents = c("Pvalb", "Vip", "Sst", "Lamp5", "Lamp5_Lhx6", "Sncg", "Pvalb_Vipr2", "Sst_Chodl"), invert = T)


### Save your new annotated and subset objects
saveRDS(merged.mouse, file = "Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")

saveRDS(in.mouse, file = "Data_EC_atlas/Seurat_obj/mouse/merged.mouseIN_annot.rds")

saveRDS(exn.mouse, file = "Data_EC_atlas/Seurat_obj/mouse/merged.mouseEXN_annot.rds")


