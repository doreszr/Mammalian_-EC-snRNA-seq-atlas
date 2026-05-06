########################################################################################################################

### Author: Dorottya Ralbovszki

### Date: 11-24

### Cross-species integration of baboon (renamed to human orthologues) and human samples using RPCA and other integration methods

### Renaming human and baboon (renamed to human orthologues) genes to their mouse orthologues in neuronal only data sets

### Cross-species integration of baboon and human (renamed to mouse orthologues) neuronal cells using RPCA

### Label transfer from Allen EC snRNA-seq neuronal only reference

### Neuronal subtype annotation by Allen EC snRNA-seq neuronal only reference and marker genes

### Saving Seurat objects for further analysis

########################################################################################################################


library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(cowplot)
library(sctransform)
library(glmGamPoi)
library(ROCR)
library(parallel)
library(fields)
library(presto)
library(SeuratWrappers)

##### Load preprocessed and cell typed stripped down samples #####
human <- readRDS("Data_EC_atlas/Seurat_obj/human/merged.human_annot_diet.rds")

baboon <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamed.merged.baboon_annot_diet.rds")

##### Merging objects #####

# merge for integrative layers
merged.primate <- merge(x = human, y = baboon, project = "PRIMATE_EC")


##### Rerunning normalization with SCTransform #####
merged.primate <- SCTransform(merged.primate, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.
DefaultAssay(merged.primate) <- "SCT"

##### Dimensional reduction before integration #####
merged.primate <- RunPCA(merged.primate)
merged.primate <- FindNeighbors(merged.primate, dims = 1:30, reduction = "pca")
merged.primate <- FindClusters(merged.primate, resolution = 2, cluster.name = "unintegrated_clusters")

merged.primate <- RunUMAP(merged.primate, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


##### Integration #####

### RPCA ###
merged.primate <- IntegrateLayers(object = merged.primate, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                                  normalization.method = "SCT", verbose = TRUE)

##### Find neighbors and clusters
merged.primate <- FindNeighbors(merged.primate, reduction = "integrated.rpca", dims = 1:30)

merged.primate <- FindClusters(merged.primate, resolution = 0.1, cluster.name = "rpca_cluster0.1")

merged.primate <- FindClusters(merged.primate, resolution = 0.2, cluster.name = "rpca_cluster0.2")

merged.primate <- FindClusters(merged.primate, resolution = 0.3, cluster.name = "rpca_cluster0.3")

## run UMAP on RPCA integration
merged.primate <- RunUMAP(merged.primate, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


### CCA ###
merged.primate <- IntegrateLayers(object = merged.primate, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method = "SCT", verbose = T)


##### Find neighbors and clusters

merged.primate <- FindNeighbors(merged.primate, reduction = "integrated.cca", dims = 1:30)

# We can get a wide range of resolutions and they will all be stored in the
# object metadata.

merged.primate <- FindClusters(merged.primate, resolution = 0.6, cluster.name = "cca_clusters", verbose = T)

## run UMAP on CCA integration

merged.primate <- RunUMAP(merged.primate, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")


### Harmony ###
merged.primate <- IntegrateLayers(object = merged.primate, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", assay = "SCT", verbose = FALSE)

##### Find neighbors and clusters

merged.primate <- FindNeighbors(merged.primate, reduction = "harmony", dims = 1:30)

# use the same resolution as in CCA for better comparison

merged.primate <- FindClusters(merged.primate, resolution = 0.6, cluster.name = "harmony_cluster")

## run UMAP on Harmony integration

merged.primate <- RunUMAP(merged.primate, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

# Joining layers again
merged.primate[["RNA"]] <- JoinLayers(merged.primate[["RNA"]])


###### Save Seurat objects for further analyses ######

### Save your new clustered object
saveRDS(merged.primate, file = "merged.primate_int.rds")


##### Renaming neuronal to mouse #####

### Load neuronal only samples 
human <- readRDS("merged.human_neu_diet.rds")

baboon <- readRDS("renamedmerged.baboon_neudiet.rds")


### Subset each sample for renaming
Idents(baboon) <- baboon$orig.ident
Idents(human) <- human$orig.ident

baboon1 <- subset(x = baboon, idents = "baboon1")
baboon2 <- subset(x = baboon, idents = "baboon2")
baboon3 <- subset(x = baboon, idents = "baboon3")
human1 <- subset(x = human, idents = "human1")
human2 <- subset(x = human, idents = "human2")
human3 <- subset(x = human, idents = "human3")


### Rename human gene names to mouse

# Load ortho df mouse vs human
allorthohm <- read.csv("Data_EC_atlas/orthologues/mappingorthologs/human_vs_mouse/allorthohm_primate.csv", sep = ";")


# Create a count matrix from RNA assaq
exp_mtx.baboon1 <- as.matrix(baboon1[["RNA"]]$counts)
exp_mtx.baboon2 <- as.matrix(baboon2[["RNA"]]$counts)
exp_mtx.baboon3 <- as.matrix(baboon3[["RNA"]]$counts)
exp_mtx.human1 <- as.matrix(human1[["RNA"]]$counts)
exp_mtx.human2 <- as.matrix(human2[["RNA"]]$counts)
exp_mtx.human3 <- as.matrix(human3[["RNA"]]$counts)

# get the genes that are present in the count matrix
allorthohmb1 <- filter(allorthohm,
                       gene_h %in% rownames(exp_mtx.baboon1))
allorthohmb2 <- filter(allorthohm,
                       gene_h %in% rownames(exp_mtx.baboon2))
allorthohmb3 <- filter(allorthohm,
                       gene_h %in% rownames(exp_mtx.baboon3))
allorthohmh1 <- filter(allorthohm,
                       gene_h %in% rownames(exp_mtx.human1))
allorthohmh2 <- filter(allorthohm,
                       gene_h %in% rownames(exp_mtx.human2))
allorthohmh3 <- filter(allorthohm,
                       gene_h %in% rownames(exp_mtx.human3))

## Filter the expression matrix for genes which a human counterpart is available
exp_mtx.baboon1 <- exp_mtx.baboon1[allorthohmb1$gene_h,]
exp_mtx.baboon2 <- exp_mtx.baboon2[allorthohmb2$gene_h,]
exp_mtx.baboon3 <- exp_mtx.baboon3[allorthohmb3$gene_h,]
exp_mtx.human1 <- exp_mtx.human1[allorthohmh1$gene_h,]
exp_mtx.human2 <- exp_mtx.human2[allorthohmh2$gene_h,]
exp_mtx.human3 <- exp_mtx.human3[allorthohmh3$gene_h,]


## Now chnage the rownames of the matrix to the human gene names
rownames(exp_mtx.baboon1) <- allorthohmb1$gene_m
rownames(exp_mtx.baboon2) <- allorthohmb2$gene_m
rownames(exp_mtx.baboon3) <- allorthohmb3$gene_m
rownames(exp_mtx.human1) <- allorthohmh1$gene_m
rownames(exp_mtx.human2) <- allorthohmh2$gene_m
rownames(exp_mtx.human3) <- allorthohmh3$gene_m

## Create the seurat object with human genes
baboon1.renamed <- CreateSeuratObject(counts = exp_mtx.baboon1, meta.data = baboon1@meta.data, project = "baboon1" )
baboon2.renamed <- CreateSeuratObject(counts = exp_mtx.baboon2, meta.data = baboon2@meta.data, project = "baboon2" )
baboon3.renamed <- CreateSeuratObject(counts = exp_mtx.baboon3, meta.data = baboon3@meta.data, project = "baboon3" )
human1.renamed <- CreateSeuratObject(counts = exp_mtx.human1, meta.data = human1@meta.data, project = "human1" )
human2.renamed <- CreateSeuratObject(counts = exp_mtx.human2, meta.data = human2@meta.data, project = "human2" )
human3.renamed <- CreateSeuratObject(counts = exp_mtx.human3, meta.data = human3@meta.data, project = "human3" )

##### Merging objects #####
renamed.primate <- merge(baboon1.renamed, y = c(baboon2.renamed, baboon3.renamed, human1.renamed, human2.renamed, human3.renamed),  project = "PRIMATE_NEU_EC_mus")


##### Integration #####

### Rerunning normalization with SCTransform ###

renamed.primate <- SCTransform(renamed.primate, method = "glmGamPoi",verbose = T)

# Before you start integration, make sure you are in the correct assay, which is
# "SCT" as you used SCTransform to normalize the count data.

DefaultAssay(renamed.primate) <- "SCT"

### Dimensional reduction before integration ###
renamed.primate <- RunPCA(renamed.primatee)
renamed.primate <- FindNeighbors(renamed.primate, dims = 1:30, reduction = "pca")
renamed.primate <- FindClusters(renamed.primate, resolution = 2, cluster.name = "unintegrated_clusters")
renamed.primate <- RunUMAP(renamed.primate, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")


### RPCA ###
renamed.primate <- IntegrateLayers(object = renamed.primate, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca",
                           normalization.method = "SCT", verbose = TRUE)

### Re-join layers - nex step
renamed.primate[["RNA"]] <- JoinLayers(renamed.primate[["RNA"]])


##### Find neighbors and clusters
renamed.primate <- FindNeighbors(renamed.primate, reduction = "integrated.rpca", dims = 1:30)

renamed.primate <- FindClusters(renamed.primate, resolution = 0.1, cluster.name = "rpca_cluster0.1")
renamed.primate <- FindClusters(renamed.primate, resolution = 0.2, cluster.name = "rpca_cluster0.2")
renamed.primate <- FindClusters(renamed.primate, resolution = 0.3, cluster.name = "rpca_cluster0.3")
renamed.primate <- FindClusters(renamed.primate, resolution = 0.8, cluster.name = "rpca_cluster0.8")

## run UMAP on RPCA integration
renamed.primate <- RunUMAP(renamed.primate, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")


##### Label renamed primate neuronal data set with Allen mouse brain atlas EC #####

### Annotation to Allen mouse brain atlas EC reference ###

# import reference.
mouse_atlas.neu <- readRDS("Data_EC_atlas/Seurat_obj/mouseAllen/mouse_atlasv5_neu.rds")

#In this example, the qyuery seurat object is named re
# Make sure your reference and query are in the SCT assay
DefaultAssay(renamed.primate) <- "SCT"
DefaultAssay(mouse_atlas.neu) <- "SCT"

# Find transfer anchors
re.anchors <- FindTransferAnchors(reference = mouse_atlas.neu, 
                                  query = renamed.primate,
                                  normalization.method = 'SCT',
                                  dims = 1:30)

# below is the workflow for generating prediction scores for each cell in your query
# object. I have done this 3 times as I want my query object you have 3 sets of cell type ids.

# metadata #1 : class_label
predictions.class <- TransferData(anchorset = re.anchors, 
                                  refdata =  mouse_atlas.neu$class_label,
                                  dims = 1:30)

renamed.primate <- AddMetaData(renamed.primate, metadata = predictions.class)

renamed.primate$class_label <- renamed.primate$predicted.id

# metadata #2 : subclass_label
predictions.subclass <- TransferData(anchorset = re.anchors, 
                                     refdata = mouse_atlas.neu$subclass_label,
                                     dims = 1:30)

renamed.primate <- AddMetaData(renamed.primate, metadata = predictions.subclass)

primate$subclass_label <- primate$predicted.id

### Save Seurat object
saveRDS(renamed.primate, file = "Data_EC_atlas/Seurat_obj/primate/renamed.primate.neu_annot.rds")


##### Save stripped down object for further integration #####

merged.primate$evolution <- "PRIMATE"

# switch to RNA assay for stripping data set down
DefaultAssay(merged.primate) <- "RNA"

# Strip down your data set
renamed.primate.diet <- DietSeurat(merged.primate, assays = "RNA", dimreducs = NULL)
renamed.primate.diet[["RNA"]] <- split(renamed.primate.diet[["RNA"]], f = renamed.primate.diet$orig.ident)

# Subset species
Idents(renamed.primate.diet) <- renamed.primate.diet$species

baboon.diet <- subset(renamed.primate.diet, idents = "BABOON")
human.diet <- subset(renamed.primate.diet, idents = "HUMAN")


# Save your stripped data set
saveRDS(baboon.diet, file = "Data_EC_atlas/Seurat_obj/primate/rerenamed.baboon.neu_annot_diet.rds")
saveRDS(human.diet, file = "Data_EC_atlas/Seurat_obj/primate/renamed.human.neu_annot_diet.rds")


##### Save final annotated object #####
merged.primate@meta.data <- merged.primate@meta.data %>%
  select(orig.ident, percent.mt, nCount_SCT, nFeature_SCT, sample, species, major_cell_types, neuronal_celltype_species, evolution)

saveRDS(merged.primate, file = "Data_EC_atlas/Seurat_obj/primate/renamed_primate_annot.rds")

# Subset EXN and IN
exn.primate <- subset(merged.primate, idents = c("Pvalb", "Vip", "Sst", "Lamp5", "Lamp5_Lhx6", "Pax6_Adarb2", "Pvalb_Vipr2", "Sst_Chodl"), invert = T)


in.primate <- subset(merged.primate, idents = c("Pvalb", "Vip", "Sst", "Lamp5", "Lamp5_Lhx6", "Pax6_Adarb2", "Pvalb_Vipr2", "Sst_Chodl"))


saveRDS(exn.primate, file = "Data_EC_atlas/Seurat_obj/primate/renamed.primateEXN_annot.rds")

saveRDS(in.primate, file = "Data_EC_atlas/Seurat_obj/primate/renamed.primateIN_annot.rds")



