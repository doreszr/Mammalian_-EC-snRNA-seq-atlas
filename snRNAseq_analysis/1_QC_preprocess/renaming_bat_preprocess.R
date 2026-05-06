########################################################################################################################

### Author: Dorottya Ralbovszki
### Date:01-25

### QC and doublet detection by Doubletfinder

### Rename fruit bat genes to mouse orthologues in all fruit bat samples

### QC and doublet detection by Doubletfinder

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
library(DoubletFinder)
library(reticulate)


#### QC of fruit bat samples before renaming and saving objects for renaming ####

### Import your data in the form of GEX count matricies
bat3.data <- Read10X("Data_EC_atlas/count_matrix/fruitbat/EFB3/filtered_feature_bc_matrix")
bat4.data <- Read10X("Data_EC_atlas/count_matrix/fruitbat/EFB4/filtered_feature_bc_matrix")
bat5.data <- Read10X("Data_EC_atlas/count_matrix/fruitbat/EFB5/filtered_feature_bc_matrix")

### Create your Seurat objects
bat3 <- CreateSeuratObject(counts = bat3.data,
                           project = "bat3",
                           min.cells = 3,
                           min.features = 200)
bat4 <- CreateSeuratObject(counts = bat4.data,
                           project = "bat4",
                           min.cells = 3,
                           min.features = 200)
bat5 <- CreateSeuratObject(counts = bat5.data,
                           project = "bat5",
                           min.cells = 3,
                           min.features = 200)

# save for renaming
saveRDS(bat3, file = "Data_EC_atlas/Seurat_obj/fruitbat/bat3raw.rds")
saveRDS(bat4, file = "Data_EC_atlas/Seurat_obj/fruitbat/bat4raw.rds")
saveRDS(bat5, file = "Data_EC_atlas/Seurat_obj/fruitbat/bat5raw.rds")

### Find % mitochondrial reads in each cell
# Create the set of mitochondrial genes
mito.list5 <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5")
mito.list34 <- c("ATP6", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5")

bat3[["percent.mt"]] <- PercentageFeatureSet(bat3, features = mito.list34)
bat4[["percent.mt"]] <- PercentageFeatureSet(bat4, features = mito.list34)
bat5[["percent.mt"]] <- PercentageFeatureSet(bat5, features = mito.list5)

### Visualize QC metrics
VlnPlot(bat3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(bat4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(bat5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

### Subset Seurat objects

# Insert the values you decided on in the last step into the function for each sample
bat3 <- subset(
  x = bat3,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 800 &
    nFeature_RNA < 4000 &
    percent.mt < 2
)
bat4 <- subset(
  x = bat4,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 800 &
    nFeature_RNA < 3000 &
    percent.mt < 2
)
bat5 <- subset(
  x = bat5,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 3000 &
    percent.mt < 2
)

### Normalize with SCTransform
bat3 <- SCTransform(bat3, method = "glmGamPoi",verbose = T)
bat4 <- SCTransform(bat4, method = "glmGamPoi",verbose = T)
bat5 <- SCTransform(bat5, method = "glmGamPoi",verbose = T)

### Do preemptive dimensional reduction for Doubletfinder
bat3 <- RunPCA(bat3,verbose = T)
bat4 <- RunPCA(bat4,verbose = T)
bat5 <- RunPCA(bat5,verbose = T)

bat3 <- RunUMAP(bat3, dims = 1:30, verbose = T)
bat4 <- RunUMAP(bat4, dims = 1:30, verbose = T)
bat5 <- RunUMAP(bat5, dims = 1:30, verbose = T)

#check your UMAP
DimPlot(bat3)
DimPlot(bat4)
DimPlot(bat5)

### Doubletfinder workflow

# Note: This workflow can be put into an pair of lapply functions, but that makes
# it slightly less plug-and-play with find and replace

# for more info: https://github.com/chris-mcginnis-ucsf/DoubletFinder

## pK Identification (no ground-truth)
sweep.res.bat3 <- paramSweep(bat3, PCs = 1:30, sct = T)
sweep.res.bat4 <- paramSweep(bat4, PCs = 1:30, sct = T)
sweep.res.bat5 <- paramSweep(bat5, PCs = 1:30, sct = T)

sweep.stats.bat3 <- summarizeSweep(sweep.res.bat3, GT = FALSE)
sweep.stats.bat4 <- summarizeSweep(sweep.res.bat4, GT = FALSE)
sweep.stats.bat5 <- summarizeSweep(sweep.res.bat5, GT = FALSE)

bcmvn_bat3 <- find.pK(sweep.stats.bat3)
bcmvn_bat4 <- find.pK(sweep.stats.bat4)
bcmvn_bat5 <- find.pK(sweep.stats.bat5)

## save pK values to be used
bat3pk <- as.numeric(as.vector(top_n(bcmvn_bat3,1,BCmetric)[["pK"]]))
bat4pk <- as.numeric(as.vector(top_n(bcmvn_bat4,1,BCmetric)[["pK"]]))
bat5pk <- as.numeric(as.vector(top_n(bcmvn_bat5,1,BCmetric)[["pK"]]))


# Note: Below, change the 0.075 to a value that represents the expected doublet
# rate in the appropriate 10X Genomics protocol. 7.5% is for 3' targeting 10k cells.
nExp_poi.bat3 <- round(0.075*nrow(bat3@meta.data))
nExp_poi.bat4 <- round(0.075*nrow(bat4@meta.data))
nExp_poi.bat5 <- round(0.075*nrow(bat5@meta.data))

## run doubletFinder
bat3 <- doubletFinder(bat3, PCs = 1:30, pN = 0.25, pK = bat3pk, nExp = nExp_poi.bat3, reuse.pANN = FALSE, sct = T)
bat4 <- doubletFinder(bat4, PCs = 1:30, pN = 0.25, pK = bat4pk, nExp = nExp_poi.bat4, reuse.pANN = FALSE, sct = T)
bat5 <- doubletFinder(bat5, PCs = 1:30, pN = 0.25, pK = bat5pk, nExp = nExp_poi.bat5, reuse.pANN = FALSE, sct = T)

## change active identity for subsetting
Idents(bat3) <- bat3@meta.data[[paste0("DF.classifications_0.25_",bat3pk,"_",nExp_poi.bat3)]]
Idents(bat4) <- bat4@meta.data[[paste0("DF.classifications_0.25_",bat4pk,"_",nExp_poi.bat4)]]
Idents(bat5) <- bat5@meta.data[[paste0("DF.classifications_0.25_",bat5pk,"_",nExp_poi.bat5)]]

## Visualize doublet detection
DimPlot(bat3)
DimPlot(bat4)
DimPlot(bat5)

#### RENAMING from bat to mouse gene names ####

##### Import your data in the form of GEX count matricies #####
bat3.data <- Read10X("Data_EC_atlas/count_matrix/fruitbat/EFB3/filtered_feature_bc_matrix")
bat4.data <- Read10X("Data_EC_atlas/count_matrix/fruitbat/EFB4/filtered_feature_bc_matrix")
bat5.data <- Read10X("Data_EC_atlas/count_matrix/fruitbat/EFB5/filtered_feature_bc_matrix")

##### Create your Seurat objects #####
bat3 <- CreateSeuratObject(counts = bat3.data,
                           project = "bat3",
                           min.cells = 3,
                           min.features = 200)
bat4 <- CreateSeuratObject(counts = bat4.data,
                           project = "bat4",
                           min.cells = 3,
                           min.features = 200)
bat5 <- CreateSeuratObject(counts = bat5.data,
                           project = "bat5",
                           min.cells = 3,
                           min.features = 200)


# load ortho df bat vs mouse
allorthobatm <- read.csv("Data_EC_atlas/orthologues/mappingortho/bat_vs_mouse/allorthobatm.csv", sep = ";")


# create a count matrix from RNA assaq
exp_mtx.bat3 <- as.matrix(bat3[["RNA"]]$counts)
exp_mtx.bat4 <- as.matrix(bat4[["RNA"]]$counts)
exp_mtx.bat5 <- as.matrix(bat5[["RNA"]]$counts)

# get the genes that are present in the count matrix
allorthobatm3 <- filter(allorthobatm,
                        gene_bat %in% rownames(exp_mtx.bat3))
allorthobatm4 <- filter(allorthobatm,
                        gene_bat %in% rownames(exp_mtx.bat4))
allorthobatm5 <- filter(allorthobatm,
                        gene_bat %in% rownames(exp_mtx.bat5))


## Filter the expression matrix for genes which a mouse counterpart is available
exp_mtx.bat3 <- exp_mtx.bat3[allorthobatm3$gene_bat,]
exp_mtx.bat4 <- exp_mtx.bat4[allorthobatm4$gene_bat,]
exp_mtx.bat5 <- exp_mtx.bat5[allorthobatm5$gene_bat,]



## Now chnage the rownames of the matrix to the mouse gene names
rownames(exp_mtx.bat3) <- allorthobatm3$gene_m
rownames(exp_mtx.bat4) <- allorthobatm4$gene_m
rownames(exp_mtx.bat5) <- allorthobatm5$gene_m



## Create the seurat object with mouse genes
bat3.renamed <- CreateSeuratObject(counts = exp_mtx.bat3, meta.data = bat3@meta.data, project = "bat3" )
bat4.renamed <- CreateSeuratObject(counts = exp_mtx.bat4, meta.data = bat4@meta.data, project = "bat4" )
bat5.renamed <- CreateSeuratObject(counts = exp_mtx.bat5, meta.data = bat5@meta.data, project = "bat5" )

# rename the mousezed bat seurat objects so the bat workflow can be run on them
bat3 <- bat3.renamed
bat4 <- bat4.renamed
bat5 <- bat5.renamed

##### Find % mitochondrial reads in each cell #####
bat3[["percent.mt"]] <- PercentageFeatureSet(bat3, pattern = "^mt-")
bat4[["percent.mt"]] <- PercentageFeatureSet(bat4, pattern = "^mt-")
bat5[["percent.mt"]] <- PercentageFeatureSet(bat5, pattern = "^mt-")

##### Visualize QC metrics #####
VlnPlot(bat3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(bat4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(bat5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##### Subset Seurat objects #####

# Insert the values you decided on in the last step into the function for each sample
bat3 <- subset(
  x = bat3,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 800 &
    nFeature_RNA < 4000 &
    percent.mt < 2
)
bat4 <- subset(
  x = bat4,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 800 &
    nFeature_RNA < 3000 &
    percent.mt < 2
)
bat5 <- subset(
  x = bat5,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 3000 &
    percent.mt < 2
)

##### Normalize with SCTransform #####
bat3 <- SCTransform(bat3, method = "glmGamPoi",verbose = T)
bat4 <- SCTransform(bat4, method = "glmGamPoi",verbose = T)
bat5 <- SCTransform(bat5, method = "glmGamPoi",verbose = T)

##### Do preemptive dimensional reduction for Doubletfinder #####
bat3 <- RunPCA(bat3,verbose = T)
bat4 <- RunPCA(bat4,verbose = T)
bat5 <- RunPCA(bat5,verbose = T)

bat3 <- RunUMAP(bat3, dims = 1:30, verbose = T)
bat4 <- RunUMAP(bat4, dims = 1:30, verbose = T)
bat5 <- RunUMAP(bat5, dims = 1:30, verbose = T)

#check your UMAP
DimPlot(bat3)
DimPlot(bat4)
DimPlot(bat5)

##### Doubletfinder workflow #####

# Note: This workflow can be put into an pair of lapply functions, but that makes
# it slightly less plug-and-play with find and replace

# for more info: https://github.com/chris-mcginnis-ucsf/DoubletFinder

## pK Identification (no ground-truth)
sweep.res.bat3 <- paramSweep(bat3, PCs = 1:30, sct = T)
sweep.res.bat4 <- paramSweep(bat4, PCs = 1:30, sct = T)
sweep.res.bat5 <- paramSweep(bat5, PCs = 1:30, sct = T)

sweep.stats.bat3 <- summarizeSweep(sweep.res.bat3, GT = FALSE)
sweep.stats.bat4 <- summarizeSweep(sweep.res.bat4, GT = FALSE)
sweep.stats.bat5 <- summarizeSweep(sweep.res.bat5, GT = FALSE)

bcmvn_bat3 <- find.pK(sweep.stats.bat3)
bcmvn_bat4 <- find.pK(sweep.stats.bat4)
bcmvn_bat5 <- find.pK(sweep.stats.bat5)

## save pK values to be used
bat3pk <- as.numeric(as.vector(top_n(bcmvn_bat3,1,BCmetric)[["pK"]]))
bat4pk <- as.numeric(as.vector(top_n(bcmvn_bat4,1,BCmetric)[["pK"]]))
bat5pk <- as.numeric(as.vector(top_n(bcmvn_bat5,1,BCmetric)[["pK"]]))


# Note: Below, change the 0.075 to a value that represents the expected doublet
# rate in the appropriate 10X Genomics protocol. 7.5% is for 3' targeting 10k cells.
nExp_poi.bat3 <- round(0.075*nrow(bat3@meta.data))
nExp_poi.bat4 <- round(0.075*nrow(bat4@meta.data))
nExp_poi.bat5 <- round(0.075*nrow(bat5@meta.data))

## run doubletFinder
bat3 <- doubletFinder(bat3, PCs = 1:30, pN = 0.25, pK = bat3pk, nExp = nExp_poi.bat3, reuse.pANN = FALSE, sct = T)
bat4 <- doubletFinder(bat4, PCs = 1:30, pN = 0.25, pK = bat4pk, nExp = nExp_poi.bat4, reuse.pANN = FALSE, sct = T)
bat5 <- doubletFinder(bat5, PCs = 1:30, pN = 0.25, pK = bat5pk, nExp = nExp_poi.bat5, reuse.pANN = FALSE, sct = T)

## change active identity for subsetting
Idents(bat3) <- bat3@meta.data[[paste0("DF.classifications_0.25_",bat3pk,"_",nExp_poi.bat3)]]
Idents(bat4) <- bat4@meta.data[[paste0("DF.classifications_0.25_",bat4pk,"_",nExp_poi.bat4)]]
Idents(bat5) <- bat5@meta.data[[paste0("DF.classifications_0.25_",bat5pk,"_",nExp_poi.bat5)]]

## Visualize doublet detection
DimPlot(bat3)
DimPlot(bat4)
DimPlot(bat5)

FeaturePlot(bat3, features = "nFeature_SCT")
FeaturePlot(bat4, features = "nFeature_SCT")
FeaturePlot(bat5, features = "nFeature_SCT")

# remove duplets
bat3 <- subset(bat3, idents = "Singlet")
bat4 <- subset(bat4, idents = "Singlet")
bat5 <- subset(bat5, idents = "Singlet")


## Visualize singlets
DimPlot(bat3)
DimPlot(bat4)
DimPlot(bat5)


saveRDS(bat3, file = "Data_EC_atlas/Seurat_obj/fruitbat/forplotrenamedbat3.rds")
saveRDS(bat4, file = "Data_EC_atlas/Seurat_obj/fruitbat/forplotrenamedbat4.rds")
saveRDS(bat5, file = "Data_EC_atlas/Seurat_obj/fruitbat/forplotrenamedbat5.rds")


##### Setting the active assay back to RNA #####
DefaultAssay(bat3) <- "RNA"
DefaultAssay(bat4) <- "RNA"
DefaultAssay(bat5) <- "RNA"

##### Deleting SCT assay #####
bat3@assays[["SCT"]] <- NULL
bat4@assays[["SCT"]] <- NULL
bat5@assays[["SCT"]] <- NULL

bat3@reductions[["pca"]] <- NULL
bat3@reductions[["umap"]] <- NULL
bat4@reductions[["pca"]] <- NULL
bat5@reductions[["pca"]] <- NULL
bat4@reductions[["umap"]] <- NULL
bat5@reductions[["umap"]] <- NULL


##### save progress before merging #####

saveRDS(bat3, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamedbat3.rds")
saveRDS(bat4, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamedbat4.rds")
saveRDS(bat5, file = "Data_EC_atlas/Seurat_obj/fruitbat/renamedbat5.rds")
