########################################################################################################################

### Author: Dorottya Ralbovszki
### Date: 01-25

### QC and doublet detection by Doubletfinder

### Rename baboon genes to human orthologues in all baboon samples

### QC and doublet detection by Doubletfinder of renamed baboon samples

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


#### QC of baboon samples before renaming and saving objects for renaming ####

### Import your data in the form of GEX count matricies 
baboon1.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B1/filtered_feature_bc_matrix")
baboon2.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B2/filtered_feature_bc_matrix")
baboon3.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B3/filtered_feature_bc_matrix")
baboon5.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B5/filtered_feature_bc_matrix")


### Create your Seurat objects
baboon1 <- CreateSeuratObject(counts = baboon1.data,
                              project = "baboon1",
                              min.cells = 3,
                              min.features = 200)
baboon2 <- CreateSeuratObject(counts = baboon2.data,
                              project = "baboon2",
                              min.cells = 3,
                              min.features = 200)
baboon3 <- CreateSeuratObject(counts = baboon3.data,
                              project = "baboon3",
                              min.cells = 3,
                              min.features = 200)
baboon5 <- CreateSeuratObject(counts = baboon5.data,
                              project = "baboon5",
                              min.cells = 3,
                              min.features = 200)


# save for renaming
saveRDS(baboon1, file = "Data_EC_atlas/Seurat_obj/baboon/rawbaboon1.rds")
saveRDS(baboon2, file = "Data_EC_atlas/Seurat_obj/baboon/rawbaboon2.rds")
saveRDS(baboon3, file = "Data_EC_atlas/Seurat_obj/baboon/rawbaboon3.rds")
saveRDS(baboon5, file = "Data_EC_atlas/Seurat_obj/baboon/rawbaboon5.rds")


##### Find % mitochondrial reads in each cell #####
# Create the set of mitochondrial genes
mito.list <- c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6")

baboon1[["percent.mt"]] <- PercentageFeatureSet(baboon1, features = mito.list)
baboon2[["percent.mt"]] <- PercentageFeatureSet(baboon2, features = mito.list)
baboon3[["percent.mt"]] <- PercentageFeatureSet(baboon3, features = mito.list)
baboon5[["percent.mt"]] <- PercentageFeatureSet(baboon5, features = mito.list)

### Visualize QC metrics 
VlnPlot(baboon1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(baboon2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(baboon3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(baboon5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##### Subset Seurat objects #####

# Insert the values you decided on in the last step into the function for each sample
baboon1 <- subset(
  x = baboon1,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 6000 &
    percent.mt < 5
)
baboon2 <- subset(
  x = baboon2,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 6000 &
    percent.mt < 5
)
baboon3 <- subset(
  x = baboon3,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 5000 &
    percent.mt < 5
)
baboon5 <- subset(
  x = baboon5,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 5000 &
    percent.mt < 5
)

### Normalize with SCTransform 
baboon1 <- SCTransform(baboon1, method = "glmGamPoi",verbose = T)
baboon2 <- SCTransform(baboon2, method = "glmGamPoi",verbose = T)
baboon3 <- SCTransform(baboon3, method = "glmGamPoi",verbose = T)
baboon5 <- SCTransform(baboon5, method = "glmGamPoi",verbose = T)

### Do preemptive dimensional reduction for Doubletfinder 
baboon1 <- RunPCA(baboon1,verbose = T)
baboon2 <- RunPCA(baboon2,verbose = T)
baboon3 <- RunPCA(baboon3,verbose = T)
baboon5 <- RunPCA(baboon5,verbose = T)

baboon1 <- RunUMAP(baboon1, dims = 1:30, verbose = T)
baboon2 <- RunUMAP(baboon2, dims = 1:30, verbose = T)
baboon3 <- RunUMAP(baboon3, dims = 1:30, verbose = T)
baboon5 <- RunUMAP(baboon5, dims = 1:30, verbose = T)

#check your UMAP
DimPlot(baboon1)
DimPlot(baboon2)
DimPlot(baboon3)
DimPlot(baboon5)

### Doubletfinder workflow 

# Note: This workflow can be put into an pair of lapply functions, but that makes
# it slightly less plug-and-play with find and replace

# for more info: https://github.com/chris-mcginnis-ucsf/DoubletFinder

## pK Identification (no ground-truth)
sweep.res.baboon1 <- paramSweep(baboon1, PCs = 1:30, sct = T)
sweep.res.baboon2 <- paramSweep(baboon2, PCs = 1:30, sct = T)
sweep.res.baboon3 <- paramSweep(baboon3, PCs = 1:30, sct = T)
sweep.res.baboon5 <- paramSweep(baboon5, PCs = 1:30, sct = T)

sweep.stats.baboon1 <- summarizeSweep(sweep.res.baboon1, GT = FALSE)
sweep.stats.baboon2 <- summarizeSweep(sweep.res.baboon2, GT = FALSE)
sweep.stats.baboon3 <- summarizeSweep(sweep.res.baboon3, GT = FALSE)
sweep.stats.baboon5 <- summarizeSweep(sweep.res.baboon5, GT = FALSE)

bcmvn_baboon1 <- find.pK(sweep.stats.baboon1)
bcmvn_baboon2 <- find.pK(sweep.stats.baboon2)
bcmvn_baboon3 <- find.pK(sweep.stats.baboon3)
bcmvn_baboon5 <- find.pK(sweep.stats.baboon5)

## save pK values to be used
baboon1pk <- as.numeric(as.vector(top_n(bcmvn_baboon1,1,BCmetric)[["pK"]]))
baboon2pk <- as.numeric(as.vector(top_n(bcmvn_baboon2,1,BCmetric)[["pK"]]))
baboon3pk <- as.numeric(as.vector(top_n(bcmvn_baboon3,1,BCmetric)[["pK"]]))
baboon5pk <- as.numeric(as.vector(top_n(bcmvn_baboon5,1,BCmetric)[["pK"]]))

# Note: Below, change the 0.075 to a value that represents the expected doublet
# rate in the appropriate 10X Genomics protocol. 7.5% is for 3' targeting 10k cells.
nExp_poi.baboon1 <- round(0.075*nrow(baboon1@meta.data))
nExp_poi.baboon2 <- round(0.075*nrow(baboon2@meta.data))
nExp_poi.baboon3 <- round(0.075*nrow(baboon3@meta.data))
nExp_poi.baboon5 <- round(0.075*nrow(baboon5@meta.data))

## run doubletFinder
baboon1 <- doubletFinder(baboon1, PCs = 1:30, pN = 0.25, pK = baboon1pk, nExp = nExp_poi.baboon1, reuse.pANN = FALSE, sct = T)
baboon2 <- doubletFinder(baboon2, PCs = 1:30, pN = 0.25, pK = baboon2pk, nExp = nExp_poi.baboon2, reuse.pANN = FALSE, sct = T)
baboon3 <- doubletFinder(baboon3, PCs = 1:30, pN = 0.25, pK = baboon3pk, nExp = nExp_poi.baboon3, reuse.pANN = FALSE, sct = T)
baboon5 <- doubletFinder(baboon5, PCs = 1:30, pN = 0.25, pK = baboon5pk, nExp = nExp_poi.baboon5, reuse.pANN = FALSE, sct = T)

## change active identity for subsetting
Idents(baboon1) <- baboon1@meta.data[[paste0("DF.classifications_0.25_",baboon1pk,"_",nExp_poi.baboon1)]]
Idents(baboon2) <- baboon2@meta.data[[paste0("DF.classifications_0.25_",baboon2pk,"_",nExp_poi.baboon2)]]
Idents(baboon3) <- baboon3@meta.data[[paste0("DF.classifications_0.25_",baboon3pk,"_",nExp_poi.baboon3)]]
Idents(baboon5) <- baboon5@meta.data[[paste0("DF.classifications_0.25_",baboon5pk,"_",nExp_poi.baboon5)]]

## Visualize doublet detection
DimPlot(baboon1)
DimPlot(baboon2)
DimPlot(baboon3)
DimPlot(baboon5)


##### Import your data in the form of GEX count matricies #####
baboon1.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B1/filtered_feature_bc_matrix")
baboon2.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B2/filtered_feature_bc_matrix")
baboon3.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B3/filtered_feature_bc_matrix")
baboon5.data <- Read10X("Data_EC_atlas/count_matrix/baboon/B5/filtered_feature_bc_matrix")


##### Create your Seurat objects #####
baboon1 <- CreateSeuratObject(counts = baboon1.data,
                              project = "baboon1",
                              min.cells = 3,
                              min.features = 200)
baboon2 <- CreateSeuratObject(counts = baboon2.data,
                              project = "baboon2",
                              min.cells = 3,
                              min.features = 200)
baboon3 <- CreateSeuratObject(counts = baboon3.data,
                              project = "baboon3",
                              min.cells = 3,
                              min.features = 200)
baboon5 <- CreateSeuratObject(counts = baboon5.data,
                              project = "baboon5",
                              min.cells = 3,
                              min.features = 200)


# load ortho df baboon vs human
allorthobh <- read.csv("Data_EC_atlas/orthologues/mappingortho/baboon_vs_humanallorthobh.csv", sep = ";")


# create a count matrix from RNA assaq
exp_mtx.baboon1 <- as.matrix(baboon1[["RNA"]]$counts)
exp_mtx.baboon2 <- as.matrix(baboon2[["RNA"]]$counts)
exp_mtx.baboon3 <- as.matrix(baboon3[["RNA"]]$counts)
exp_mtx.baboon5 <- as.matrix(baboon5[["RNA"]]$counts)

# get the genes that are present in the count matrix
allorthobh1 <- filter(allorthobh,
                      gene_b %in% rownames(exp_mtx.baboon1))
allorthobh2 <- filter(allorthobh,
                      gene_b %in% rownames(exp_mtx.baboon2))
allorthobh3 <- filter(allorthobh,
                      gene_b %in% rownames(exp_mtx.baboon3))
allorthobh5 <- filter(allorthobh,
                      gene_b %in% rownames(exp_mtx.baboon5))

## Filter the expression matrix for genes which a human counterpart is available
exp_mtx.baboon1 <- exp_mtx.baboon1[allorthobh1$gene_b,]
exp_mtx.baboon2 <- exp_mtx.baboon2[allorthobh2$gene_b,]
exp_mtx.baboon3 <- exp_mtx.baboon3[allorthobh3$gene_b,]
exp_mtx.baboon5 <- exp_mtx.baboon5[allorthobh5$gene_b,]


## Now chnage the rownames of the matrix to the human gene names
rownames(exp_mtx.baboon1) <- allorthobh1$gene_h
rownames(exp_mtx.baboon2) <- allorthobh2$gene_h
rownames(exp_mtx.baboon3) <- allorthobh3$gene_h
rownames(exp_mtx.baboon5) <- allorthobh5$gene_h


## Create the seurat object with human genes
baboon1.renamed <- CreateSeuratObject(counts = exp_mtx.baboon1, meta.data = baboon1@meta.data, project = "baboon1" )
baboon2.renamed <- CreateSeuratObject(counts = exp_mtx.baboon2, meta.data = baboon2@meta.data, project = "baboon2" )
baboon3.renamed <- CreateSeuratObject(counts = exp_mtx.baboon3, meta.data = baboon3@meta.data, project = "baboon3" )
baboon5.renamed <- CreateSeuratObject(counts = exp_mtx.baboon5, meta.data = baboon5@meta.data, project = "baboon5" )

# rename the humanized baboon seurat objects so the baboon workflow can be run on them
baboon1 <- baboon1.renamed
baboon2 <- baboon2.renamed
baboon3 <- baboon3.renamed
baboon5 <- baboon5.renamed

##### Find % mitochondrial reads in each cell #####
baboon1[["percent.mt"]] <- PercentageFeatureSet(baboon1, pattern = "^MT-")
baboon2[["percent.mt"]] <- PercentageFeatureSet(baboon2, pattern = "^MT-")
baboon3[["percent.mt"]] <- PercentageFeatureSet(baboon3, pattern = "^MT-")
baboon5[["percent.mt"]] <- PercentageFeatureSet(baboon5, pattern = "^MT-")


##### Visualize QC metrics #####
VlnPlot(baboon1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(baboon2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(baboon3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(baboon5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##### Subset Seurat objects #####

# Insert the values you decided on in the last step into the function for each sample
baboon1 <- subset(
  x = baboon1,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 6000 &
    percent.mt < 5
)
baboon2 <- subset(
  x = baboon2,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 6000 &
    percent.mt < 5
)
baboon3 <- subset(
  x = baboon3,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 5000 &
    percent.mt < 5
)
baboon5 <- subset(
  x = baboon5,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 5000 &
    percent.mt < 5
)

##### Normalize with SCTransform #####
baboon1 <- SCTransform(baboon1, method = "glmGamPoi",verbose = T)
baboon2 <- SCTransform(baboon2, method = "glmGamPoi",verbose = T)
baboon3 <- SCTransform(baboon3, method = "glmGamPoi",verbose = T)
baboon5 <- SCTransform(baboon5, method = "glmGamPoi",verbose = T)

##### Do preemptive dimensional reduction for Doubletfinder #####
baboon1 <- RunPCA(baboon1,verbose = T)
baboon2 <- RunPCA(baboon2,verbose = T)
baboon3 <- RunPCA(baboon3,verbose = T)
baboon5 <- RunPCA(baboon5,verbose = T)

baboon1 <- RunUMAP(baboon1, dims = 1:30, verbose = T)
baboon2 <- RunUMAP(baboon2, dims = 1:30, verbose = T)
baboon3 <- RunUMAP(baboon3, dims = 1:30, verbose = T)
baboon5 <- RunUMAP(baboon5, dims = 1:30, verbose = T)

#check your UMAP
DimPlot(baboon1)
DimPlot(baboon2)
DimPlot(baboon3)
DimPlot(baboon5)

##### Doubletfinder workflow #####

# Note: This workflow can be put into an pair of lapply functions, but that makes
# it slightly less plug-and-play with find and replace

# for more info: https://github.com/chris-mcginnis-ucsf/DoubletFinder

## pK Identification (no ground-truth)
sweep.res.baboon1 <- paramSweep(baboon1, PCs = 1:30, sct = T)
sweep.res.baboon2 <- paramSweep(baboon2, PCs = 1:30, sct = T)
sweep.res.baboon3 <- paramSweep(baboon3, PCs = 1:30, sct = T)
sweep.res.baboon5 <- paramSweep(baboon5, PCs = 1:30, sct = T)

sweep.stats.baboon1 <- summarizeSweep(sweep.res.baboon1, GT = FALSE)
sweep.stats.baboon2 <- summarizeSweep(sweep.res.baboon2, GT = FALSE)
sweep.stats.baboon3 <- summarizeSweep(sweep.res.baboon3, GT = FALSE)
sweep.stats.baboon5 <- summarizeSweep(sweep.res.baboon5, GT = FALSE)

bcmvn_baboon1 <- find.pK(sweep.stats.baboon1)
bcmvn_baboon2 <- find.pK(sweep.stats.baboon2)
bcmvn_baboon3 <- find.pK(sweep.stats.baboon3)
bcmvn_baboon5 <- find.pK(sweep.stats.baboon5)

## save pK values to be used
baboon1pk <- as.numeric(as.vector(top_n(bcmvn_baboon1,1,BCmetric)[["pK"]]))
baboon2pk <- as.numeric(as.vector(top_n(bcmvn_baboon2,1,BCmetric)[["pK"]]))
baboon3pk <- as.numeric(as.vector(top_n(bcmvn_baboon3,1,BCmetric)[["pK"]]))
baboon5pk <- as.numeric(as.vector(top_n(bcmvn_baboon5,1,BCmetric)[["pK"]]))

# Note: Below, change the 0.075 to a value that represents the expected doublet
# rate in the appropriate 10X Genomics protocol. 7.5% is for 3' targeting 10k cells.
nExp_poi.baboon1 <- round(0.075*nrow(baboon1@meta.data))
nExp_poi.baboon2 <- round(0.075*nrow(baboon2@meta.data))
nExp_poi.baboon3 <- round(0.075*nrow(baboon3@meta.data))
nExp_poi.baboon5 <- round(0.075*nrow(baboon5@meta.data))

## run doubletFinder
baboon1 <- doubletFinder(baboon1, PCs = 1:30, pN = 0.25, pK = baboon1pk, nExp = nExp_poi.baboon1, reuse.pANN = FALSE, sct = T)
baboon2 <- doubletFinder(baboon2, PCs = 1:30, pN = 0.25, pK = baboon2pk, nExp = nExp_poi.baboon2, reuse.pANN = FALSE, sct = T)
baboon3 <- doubletFinder(baboon3, PCs = 1:30, pN = 0.25, pK = baboon3pk, nExp = nExp_poi.baboon3, reuse.pANN = FALSE, sct = T)
baboon5 <- doubletFinder(baboon5, PCs = 1:30, pN = 0.25, pK = baboon5pk, nExp = nExp_poi.baboon5, reuse.pANN = FALSE, sct = T)

## change active identity for subsetting
Idents(baboon1) <- baboon1@meta.data[[paste0("DF.classifications_0.25_",baboon1pk,"_",nExp_poi.baboon1)]]
Idents(baboon2) <- baboon2@meta.data[[paste0("DF.classifications_0.25_",baboon2pk,"_",nExp_poi.baboon2)]]
Idents(baboon3) <- baboon3@meta.data[[paste0("DF.classifications_0.25_",baboon3pk,"_",nExp_poi.baboon3)]]
Idents(baboon5) <- baboon5@meta.data[[paste0("DF.classifications_0.25_",baboon5pk,"_",nExp_poi.baboon5)]]

## Visualize doublet detection
DimPlot(baboon1)
DimPlot(baboon2)
DimPlot(baboon3)
DimPlot(baboon5)

FeaturePlot(baboon1, features = "nFeature_SCT")
FeaturePlot(baboon2, features = "nFeature_SCT")
FeaturePlot(baboon3, features = "nFeature_SCT")
FeaturePlot(baboon5, features = "nFeature_SCT")

baboon1 <- subset(baboon1, idents = "Singlet")
baboon2 <- subset(baboon2, idents = "Singlet")
baboon3 <- subset(baboon3, idents = "Singlet")
baboon5 <- subset(baboon5, idents = "Singlet")


## Visualize singlets
DimPlot(baboon1)
DimPlot(baboon2)
DimPlot(baboon3)
DimPlot(baboon5)

## Save these objects for visualization
saveRDS(baboon1, file = "Data_EC_atlas/Seurat_obj/baboon/forplotrenamedbaboon1.rds")
saveRDS(baboon2, file = "Data_EC_atlas/Seurat_obj/baboon/forplotrenamedbaboon2.rds")
saveRDS(baboon3, file = "Data_EC_atlas/Seurat_obj/baboon/forplotrenamedbaboon3.rds")



##### Setting the active assay back to RNA #####
DefaultAssay(baboon1) <- "RNA"
DefaultAssay(baboon2) <- "RNA"
DefaultAssay(baboon3) <- "RNA"
DefaultAssay(baboon5) <- "RNA"

##### Deleting SCT assay #####
baboon1@assays[["SCT"]] <- NULL
baboon2@assays[["SCT"]] <- NULL
baboon3@assays[["SCT"]] <- NULL
baboon5@assays[["SCT"]] <- NULL


##### Deleting dimensional reductions made for duplet detection #####
baboon1@reductions[["pca"]] <- NULL
baboon1@reductions[["umap"]] <- NULL
baboon2@reductions[["pca"]] <- NULL
baboon2@reductions[["umap"]] <- NULL
baboon3@reductions[["pca"]] <- NULL
baboon3@reductions[["umap"]] <- NULL
baboon5@reductions[["pca"]] <- NULL
baboon5@reductions[["umap"]] <- NULL


##### Save objects for further analysis #####

saveRDS(baboon1, file = "Data_EC_atlas/Seurat_obj/baboon/renamedbaboon1.rds")
saveRDS(baboon2, file = "Data_EC_atlas/Seurat_obj/baboon/renamedbaboon2.rds")
saveRDS(baboon3, file = "Data_EC_atlas/Seurat_obj/baboon/renamedbaboon3.rds")
saveRDS(baboon5, file = "Data_EC_atlas/Seurat_obj/baboon/renamedbaboon5.rds")
