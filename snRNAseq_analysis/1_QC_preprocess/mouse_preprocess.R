########################################################################################################################

### Author: Dorottya Ralbovszki
### Date: 11-24

### Preprocessing filtered count matrix from published mouse samples

### QC by counts, mito reads, features and doublet detection

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



##### Import your data in the form of GEX count matricies #####
mouse1.data <- Read10X("Data_EC_atlas/count_matrix/mouse/SRR16409869_78/filtered_feature_bc_matrix")
mouse2.data <- Read10X("Data_EC_atlas/count_matrix/mouse/SRR16409870_80/filtered_feature_bc_matrix")
mouse3.data <- Read10X("Data_EC_atlas/count_matrix/mouse/SRR16409871_81/filtered_feature_bc_matrix")

##### Create your Seurat objects #####
mouse1 <- CreateSeuratObject(counts = mouse1.data,
                             project = "mouse1",
                             min.cells = 3,
                             min.features = 200)
mouse2 <- CreateSeuratObject(counts = mouse2.data,
                             project = "mouse2",
                             min.cells = 3,
                             min.features = 200)
mouse3 <- CreateSeuratObject(counts = mouse3.data,
                             project = "mouse3",
                             min.cells = 3,
                             min.features = 200)

# save for renaming
saveRDS(mouse1, file = "Data_EC_atlas/Seurat_obj/mouse/mouseraw1.rds")
saveRDS(mouse2, file = "Data_EC_atlas/Seurat_obj/mouse/mouseraw2.rds")
saveRDS(mouse3, file = "Data_EC_atlas/Seurat_obj/mouse/mouseraw3.rds")

##### Find % mitochondrial reads in each cell #####
mouse1[["percent.mt"]] <- PercentageFeatureSet(mouse1, pattern = "^mt-")
mouse2[["percent.mt"]] <- PercentageFeatureSet(mouse2, pattern = "^mt-")
mouse3[["percent.mt"]] <- PercentageFeatureSet(mouse3, pattern = "^mt-")

##### Visualize QC metrics #####
VlnPlot(mouse1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mouse2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(mouse3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##### Subset Seurat objects #####

# Insert the values you decided on in the last step into the function for each sample
mouse1 <- subset(
  x = mouse1,
  subset = nCount_RNA < 40000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 8000 &
    percent.mt < 5
)
mouse2 <- subset(
  x = mouse2,
  subset = nCount_RNA < 40000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 8000 &
    percent.mt < 5
)
mouse3 <- subset(
  x = mouse3,
  subset = nCount_RNA < 40000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 8000 &
    percent.mt < 5
)

##### Normalize with SCTransform #####
mouse1 <- SCTransform(mouse1, method = "glmGamPoi",verbose = T)
mouse2 <- SCTransform(mouse2, method = "glmGamPoi",verbose = T)
mouse3 <- SCTransform(mouse3, method = "glmGamPoi",verbose = T)

##### Do preemptive dimensional reduction for Doubletfinder #####
mouse1 <- RunPCA(mouse1,verbose = T)
mouse2 <- RunPCA(mouse2,verbose = T)
mouse3 <- RunPCA(mouse3,verbose = T)

mouse1 <- RunUMAP(mouse1, dims = 1:30, verbose = T)
mouse2 <- RunUMAP(mouse2, dims = 1:30, verbose = T)
mouse3 <- RunUMAP(mouse3, dims = 1:30, verbose = T)

#check your UMAP
DimPlot(mouse1)
DimPlot(mouse2)
DimPlot(mouse3)

##### Doubletfinder workflow #####

# Note: This workflow can be put into an pair of lapply functions, but that makes
# it slightly less plug-and-play with find and replace

# for more info: https://github.com/chris-mcginnis-ucsf/DoubletFinder

## pK Identification (no ground-truth)
sweep.res.mouse1 <- paramSweep(mouse1, PCs = 1:30, sct = T)
sweep.res.mouse2 <- paramSweep(mouse2, PCs = 1:30, sct = T)
sweep.res.mouse3 <- paramSweep(mouse3, PCs = 1:30, sct = T)

sweep.stats.mouse1 <- summarizeSweep(sweep.res.mouse1, GT = FALSE)
sweep.stats.mouse2 <- summarizeSweep(sweep.res.mouse2, GT = FALSE)
sweep.stats.mouse3 <- summarizeSweep(sweep.res.mouse3, GT = FALSE)

bcmvn_mouse1 <- find.pK(sweep.stats.mouse1)
bcmvn_mouse2 <- find.pK(sweep.stats.mouse2)
bcmvn_mouse3 <- find.pK(sweep.stats.mouse3)

## save pK values to be used
mouse1pk <- as.numeric(as.vector(top_n(bcmvn_mouse1,1,BCmetric)[["pK"]]))
mouse2pk <- as.numeric(as.vector(top_n(bcmvn_mouse2,1,BCmetric)[["pK"]]))
mouse3pk <- as.numeric(as.vector(top_n(bcmvn_mouse3,1,BCmetric)[["pK"]]))


# Note: Below, change the 0.075 to a value that represents the expected doublet
# rate in the appropriate 10X Genomics protocol. 7.5% is for 3' targeting 10k cells.
nExp_poi.mouse1 <- round(0.075*nrow(mouse1@meta.data))
nExp_poi.mouse2 <- round(0.075*nrow(mouse2@meta.data))
nExp_poi.mouse3 <- round(0.075*nrow(mouse3@meta.data))

## run doubletFinder
mouse1 <- doubletFinder(mouse1, PCs = 1:30, pN = 0.25, pK = mouse1pk, nExp = nExp_poi.mouse1, reuse.pANN = FALSE, sct = T)
mouse2 <- doubletFinder(mouse2, PCs = 1:30, pN = 0.25, pK = mouse2pk, nExp = nExp_poi.mouse2, reuse.pANN = FALSE, sct = T)
mouse3 <- doubletFinder(mouse3, PCs = 1:30, pN = 0.25, pK = mouse3pk, nExp = nExp_poi.mouse3, reuse.pANN = FALSE, sct = T)


## change active identity for subsetting
Idents(mouse1) <- mouse1@meta.data[[paste0("DF.classifications_0.25_",mouse1pk,"_",nExp_poi.mouse1)]]
Idents(mouse2) <- mouse2@meta.data[[paste0("DF.classifications_0.25_",mouse2pk,"_",nExp_poi.mouse2)]]
Idents(mouse3) <- mouse3@meta.data[[paste0("DF.classifications_0.25_",mouse3pk,"_",nExp_poi.mouse3)]]

## Visualize doublet detection
DimPlot(mouse1)
DimPlot(mouse2)
DimPlot(mouse3)

FeaturePlot(mouse1, features = "nFeature_SCT")
FeaturePlot(mouse2, features = "nFeature_SCT")
FeaturePlot(mouse3, features = "nFeature_SCT")

# remove duplets
mouse1 <- subset(mouse1, idents = "Singlet")
mouse2 <- subset(mouse2, idents = "Singlet")
mouse3 <- subset(mouse3, idents = "Singlet")

##### Setting the active assay back to RNA #####
DefaultAssay(mouse1) <- "RNA"
DefaultAssay(mouse2) <- "RNA"
DefaultAssay(mouse3) <- "RNA"

##### Deleting SCT assay #####
mouse1@assays[["SCT"]] <- NULL
mouse2@assays[["SCT"]] <- NULL
mouse3@assays[["SCT"]] <- NULL

##### Deleting reductions made for duplet detection #####
mouse1@reductions[["pca"]] <- NULL
mouse1@reductions[["umap"]] <- NULL
mouse2@reductions[["pca"]] <- NULL
mouse2@reductions[["umap"]] <- NULL
mouse3@reductions[["pca"]] <- NULL
mouse3@reductions[["umap"]] <- NULL

##### save progress before merging #####

saveRDS(mouse1, file = "Data_EC_atlas/Seurat_obj/mouse/mouse1.rds")
saveRDS(mouse2, file = "Data_EC_atlas/Seurat_obj/mouse/mouse2.rds")
saveRDS(mouse3, file = "Data_EC_atlas/Seurat_obj/mouse/mouse3.rds")