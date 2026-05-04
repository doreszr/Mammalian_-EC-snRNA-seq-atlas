
########################################################################################################################

### Author: Dorottya Ralbovszki
### Date: 11-24

### Preprocessing filtered count matrix from published human samples

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


##### Import your data in the form of GEX count matricies #####
human1.data <- Read10X("Data_EC_atlas/count_matrix/human/SRR16928876/filtered_feature_bc_matrix")
human2.data <- Read10X("Data_EC_atlas/count_matrix/human/SRR16928884/filtered_feature_bc_matrix")
human3.data <- Read10X("Data_EC_atlas/count_matrix/human/SRR16928886/filtered_feature_bc_matrix")

##### Create your Seurat objects #####
human1 <- CreateSeuratObject(counts = human1.data,
                             project = "human1",
                             min.cells = 3,
                             min.features = 200)
human2 <- CreateSeuratObject(counts = human2.data,
                             project = "human2",
                             min.cells = 3,
                             min.features = 200)
human3 <- CreateSeuratObject(counts = human3.data,
                             project = "human3",
                             min.cells = 3,
                             min.features = 200)

# save for ortholog mapping
saveRDS(human1, file = "Data_EC_atlas/Seurat_obj/human/human1raw.rds")
saveRDS(human2, file = "Data_EC_atlas/Seurat_obj/human/human2raw.rds")
saveRDS(human3, file = "Data_EC_atlas/Seurat_obj/human/human3raw.rds")


##### Find % mitochondrial reads in each cell #####
human1[["percent.mt"]] <- PercentageFeatureSet(human1, pattern = "^MT-")
human2[["percent.mt"]] <- PercentageFeatureSet(human2, pattern = "^MT-")
human3[["percent.mt"]] <- PercentageFeatureSet(human3, pattern = "^MT-")

##### Visualize QC metrics #####
VlnPlot(human1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(human2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(human3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# This is where you need to take a look at the violin plots to make sure you are
# including data you want and excluding data that you don't want. Pick a cut-off
# for counts/cell, a min. and max. for # of features/cell, and 

##### Subset Seurat objects #####

# Insert the values you decided on in the last step into the function for each sample
human1 <- subset(
  x = human1,
  subset = nCount_RNA < 20000 &
    nFeature_RNA > 1200 &
    nFeature_RNA < 6500 &
    percent.mt < 5
)
human2 <- subset(
  x = human2,
  subset = nCount_RNA < 25000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 8000 &
    percent.mt < 5
)
human3 <- subset(
  x = human3,
  subset = nCount_RNA < 35000 &
    nFeature_RNA > 1000 &
    nFeature_RNA < 8000 &
    percent.mt < 5
)

##### Normalize with SCTransform #####
human1 <- SCTransform(human1, method = "glmGamPoi",verbose = T)
human2 <- SCTransform(human2, method = "glmGamPoi",verbose = T)
human3 <- SCTransform(human3, method = "glmGamPoi",verbose = T)

##### Do preemptive dimensional reduction for Doubletfinder #####
human1 <- RunPCA(human1,verbose = T)
human2 <- RunPCA(human2,verbose = T)
human3 <- RunPCA(human3,verbose = T)

human1 <- RunUMAP(human1, dims = 1:30, verbose = T)
human2 <- RunUMAP(human2, dims = 1:30, verbose = T)
human3 <- RunUMAP(human3, dims = 1:30, verbose = T)

#check your UMAP
DimPlot(human1)
DimPlot(human2)
DimPlot(human3)

##### Doubletfinder workflow #####

# Note: This workflow can be put into an pair of lapply functions, but that makes
# it slightly less plug-and-play with find and replace

# for more info: https://github.com/chris-mcginnis-ucsf/DoubletFinder

## pK Identification (no ground-truth)
sweep.res.human1 <- paramSweep(human1, PCs = 1:30, sct = T)
sweep.res.human2 <- paramSweep(human2, PCs = 1:30, sct = T)
sweep.res.human3 <- paramSweep(human3, PCs = 1:30, sct = T)

sweep.stats.human1 <- summarizeSweep(sweep.res.human1, GT = FALSE)
sweep.stats.human2 <- summarizeSweep(sweep.res.human2, GT = FALSE)
sweep.stats.human3 <- summarizeSweep(sweep.res.human3, GT = FALSE)

bcmvn_human1 <- find.pK(sweep.stats.human1)
bcmvn_human2 <- find.pK(sweep.stats.human2)
bcmvn_human3 <- find.pK(sweep.stats.human3)

## save pK values to be used
human1pk <- as.numeric(as.vector(top_n(bcmvn_human1,1,BCmetric)[["pK"]]))
human2pk <- as.numeric(as.vector(top_n(bcmvn_human2,1,BCmetric)[["pK"]]))
human3pk <- as.numeric(as.vector(top_n(bcmvn_human3,1,BCmetric)[["pK"]]))


# Note: Below, change the 0.075 to a value that represents the expected doublet
# rate in the appropriate 10X Genomics protocol. 7.5% is for 3' targeting 10k cells.
nExp_poi.human1 <- round(0.075*nrow(human1@meta.data))
nExp_poi.human2 <- round(0.075*nrow(human2@meta.data))
nExp_poi.human3 <- round(0.075*nrow(human3@meta.data))

## run doubletFinder
human1 <- doubletFinder(human1, PCs = 1:30, pN = 0.25, pK = human1pk, nExp = nExp_poi.human1, reuse.pANN = FALSE, sct = T)
human2 <- doubletFinder(human2, PCs = 1:30, pN = 0.25, pK = human2pk, nExp = nExp_poi.human2, reuse.pANN = FALSE, sct = T)
human3 <- doubletFinder(human3, PCs = 1:30, pN = 0.25, pK = human3pk, nExp = nExp_poi.human3, reuse.pANN = FALSE, sct = T)

## change active identity for subsetting
Idents(human1) <- human1@meta.data[[paste0("DF.classifications_0.25_",human1pk,"_",nExp_poi.human1)]]
Idents(human2) <- human2@meta.data[[paste0("DF.classifications_0.25_",human2pk,"_",nExp_poi.human2)]]
Idents(human3) <- human3@meta.data[[paste0("DF.classifications_0.25_",human3pk,"_",nExp_poi.human3)]]

## Visualize doublet detection
DimPlot(human1)
DimPlot(human2)
DimPlot(human3)

FeaturePlot(human1, features = "nFeature_SCT")
FeaturePlot(human2, features = "nFeature_SCT")
FeaturePlot(human3, features = "nFeature_SCT")


human1 <- subset(human1, idents = "Singlet")
human2 <- subset(human2, idents = "Singlet")
human3 <- subset(human3, idents = "Singlet")


##### Setting the active assay back to RNA #####
DefaultAssay(human1) <- "RNA"
DefaultAssay(human2) <- "RNA"
DefaultAssay(human3) <- "RNA"

##### Deleting SCT assay #####
human1@assays[["SCT"]] <- NULL
human2@assays[["SCT"]] <- NULL
human3@assays[["SCT"]] <- NULL

##### Deleting dimensional reductions made for duplet detection #####
human1@reductions[["pca"]] <- NULL
human1@reductions[["umap"]] <- NULL
human2@reductions[["pca"]] <- NULL
human2@reductions[["umap"]] <- NULL
human3@reductions[["pca"]] <- NULL
human3@reductions[["umap"]] <- NULL

##### save progress before merging #####

# It's a great idea to save your individual objects in case you want to make an
# merged object with a specific set of samples later.

saveRDS(human1, file = "Data_EC_atlas/Seurat_obj/human/human1.rds")
saveRDS(human2, file = "Data_EC_atlas/Seurat_obj/human/human2.rds")
saveRDS(human3, file = "Data_EC_atlas/Seurat_obj/human/human3.rds")
