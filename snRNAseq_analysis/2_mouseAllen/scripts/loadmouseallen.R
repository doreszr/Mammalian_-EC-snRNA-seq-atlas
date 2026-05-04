## SCRIPT FOR READING LARGE hdf5 DATA SET INTO R

## Load necessary libraries
library(rhdf5)      # For reading hdf5; https://www.bioconductor.org/packages/devel/bioc/vignettes/rhdf5/inst/doc/rhdf5.html
library(HDF5Array)  # Alternatively option for reading the data matrix with less memory: https://rdrr.io/github/Bioconductor/HDF5Array/
library(data.table) # For fast reading of csv files
library(Seurat)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(ROCR)
#library(parallel)
#library(fields)
library(patchwork)
library(tidyverse)
#library(clustree)
#library(cowplot)
library(Matrix)

## Read in sample and gene names for use later
genes   <- h5read("mouseAllen/expression_matrix.hdf5","/data/gene")
samples <- h5read("mouseAllen/expression_matrix.hdf5","/data/samples")


## Read in the metadata using fread (fast!)
metadata <- fread("mouseAllenatlas/metadata_mouse.csv")
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sample_name
# Note that the order of metadata and counts and the number of cells are DIFFERENT 
#   (there are 107 cells with no metadata). This will be critical later!


## Subsample your data
## -- If you want to analyze these data in R and especially using Seurat, you will
## --   have a better chance of success if you only work on parts of the data at a
## --   time.  We suggest, using information in the metadata (e.g., cell type or
## --   brain region columns) to subset.   
## -- For this example code, I'm selecting 1000 random cells checking overlap with samples

metaEC <- filter(metadata,
                 region_label == "ENT")

use_samples <- intersect(rownames(metaEC), samples)
read_samples <- sort(match(use_samples,samples))
## Read in the count matrix in one of three ways


## Strategy #2: Read in only a relevant subset of data using h5read. This is the 
##   method that probably works best in most situations.
system.time({
  counts1 <- h5read("expression_matrix.hdf5", "/data/counts", index = list(read_samples,NULL))
  counts1 <- t(counts1)
  subcounts1 <- as(counts1, "dgCMatrix")
})  


## Add gene and sample names to the data matrix
## -- Note: We'll works with the recommended strategy (#2) for the rest of the script
rownames(subcounts1) <- as.character(genes)
colnames(subcounts1) <- as.character(samples) [read_samples]


## Read the subsetted data and metadata into Seurat
## -- Note that this **WILL NOT WORK** for the full 1 Million+ cell data set!
## -- I would strongly encourage other methods for analysis when dealing with more than ~100,000 cells at a time
## -- Also recall that the order of data and meta do not match, so reorder here
options(Seurat.object.assay.version = "v5")           # Make sure to create v5 object
seu <- CreateSeuratObject(counts=subcounts1)          # Put in a Seurat object
met <- as.data.frame(metadata[colnames(subcounts1),]) # Format the metadata
seu <- AddMetaData(seu,met)                           # Add the metadata


##### Normalize with SCTransform #####
seu <- SCTransform(seu, method = "glmGamPoi",verbose = T)

##### Do dimensional reduction #####
seu <- RunPCA(seu, verbose = T)
seu <- FindNeighbors(object = seu, dims = 1:30)
seu <- FindClusters(object = seu)
seu <- RunUMAP(seu, dims = 1:30, verbose = T)

saveRDS(seu, file = "Data_EC_atlas/Seurat_obj/2_mouseAllen/mouse_atlasv5.rds")
