########################################################################################################################

### Marker gene expression plots

########################################################################################################################


library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

#### Vip species-specific markers ####

# Load objects
merged.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed_primate_annot.rds")

merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")


# Vip subtype Mouse vs Primate plotting function
save_featureplots <- function(seurat_obj, features,
                              outdir = "YOURPATH",
                              min.cutoff = NULL,
                              max.cutoff = NULL,
                              colors = c("lightgrey", "#C23A22"),
                              width = 6,
                              height = 5,
                              dpi = 300) {
  
  
  # capture Seurat object name
  obj_name <- deparse(substitute(seurat_obj))
  
  # compute global expression range
  expr <- FetchData(seurat_obj, vars = features)
  
  if (is.null(min.cutoff)) min.cutoff <- min(expr)
  if (is.null(max.cutoff)) max.cutoff <- max(expr)
  
  for (gene in features) {
    
    p <- FeaturePlot(
      seurat_obj,
      features = gene,
      reduction = "umap.rpca",
      min.cutoff = min.cutoff,
      max.cutoff = max.cutoff
    ) +
      scale_color_gradientn(
        colours = colors,
        limits = c(min.cutoff, max.cutoff)
      )
    
    filename <- paste0(outdir, "/", obj_name, gene, ".pdf")
    
    ggsave(
      filename,
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
}

# Marker genes
genes <- c("Pld5", "Gpd1", "Vip", "Calb2")

# Produce plots
save_featureplots(merged.primate, genes)
save_featureplots(merged.mouse, genes)






#### Supplementary gene expression plots ####

#Load objects
merged.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_annot.rds")

merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

merged.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed_primate_annot.rds")

merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")

# S2A
FeaturePlot(merged.primate, reduction = "umap.rpca", features = "Cux2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.primate, reduction = "umap.rpca", features = "Foxp2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.primate, reduction = "umap.rpca", features = "Rorb", cols = c("lightgrey", "#C23A22"))

# S2B
FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Cux2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Foxp2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Rorb", cols = c("lightgrey", "#C23A22"))

# S2C
FeaturePlot(merged.mouse, reduction = "umap.rpca", features = "Cux2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.mouse, reduction = "umap.rpca", features = "Foxp2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.mouse, reduction = "umap.rpca", features = "Rorb", cols = c("lightgrey", "#C23A22"))

# S2D
FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Maml2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Zbtb20", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.bat, reduction = "umap.rpca", features = "Cnr1", cols = c("lightgrey", "#C23A22")) # CGE marker based on allen

# S2E
FeaturePlot(merged.mammal, reduction = "umap.rpca", features = "Grik1", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.mammal, reduction = "umap.rpca", features = "Kcnip1", cols = c("lightgrey", "#C23A22"))

# S2F
FeaturePlot(merged.mammal, reduction = "umap.rpca", features = "Maml2", cols = c("lightgrey", "#C23A22"))

FeaturePlot(merged.mammal, reduction = "umap.rpca", features = "Zbtb20", cols = c("lightgrey", "#C23A22"))

