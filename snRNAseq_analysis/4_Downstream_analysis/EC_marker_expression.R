################################################################################

# EC marker expression in species-specific and cross-species data sets

################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)




#### Species EC markers dot plots ####
EC.markers <- c("Cux2", "Pcdh8", "Lef1", "Chn2", "Plch1", "Fign",  "Reln", "Calb1", "Grik1", "Rorb", "Foxp2")

# Load and integrated species subsets
exn.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_EXN_annot.rds")
exn.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed.primateEXN_annot.rds")
exn.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouseEXN_annot.rds")
exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")

Idents(exn.bat) <- exn.bat$neuronal_celltype_species
Idents(exn.primate) <- exn.primate$neuronal_celltype_species
Idents(exn.mouse) <- exn.mouse$neuronal_celltype_species
Idents(exn.mammal) <- exn.mammal$neuronal_celltype

# BAT
desired_order_bat <- c("L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_Maml2", "L2/3_IT_ENTL_Cux2", "L2/3_IT_ENTm", "L2/3_IT", "L2/3_IT_Cald1",  "L2/3_IT_PPP", "L4/5_IT",
                       "L4-6", "L6b/CT_ENT", "L6b", "L2-6", "L2/3_L6b")
exn.bat@active.ident <- factor(Idents(exn.bat), levels = desired_order_bat)

dotplot <- DotPlot(exn.bat, features = EC.markers) +  coord_flip() + theme_minimal() + RotatedAxis() 
p3 <- dotplot + scale_color_gradient2(
  low = "blue",
  mid = "grey",
  high = "red",
  midpoint = 0
) + theme(axis.text.y = element_text(face = "italic", size = 20), axis.text.x = element_text(size = 15))


# MOUSE
desired_order_mouse <- c("L2_IT_ENTm", "L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_PPP", "L3_IT_ENT", "L5_PPP",
                         "L5/6_IT_TPE-ENT", "L6b/CT_ENT", "L6_IT_ENTL", "L6b_CTX", "NP_PPP", "CR")
exn.mouse@active.ident <- factor(Idents(exn.mouse), levels = desired_order_mouse)

dotplot <- DotPlot(exn.mouse, features = EC.markers) +  coord_flip() + theme_minimal() + RotatedAxis() 
p2 <- dotplot + scale_color_gradient2(
  low = "blue",
  mid = "grey",
  high = "red",
  midpoint = 0
) + theme(axis.text.y = element_text(face = "italic", size = 20), axis.text.x = element_text(size = 15))


# PRIMATE
desired_order_primate <- c("L2/3_IT_ENTL", "L2/3_IT_ENT", "L2/3_IT", "L2/3_IT_PPP", "L5/6_IT_TPE-ENT_Cux2", "L5/6_IT_TPE-ENT", "L5/6_IT",
                           "L5/6_NP", "L5/6_NP_CTX", "L6_IT", "L6b_ENT", "L6b/CT_ENT")

exn.primate@active.ident <- factor(Idents(exn.primate), levels = desired_order_primate)

dotplot <- DotPlot(exn.primate, features = EC.markers) +  coord_flip() + theme_minimal() + RotatedAxis() 
p1 <- dotplot + scale_color_gradient2(
  low = "blue",
  mid = "grey",
  high = "red",
  midpoint = 0
) + theme(axis.text.y = element_text(face = "italic", size = 20), axis.text.x = element_text(size = 15))

# INTEGRATED
desired_order_mammal <- c("L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_Reln", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat", "L2/3_IT_PPP_mouse", "L2/3_IT", "L3_IT_ENT",
                          "L4-6", "L5_PPP_mouse", "L5/6_IT_TPE-ENT", "L5/6_NP", "L6b/CT_ENT", "L6_IT_ENTL", "L6b/CT_ENT_mouse", "L6_IT_ENTL_mouse", "L6_IT_primate", "L6b",
                          "L2-6", "L2-6_bat", "L2-6_primate", "CR")
exn.mammal@active.ident <- factor(Idents(exn.mammal), levels = desired_order_mammal)

dotplot <- DotPlot(exn.mammal, features = EC.markers) +  coord_flip() + theme_minimal() + RotatedAxis() 
p4 <- dotplot + scale_color_gradient2(
  low = "blue",
  mid = "grey",
  high = "red",
  midpoint = 0
) + theme(axis.text.y = element_text(face = "italic", size = 20), axis.text.x = element_text(size = 15))

p1 + p2 + p3 + p4

