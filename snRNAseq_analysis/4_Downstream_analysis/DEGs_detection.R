################################################################################

# DEGs detection, heat map plotting and DEGs proportion

################################################################################


library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

#### Pairwise DEGs detection in integrated data set ####
inputdir <- c("give_your_input_folder")

merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

# changing to right meta data
Idents(merged.mammal) <- merged.mammal$evo_neuronal_celltype

### Marker detection
merged.mammal <- PrepSCTFindMarkers(merged.mammal)

# species pairs for pairwise DEGs detection
pairs <- c("PRIMATE|MOUSE","PRIMATE|BAT","MOUSE|BAT", "MOUSE|PRIMATE", "BAT|PRIMATE", "BAT|MOUSE")
clusts <- c("Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5", "Sncg_Pax6", "Vip_primate",
            "L2/3_IT_ENTL", "L6b/CT_ENT", "L6b", "L2-6", "L5/6_NP", "L2_IT_ENTL","L2/3_IT", "L2-6_primate", "L5/6_IT_TPE-ENT", "L4-6", "L6_IT_primate", "L2/3_IT_Reln")

# vector of cluster names with species
species.clusts <- rownames(table(merged.mammal$evo_neuronal_celltype))
# saving normalized count matrix
scaledata <- rownames(merged.mammal@assays[["SCT"]]@scale.data)

# looping through all the pairs
for (pair in pairs){
  # print(pair)
  for (cluster in clusts){
    # print(pair)
    sps <- strsplit(pair, "|", fixed = TRUE)[[1]]
    sp1 <- sps[1]; sp2 <- sps[2]
    
    
    cls1 <- paste0(sp1, "_", cluster);
    cls2 <- paste0(sp2, "_", cluster)
    
    if (cls1 %in% rownames(table(merged.mammal$evo_neuronal_celltype)) & cls2 %in% rownames(table(merged.mammal$evo_neuronal_celltype))) {
      marker1 <- FindMarkers(merged.mammal, ident.1 = cls1, ident.2 = cls2)
      marker1 <- marker1 %>%
        filter(avg_log2FC > 1) %>%
        filter(pct.1 > 0.25) %>%
        filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
      marker1 <- subset(marker1, gene %in% scaledata)
      cls1 <- gsub(" ", "_", cls1)
      cls1 <- gsub("/", ".", cls1)
      saveRDS(marker1, file = paste0(inputdir, cls1, sp2, "markers", ".rds"))
      
    } else {
      print(paste0(sp1, " or ", sp2, " is not part of ", cluster))
  
    }
    
    
  }
  
}

# vector with all the species
allspecies <- c("PRIMATE", "MOUSE", "BAT")

# joining the pairwise DEGs lists for each species
for (species in allspecies){
  sp <- allspecies[allspecies != species]
  sp1 <- sp[1]; sp2 <- sp[2]
  joined.marker <- setNames(
    lapply(clusts, function(clusts) {
      clusts <- gsub(" ", "_", clusts)
      clusts <- gsub("/", ".", clusts)
      marker1 <- readRDS(file = paste0(inputdir, species, "_", clusts, sp1, "markers", ".rds"))
      marker2 <- readRDS(file = paste0(inputdir, species, "_", clusts, sp2, "markers", ".rds"))
      joined <- bind_rows(marker1, marker2)  %>%
        distinct(gene, .keep_all = TRUE)  %>% arrange(p_val_adj, decreasing = T) %>% slice(1:50)
      
      return(joined)
    }), clusts )
  saveRDS(joined.marker , file = paste0(inputdir, species, "subtype_markers.rds"))
}



#### Species markers ####
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

mouse <- readRDS("inputdir/MOUSEsubtype_markers.rds")

primate <- readRDS("inputdir/PRIMATEsubtype_markers.rds")

bat <- readRDS("inputdir/BATsubtype_markers.rds")


PVALB.markers <- c(primate[["Pvalb"]][["gene"]][1:10], mouse[["Pvalb"]][["gene"]][1:10], bat[["Pvalb"]][["gene"]][1:10])


L2.3ITENTl.markers <- c(primate[["L2/3_IT_ENTL"]][["gene"]][1:10], mouse[["L2/3_IT_ENTL"]][["gene"]][1:10], bat[["L2/3_IT_ENTL"]][["gene"]][1:10])


Idents(merged.mammal) <- merged.mammal$evo_neuronal_celltype
desired_order_evo <- c("PRIMATE_Sst", "PRIMATE_Pvalb", "PRIMATE_Pvalb_Vipr2", "PRIMATE_Lamp5",
                       "PRIMATE_Sncg_Pax6","PRIMATE_Vip_primate",
                       "PRIMATE_L2_IT_ENTL", "PRIMATE_L2/3_IT_ENTL", "PRIMATE_L2/3_IT_ENTL_bat",
                       "PRIMATE_L2/3_IT_ENTL_Maml2_bat", "PRIMATE_L2/3_IT_PPP_mouse", "PRIMATE_L2/3_IT",
                       "PRIMATE_L2/3_IT_Reln", "PRIMATE_L3_IT_ENT","PRIMATE_L4-6", "PRIMATE_L5/6_IT_TPE-ENT",
                       "PRIMATE_L5/6_NP",  "PRIMATE_L6b", "PRIMATE_L6_IT_ENTL_mouse","PRIMATE_L6b/CT_ENT_mouse",
                       "PRIMATE_L6b/CT_ENT", "PRIMATE_L6_IT_primate", "PRIMATE_L2-6_bat", "PRIMATE_L2-6_primate",
                       "PRIMATE_L2-6", "PRIMATE_CR",
                       "MOUSE_Sst_Chodl_mouse","MOUSE_Sst","MOUSE_Pvalb","MOUSE_Pvalb_Vipr2","MOUSE_Lamp5",
                       "MOUSE_Sncg_Pax6","MOUSE_Vip_primate",
                       "MOUSE_L2_IT_ENTL", "MOUSE_L2/3_IT_ENTL", "MOUSE_L2/3_IT_ENTL_bat", "MOUSE_L2/3_IT_ENTL_Maml2_bat",
                       "MOUSE_L2/3_IT_PPP_mouse", "MOUSE_L2/3_IT", "MOUSE_L2/3_IT_Reln", "MOUSE_L3_IT_ENT","MOUSE_L4-6",
                       "MOUSE_L5/6_IT_TPE-ENT",  "MOUSE_L5/6_NP",  "MOUSE_L6b", "MOUSE_L6_IT_ENTL_mouse","MOUSE_L6b/CT_ENT_mouse",
                       "MOUSE_L6b/CT_ENT", "MOUSE_L6_IT_primate", "MOUSE_L5_PPP_mouse","MOUSE_L2-6_bat", "MOUSE_L2-6_primate",
                       "MOUSE_L2-6", "MOUSE_CR",
                       "BAT_Sst","BAT_Pvalb","BAT_Pvalb_Vipr2", "BAT_Lamp5","BAT_Sncg_Pax6", "BAT_Vip_primate",
                       "BAT_L2_IT_ENTL", "BAT_L2/3_IT_ENTL", "BAT_L2/3_IT_ENTL_bat", "BAT_L2/3_IT_ENTL_Maml2_bat",
                       "BAT_L2/3_IT_PPP_mouse", "BAT_L2/3_IT", "BAT_L2/3_IT_Reln", "BAT_L3_IT_ENT", "BAT_L4-6", "BAT_L5/6_IT_TPE-ENT",
                       "BAT_L5/6_NP", "BAT_L6b", "BAT_L6_IT_ENTL_mouse","BAT_L6b/CT_ENT_mouse", "BAT_L6b/CT_ENT", "BAT_L6_IT_primate",
                       "BAT_L2-6_bat", "BAT_L2-6_primate", "BAT_L2-6", "BAT_CR")


# Set identities in the specified order
merged.mammal@active.ident <- factor(Idents(merged.mammal), levels = desired_order_evo)


dotplot <- DotPlot(merged.mammal, features = L2.3ITENTl.markers) +  coord_flip() + theme_minimal() + RotatedAxis() 
dotplot + scale_color_gradientn(colors = c("blue", "grey", "red")) + theme(axis.text.y = element_text(face = "italic", size = 20), axis.text.x = element_text(size = 6))


dotplot <- DotPlot(merged.mammal, features = unique(PVALB.markers)) +  coord_flip() + theme_minimal() + RotatedAxis() 
dotplot + scale_color_gradientn(colors = c("blue", "grey", "red")) + theme(axis.text.y = element_text(face = "italic", size = 20), axis.text.x = element_text(size = 6))



#### Subtype DEGs detection in each species and in integrated data set ####
## PRIMATE ##
# loading data
merged.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed_primate_annot.rds")

inputdir <- c("primate")

# saving normalized count matrix
scaledata <- rownames(merged.primate@assays[["SCT"]]@scale.data)

# Setting up cluster names for DEGs detection
clusters <- c("L5/6_IT_TPE-ENT", "L5/6_IT_TPE-ENT_Cux2", "L2/3_IT_ENTL",
              "L2/3_IT","Vip","Sst","Pvalb","L6b/CT_ENT","Lamp5","L2/3_IT_ENT",
              "Pax6_Adarb2", "L6b_ENT","L2/3_IT_PPP", "Lamp5_Lhx6", "L5/6_NP_CTX",
              "Pvalb_Vipr2", "L5/6_NP", "L5/6_IT", "L6_IT", "Sst_Chodl")
EXN <- c("L2/3_IT_ENTL", "L5/6_IT_TPE-ENT", "L2/3_IT", "L5/6_IT_TPE-ENT_Cux2", "L6b/CT_ENT", "L2/3_IT_ENT",
         "L6b_ENT", "L2/3_IT_PPP", "L5/6_NP_CTX", "L5/6_NP", "L5/6_IT", "L6_IT")
IN <- c("Sst", "Sst_Chodl", "Pvalb", "Pvalb_Vipr2", "Lamp5", "Lamp5_Lhx6", "Vip", "Pax6_Adarb2")

# Switching to right meta data
Idents(merged.primate) <- merged.primate$neuronal_celltype_species

# DEGs detection
merged.primate <- PrepSCTFindMarkers(merged.primate)

# using EXN as background for EXN subtypes vica versa
for (cluster in rownames(table(merged.pirmate$neuronal_celltype_species))){
  if (cluster %in% EXN) {
    background <- EXN[EXN != cluster]
    marker1 <- FindMarkers(merged.pirmate, ident.1 = cluster, ident.2 = background)
  } else {
    background <- IN[IN != cluster]
    marker1 <- FindMarkers(merged.pirmate, ident.1 = cluster, ident.2 = background)
  }
  marker1 <- marker1 %>%
    filter(avg_log2FC > 1) %>%
    filter(pct.1 > 0.25) %>%
    filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
  marker1 <- subset(marker1, gene %in% scaledata) %>% slice(1:50)
  cluster <- gsub(" ", "_", cluster)
  cluster <- gsub("/", ".", cluster)
  saveRDS(marker1, file = paste0(inputdir, "PRIMATE", cluster, "markers", ".rds"))
  
}


# setting up vectors for superficial and deep layers
deeplayer <- c("L5/6_IT_TPE-ENT", "L5/6_IT_TPE-ENT_Cux2", "L6b/CT_ENT", "L6b_ENT", "L5/6_NP_CTX", "L5/6_NP", "L5/6_IT", "L6_IT")
superficiallayer <- c("L2/3_IT_ENTL", "L2/3_IT", "L2/3_IT_ENT", "L2/3_IT_PPP")

# DEGs detection between superficial and deep layers
deeplayer.markers <- FindMarkers(merged.primate, ident.1 = deeplayer, ident.2 = superficiallayer)
deeplayer.markers <- deeplayer.markers %>%
  filter(avg_log2FC > 1) %>%
  filter(pct.1 > 0.25) %>%
  filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
deeplayer.markers <- subset(deeplayer.markers, gene %in% scaledata) %>% slice(1:50)

superficiallayer.markers <- FindMarkers(merged.primate, ident.1 = superficiallayer, ident.2 = deeplayer)
superficiallayer.markers <- superficiallayer.markers %>%
  filter(avg_log2FC > 1) %>%
  filter(pct.1 > 0.25) %>%
  filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
superficiallayer.markers <- subset(superficiallayer.markers, gene %in% scaledata) %>% slice(1:50)

combined <-  list(superficiallayer.markers = superficiallayer.markers, deeplayer.markers = deeplayer.markers)

saveRDS(combined, file = file = paste0(inputdir, "PRIMATE_layermarkers.rds"))


## MOUSE ##
# loading data
merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")

inputdir <- c("mouse")

# saving normalized count matrix
scaledata <- rownames(merged.mouse@assays[["SCT"]]@scale.data)

# Setting up cluster names for DEGs detection
clusters <- c("L5_PPP", "L3_IT_ENT", "L6b/CT_ENT", "L2_IT_ENTL", "L5/6_IT_TPE-ENT", "L2/3_IT_PPP", "L2_IT_ENTm",
              "L6_IT_ENTL", "NP_PPP", "CR", "L6b_CTX", "L2/3_IT_ENTL",
              "Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5", "Lamp5_Lhx6", "Vip", "Sncg", "Sst_Chodl")
EXN <- c("L5_PPP", "L3_IT_ENT", "L6b/CT_ENT", "L2_IT_ENTL", "L5/6_IT_TPE-ENT", "L2/3_IT_PPP", "L2_IT_ENTm",
         "L6_IT_ENTL", "NP_PPP", "CR", "L6b_CTX", "L2/3_IT_ENTL")
IN <- c("Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5", "Lamp5_Lhx6", "Vip", "Sncg", "Sst_Chodl")

# Switching to right meta data
Idents(merged.mouse) <- merged.mouse$neuronal_celltype_species

# DEGs detection
merged.mouse <- PrepSCTFindMarkers(merged.mouse)

# using EXN as background for EXN subtypes vica versa
for (cluster in rownames(table(merged.mouse$neuronal_celltype_species))){
  if (cluster %in% EXN) {
    background <- EXN[EXN != cluster]
    marker1 <- FindMarkers(merged.mouse, ident.1 = cluster, ident.2 = background)
  } else {
    background <- IN[IN != cluster]
    marker1 <- FindMarkers(merged.mouse, ident.1 = cluster, ident.2 = background)
  }
  marker1 <- marker1 %>%
    filter(avg_log2FC > 1) %>%
    filter(pct.1 > 0.25) %>%
    filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
  marker1 <- subset(marker1, gene %in% scaledata) %>% slice(1:50)
  cluster <- gsub(" ", "_", cluster)
  cluster <- gsub("/", ".", cluster)
  saveRDS(marker1, file = paste0(inputdir, "MOUSE", cluster, "markers", ".rds"))
  
}



# setting up vectors for superficial and deep layers
deeplayer <- c("L5_PPP", "L6b/CT_ENT", "L5/6_IT_TPE-ENT", "L6_IT_ENTL", "NP_PPP", "L6b_CTX")
superficiallayer <- c("L3_IT_ENT", "L2_IT_ENTL", "L2/3_IT_PPP", "L2_IT_ENTm", "L2/3_IT_ENTL")

# DEGs detection between superficial and deep layers
deeplayer.markers <- FindMarkers(merged.mouse, ident.1 = deeplayer, ident.2 = superficiallayer)
deeplayer.markers <- deeplayer.markers %>%
  filter(avg_log2FC > 1) %>%
  filter(pct.1 > 0.25) %>%
  filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
deeplayer.markers <- subset(deeplayer.markers, gene %in% scaledata) %>% slice(1:50)

superficiallayer.markers <- FindMarkers(merged.mouse, ident.1 = superficiallayer, ident.2 = deeplayer)
superficiallayer.markers <- superficiallayer.markers %>%
  filter(avg_log2FC > 1) %>%
  filter(pct.1 > 0.25) %>%
  filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
superficiallayer.markers <- subset(superficiallayer.markers, gene %in% scaledata) %>% slice(1:50)

combined <-  list(superficiallayer.markers = superficiallayer.markers, deeplayer.markers = deeplayer.markers)

saveRDS(combined, file = paste0(inputdir, "MOUSE_layermarkers.rds"))


## BAT ##
# loading data
merged.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_annot.rds")

inputdir <- c("bat")

# saving normalized count matrix
scaledata <- rownames(merged.bat@assays[["SCT"]]@scale.data)

# Setting up cluster names for DEGs detection
clusters <- c("L2/3_IT_ENTL", "L4-6", "L2/3_IT_ENTL_Maml2", "L6b/CT_ENT", "L2/3_IT_ENTL_Cux2", "L2-6", "L2_IT_ENTL", "L4/5_IT",       
              "L2/3_L6b", "L2/3_IT_PPP", "L2/3_IT_ENTm", "L2/3_IT", "L6b", "L2/3_IT_Cald1",
              "Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5", "Lamp5_Lhx6", "Vip/Sncg")
EXN <- c("L2/3_IT_ENTL", "L4-6", "L2/3_IT_ENTL_Maml2", "L6b/CT_ENT", "L2/3_IT_ENTL_Cux2", "L2-6", "L2_IT_ENTL", "L4/5 _IT",       
         "L2/3 L6b", "L2/3 IT PPP", "L2/3 IT ENTm", "L2/3 IT", "L6b", "L2/3 IT Cald1")
IN <- c("Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5", "Lamp5_Lhx6", "Vip/Sncg")

# Switching to right meta data
Idents(merged.bat) <- merged.bat$neuronal_celltype_species

# DEGs detection
merged.bat <- PrepSCTFindMarkers(merged.bat)

# using EXN as background for EXN subtypes vica versa
for (cluster in rownames(table(merged.bat$neuronal_celltype_species))){
  if (cluster %in% EXN) {
    background <- EXN[EXN != cluster]
    marker1 <- FindMarkers(merged.bat, ident.1 = cluster, ident.2 = background)
  } else {
    background <- IN[IN != cluster]
    marker1 <- FindMarkers(merged.bat, ident.1 = cluster, ident.2 = background)
  }
  marker1 <- marker1 %>%
    filter(avg_log2FC > 1) %>%
    filter(pct.1 > 0.25) %>%
    filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
  marker1 <- subset(marker1, gene %in% scaledata) %>% slice(1:50)
  cluster <- gsub(" ", "_", cluster)
  cluster <- gsub("/", ".", cluster)
  saveRDS(marker1, file = paste0(inputdir, "BAT", cluster, "markers", ".rds"))
 
}


# setting up vectors for superficial and deep layers
deeplayer <- c("L4-6", "L6b/CT_ENT", "L4/5_IT", "L6b")
superficiallayer <- c("L2/3_IT_ENTL",  "L2/3_IT_ENTL_Maml2", "L2/3_IT_ENTL_Cux2", "L2_IT_ENTL", "L2/3_IT_PPP", "L2/3_IT_ENTm", "L2/3_IT", "L2/3_IT_Cald1")

# DEGs detection between superficial and deep layers
deeplayer.markers <- FindMarkers(merged.bat, ident.1 = deeplayer, ident.2 = superficiallayer)
deeplayer.markers <- deeplayer.markers %>%
  filter(avg_log2FC > 1) %>%
  filter(pct.1 > 0.25) %>%
  filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
deeplayer.markers <- subset(deeplayer.markers, gene %in% scaledata) %>% slice(1:50)

superficiallayer.markers <- FindMarkers(merged.bat, ident.1 = superficiallayer, ident.2 = deeplayer)
superficiallayer.markers <- superficiallayer.markers %>%
  filter(avg_log2FC > 1) %>%
  filter(pct.1 > 0.25) %>%
  filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
superficiallayer.markers <- subset(superficiallayer.markers, gene %in% scaledata) %>% slice(1:50)

combined <-  list(superficiallayer.markers = superficiallayer.markers, deeplayer.markers = deeplayer.markers)

saveRDS(combined, file = paste0(inputdir, "BAT_layermarkers.rds"))


# Load integrated data
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

# Switching to right meta data
Idents(merged.mammal) <- merged.mammal$neuronal_celltype

# Setting up cluster names for DEGs detection
clusters <- c("Vip_primate", "Sncg_Pax6", "Lamp5", "Pvalb_Vipr2", "Pvalb", "Sst", "Sst_Chodl_mouse",
              "L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat", "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
              "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse", "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
              "L2-6_bat", "L2-6_primate", "L2-6", "CR")
EXN <- c("L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat", "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
         "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse", "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
         "L2-6_bat", "L2-6_primate", "L2-6", "CR")
IN <- c("Vip_primate", "Sncg_Pax6", "Lamp5", "Pvalb_Vipr2", "Pvalb", "Sst", "Sst_Chodl_mouse")

# DEGs detection
merged.mammal <- PrepSCTFindMarkers(merged.mammal)

# using EXN as background for EXN subtypes vica versa
for (cluster in clusters){
  if (cluster %in% EXN) {
    background <- EXN[EXN != cluster]
    marker1 <- FindMarkers(merged.mammal, ident.1 = cluster, ident.2 = background)
  } else {
    background <- IN[IN != cluster]
    marker1 <- FindMarkers(merged.mammal, ident.1 = cluster, ident.2 = background)
  }
  marker1 <- marker1 %>%
    filter(avg_log2FC > 1) %>%
    filter(pct.1 > 0.25) %>%
    filter(p_val_adj < 0.01) %>% rownames_to_column(var = "gene")
  marker1 <- subset(marker1, gene %in% scaledata) %>% slice(1:50)
  cluster <- gsub(" ", "_", cluster)
  cluster <- gsub("/", ".", cluster)
  saveRDS(marker1, file = paste0(inputdir, "MAMMAL", cluster, "markers", ".rds"))
  
}


### Final marker gene lists
## selecting top10 DEGs from each EXN subtype ###
top10combined.exn.markers <- lapply(EXN, function(cluster){
  cluster <- gsub(" ", "_", cluster)
  cluster <- gsub("/", ".", cluster)
  markers <- readRDS(paste0("MAMMAL", cluster, "markers.rds"))
  top10 <- markers %>% slice(1:10)
}) %>% setNames(EXN)

saveRDS(top10combined.exn.markers, file = paste0(inputdir, "MAMMAL_top10combined.exn.markers", ".rds"))


# load exn subset
exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")

# switching to right meta data
Idents(exn.mammal) <- exn.mammal$neuronal_celltype



#### Deep layer heat map ####

exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")

top10combined.exn.markers <- readRDS("MAMMAL_top10combined.exn.markers.rds")

deep.layers <- c("L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP",  "L6b", "L6_IT_ENT_mouse", 
                 "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse")

top10combined.deep.markers.vector <- c()
for (cluster in deep.layers) {
  markers <- top10combined.exn.markers[[cluster]][["gene"]]
  names(markers) <- rep(cluster, length(markers))
  top10combined.deep.markers.vector <- c(top10combined.deep.markers.vector, markers)
}

# subsetting for plotting

deep.sub <- subset(exn.mammal, idents = deep.layers)

Idents(deep.sub) <- deep.sub$evo_neuronal_celltype

my_order <- c("PRIMATE_L4-6","MOUSE_L4-6", "BAT_L4-6",
              "PRIMATE_L5/6_IT_TPE-ENT", "MOUSE_L5/6_IT_TPE-ENT", "BAT_L5/6_IT_TPE-ENT",
              "PRIMATE_L5/6_NP", "MOUSE_L5/6_NP", "BAT_L5/6_NP",
              "PRIMATE_L6b", "MOUSE_L6b", "BAT_L6b",
              "PRIMATE_L6_IT_ENTL_mouse","MOUSE_L6_IT_ENTL_mouse",  
              "MOUSE_L6b/CT_ENT_mouse", "BAT_L6b/CT_ENT_mouse",
              "PRIMATE_L6b/CT_ENT", "MOUSE_L6b/CT_ENT", "BAT_L6b/CT_ENT",
              "PRIMATE_L6_IT_primate", "MOUSE_L6_IT_primate", "BAT_L6_IT_primate",
              "PRIMATE_L6b/CT_ENT_mouse","BAT_L6_IT_ENTL_mouse", "MOUSE_L5_PPP_mouse")

deep.sub@active.ident <- factor(x = deep.sub@active.ident, levels = my_order)

# plotting
DoHeatmap(deep.sub, features = top10combined.deep.markers.vector) + NoLegend() + scale_fill_viridis_c()

DoHeatmap(deep.sub, features = top10combined.deep.markers.vector)  + scale_fill_viridis_c()


#### Mixed layers heat map ####
exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")
top10combined.exn.markers <- readRDS("MAMMAL_top10combined.exn.markers.rds")

mix.layers <- c("L2-6_bat", "L2-6_primate", "L2-6")

top10combined.exn.markers
top10combined.mix.markers.vector <- c()
for (cluster in mix.layers) {
  markers <- top10combined.exn.markers[[cluster]][["gene"]]
  names(markers) <- rep(cluster, length(markers))
  top10combined.mix.markers.vector <- c(top10combined.mix.markers.vector, markers)
}


mix.sub <- subset(exn.mammal, idents = mix.layers)

Idents(mix.sub) <- mix.sub$evo_neuronal_celltype

my_order <- c("PRIMATE_L2-6_bat","MOUSE_L2-6_bat", "BAT_L2-6_bat",
              "PRIMATE_L2-6_primate", "MOUSE_L2-6_primate", "BAT_L2-6_primate",
              "PRIMATE_L2-6", "MOUSE_L2-6", "BAT_L2-6")

mix.sub@active.ident <- factor(x = mix.sub@active.ident, levels = my_order)

DoHeatmap(mix.sub, features = top10combined.mix.markers.vector) + NoLegend() + scale_fill_viridis_c()

DoHeatmap(mix.sub, features = top10combined.mix.markers.vector) + scale_fill_viridis_c()




#### Combining DEGs detected in IN subtypes from species data sets and from pairwise ####
# setting up species names
allspecies <- c("PRIMATE", "MOUSE", "BAT")

inputdir <- c("mammal")

# selecting IN clusters
clusts <- c("Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5", "Sncg_Pax6", "Vip_primate")

# taking intersections of DEGs detected in species and in integrated between species
for (species in allspecies){
  markers <- readRDS(file = paste0(inputdir, species, "INmarkers.rds"))
  markerslist <- lapply(clusts, function(clusts) {
    subtype.markers <- markers[[clusts]]
    subtype <- gsub(" ", "_", clusts)
    if (species != "PRIMATE") {
      
      if (species == "MOUSE") {
        inputdir <- c("mouse")
        
        if (subtype == "Sncg_Pax6") {
          
          subtype <- "Sncg"
          spec.markers <- readRDS(file = paste0(inputdir, species, subtype, "markers.rds"))
          common <- semi_join(subtype.markers, spec.markers, by = "gene")
          spec.markers.spec <- anti_join(spec.markers, subtype.markers, by = "gene")
          subtype.markers.spec <- anti_join(subtype.markers, spec.markers, by = "gene")
        } else {
          spec.markers <- readRDS(file = paste0(inputdir, species, subtype, "markers.rds"))
          common <- semi_join(subtype.markers, spec.markers, by = "gene")
          spec.markers.spec <- anti_join(spec.markers, subtype.markers, by = "gene")
          subtype.markers.spec <- anti_join(subtype.markers, spec.markers, by = "gene")
        }
      } else {
        inputdir <- c("bat")
        
        if (subtype == "Sncg_Pax6") {
          subtype <- "Vip.Sncg"
          spec.markers <- readRDS(file = paste0(inputdir, species, subtype, "markers.rds"))
          common <- semi_join(subtype.markers, spec.markers, by = "gene")
          spec.markers.spec <- anti_join(spec.markers, subtype.markers, by = "gene")
          subtype.markers.spec <- anti_join(subtype.markers, spec.markers, by = "gene")
          
        } else {
          if (subtype == "Vip_primate") {
            subtype <- "Vip.Sncg"
            spec.markers <- readRDS(file = paste0(inputdir, species, subtype, "markers.rds"))
            common <- semi_join(subtype.markers, spec.markers, by = "gene")
            spec.markers.spec <- anti_join(spec.markers, subtype.markers, by = "gene")
            subtype.markers.spec <- anti_join(subtype.markers, spec.markers, by = "gene")
          } else {
            spec.markers <- readRDS(file = paste0(inputdir, species, subtype, "markers.rds"))
            common <- semi_join(subtype.markers, spec.markers, by = "gene")
            spec.markers.spec <- anti_join(spec.markers, subtype.markers, by = "gene")
            subtype.markers.spec <- anti_join(subtype.markers, spec.markers, by = "gene")
          }
        }
      }
    } else {
      inputdir <- c("primate")
      
      if (subtype == "Sncg_Pax6") {
        subtype <- "Pax6_Adarb2"
        spec.markers <- readRDS(file = paste0(inputdir, species, subtype, "markers.rds"))
        common <- semi_join(subtype.markers, spec.markers, by = "gene")
        spec.markers.spec <- anti_join(spec.markers, subtype.markers, by = "gene")
        subtype.markers.spec <- anti_join(subtype.markers, spec.markers, by = "gene")
      } else {
        spec.markers <- readRDS(file = paste0(inputdir, species, subtype, "markers.rds"))
        common <- semi_join(subtype.markers, spec.markers, by = "gene")
        spec.markers.spec <- anti_join(spec.markers, subtype.markers, by = "gene")
        subtype.markers.spec <- anti_join(subtype.markers, spec.markers, by = "gene")
      }
    }
    combined <-  list(common = common, spec.markers.spec = spec.markers.spec, subtype.markers.spec = subtype.markers.spec)
 
  }) %>% setNames(., clusts)
  inputdir <- c("mammal")
  saveRDS(markerslist, file = paste0(inputdir, species, "INmarkers_intersects.rds"))
}

# loading combined DEGs lists for each species
primate.in.markers <- readRDS("PRIMATEINmarkers_intersects.rds")
mouse.in.markers <- readRDS("MOUSEINmarkers_intersects.rds")
bat.in.markers <- readRDS("BATINmarkers_intersects.rds")

# combining lists for plotting
combined.list <- list(PRIMATE = primate.in.markers, mouse = mouse.in.markers, bat = bat.in.markers)

# IN clusters
clusts.in <- c("Sst", "Vip_primate", "Lamp5", "Pvalb", "Pvalb_Vipr2", "Sncg_Pax6")

# creating combined lists of DEGs both detected in pairwise and subtype within species
combined.species.spec.markers.in <- lapply(clusts.in, function(cluster) {
  primate <- primate.in.markers[[cluster]][["common"]][["gene"]]
  names(primate) <- rep("PRIMATE", length(primate))
  mouse <- mouse.in.markers[[cluster]][["common"]][["gene"]]
  names(mouse) <- rep("MOUSE", length(mouse))
  bat <- bat.in.markers[[cluster]][["common"]][["gene"]]
  names(bat) <- rep("BAT", length(bat))
  
  c(primate, mouse, bat)
}) %>% setNames(., clusts.in)
saveRDS(combined.species.spec.markers.in, file = "INmarkers_combined_speciesspec.rds")

#### IN subytype species specific heat map ####
in.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalIN_annot.rds")

primate.in.markers <- readRDS("PRIMATEINmarkers_intersects.rds")
mouse.in.markers <- readRDS("MOUSEINmarkers_intersects.rds")
bat.in.markers <- readRDS("BATINmarkers_intersects.rds")

combined.list <- list(PRIMATE = primate.in.markers, mouse = mouse.in.markers, bat = bat.in.markers)

clusts.in <- c("Sst", "Vip_primate", "Lamp5", "Pvalb", "Pvalb_Vipr2", "Sncg_Pax6")

combined.species.spec.markers.in <- lapply(clusts.in, function(cluster) {
  primate <- primate.in.markers[[cluster]][["common"]][["gene"]]
  names(primate) <- rep("PRIMATE", length(primate))
  mouse <- mouse.in.markers[[cluster]][["common"]][["gene"]]
  names(mouse) <- rep("MOUSE", length(mouse))
  bat <- bat.in.markers[[cluster]][["common"]][["gene"]]
  names(bat) <- rep("BAT", length(bat))
  
  c(primate, mouse, bat)
}) %>% setNames(., clusts.in)

Idents(in.mammal) <- in.mammal$evo_neuronal_celltype

combined.combined <- c(combined.species.spec.markers.in[["Sst"]], combined.species.spec.markers.in[["Pvalb"]], combined.species.spec.markers.in[["Pvalb_Vipr2"]],
                       combined.species.spec.markers.in[["Lamp5"]], combined.species.spec.markers.in[["Sncg_Pax6"]], combined.species.spec.markers.in[["Vip_primate"]])
my_order <- c("PRIMATE_Sst", "MOUSE_Sst", "BAT_Sst",
              "PRIMATE_Pvalb", "MOUSE_Pvalb", "BAT_Pvalb",
              "PRIMATE_Pvalb_Vipr2", "MOUSE_Pvalb_Vipr2", "BAT_Pvalb_Vipr2",
              "PRIMATE_Lamp5", "MOUSE_Lamp5", "BAT_Lamp5",
              "PRIMATE_Sncg_Pax6", "MOUSE_Sncg_Pax6", "BAT_Sncg_Pax6",
              "PRIMATE_Vip_primate", "MOUSE_Vip_primate", "BAT_Vip_primate",
              "MOUSE_Sst_Chodl_momuse")

in.mammal@active.ident <- factor(x = in.mammal@active.ident, levels = my_order)


DoHeatmap(in.mammal, features = combined.combined, cells = NULL) + NoLegend() + scale_fill_viridis_c()


#### Species specific subtype marker detection in IN ####
# load subset data
in.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalIN_annot.rds")

# load markers object
combined.species.spec.markers.in <- readRDS("/Users/lbr214/Documents/phd/10X_snRNAseq/sequencing/Genewiz/data/salmon_bat_baboon/countmatrix/integration/mammalian_neu/mammals/rpca/clean/INmarkers_combined_speciesspec.rds")

IN_celltypes <- c("Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5_Lhx6", "Lamp5", "Sncg_Pax6", "Vip_primate")

subtype_sizes <- tibble(
  subtype = c("Sst", "Vip_primate", "Lamp5", "Pvalb", "Pvalb_Vipr2", "Sncg_Pax6"),
  n_cells = c(3343, 2722, 2427, 2879, 562, 1747)
)


df_counts <- imap_dfr(combined.species.spec.markers.in, ~{
  
  tibble(
    subtype = .y,
    species = names(.x)
  ) %>%
    count(subtype, species, name = "gene_count")
  
})

df_wide <- df_counts %>%
  pivot_wider(
    names_from = species,
    values_from = gene_count,
    values_fill = 0
  )


df_norm <- df_wide %>%
  left_join(subtype_sizes, by = "subtype") %>%
  mutate(across(c(BAT, PRIMATE, MOUSE), ~ .x / n_cells)) %>%  # normalize
  pivot_longer(
    cols = c(BAT, PRIMATE, MOUSE),
    names_to = "species",
    values_to = "n_norm"
  )


ggplot(df_norm, aes(x = subtype, y = n_norm, fill = species)) +
  geom_bar(stat = "identity", width = 0.4) +
  labs(
    title = "DEGs identified in Inhibitory Subtypes by Species",
    x = "Subtype",
    y = "DEGs per cell",
    fill = "Species"
  ) +
  theme_minimal(base_size = 15) +
  scale_fill_manual(values = c(
    "PRIMATE" = "#d5e59a",
    "MOUSE" = "#49348c",
    "BAT" = "#a22d71"
  )) +
  theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        legend.text = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        title = element_text(size = 22)
        
  )

