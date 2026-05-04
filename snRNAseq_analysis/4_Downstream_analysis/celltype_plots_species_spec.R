################################################################################

# Visualizing cell types in species-specific data sets

################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)

#### Species all cell types UMAPs ####

## BAT

merged.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_annot.rds")

DimPlot(merged.bat, reduction = "umap.rpca", label = F, cols = c("EXN" = "#F8766D", "IN" = "#A3A500",
                                                                 "OPC" = "#00B0F6", "Oligo" = "#E76BF3", "VLMC/Astrocyte/Endothelial" = "#00BFC4")) + NoLegend()

## BABOON

merged.baboon <- readRDS("Data_EC_atlas/Seurat_obj/baboon/enamed.merged.baboon_annot.rds")

DimPlot(merged.baboon, reduction = "umap.rpca", label = F, cols = c("EXN" = "#F8766D", "IN" = "#A3A500", "ASTR" = "#00BA38",
                                                                    "OPC" = "#00B0F6", "Oligo" = "#E76BF3", "Endothelial" = "#00BFC4")) + NoLegend()

## HUMAN

merged.human <- readRDS("Data_EC_atlas/Seurat_obj/human/merged.human_annot.rds")

DimPlot(merged.human, reduction = "umap.rpca", label = F, cols = c("EXN" = "#F8766D", "IN" = "#A3A500", "ASTR" = "#00BA38",
                                                                   "OPC" = "#00B0F6", "Oligo" = "#E76BF3", "Endothelial" = "#00BFC4", "Oligo/ASTR" = "#C77CFF",
                                                                   "Microglia" = "#CD9600")) + NoLegend()


### MOUSE
merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")
DimPlot(merged.mouse, reduction = "umap.rpca", group.by = "class_label")

#### Major cell types distribution bar plot ####

## MOUSE
merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")
table(merged.mouse$class_label)

df_mouse <- data.frame(species = c("mouse", "mouse"),
                          major_cell_types = c("EXN", "IN"),
                          count = c(sum(merged.mouse$class_label == "Glutamatergic"), sum(merged.mouse$class_label == "GABAergic")))



## BABOON
merged.baboon <- readRDS("Data_EC_atlas/Seurat_obj/baoon/renamed.merged.baboon_annot.rds")

df_baboon <- as.data.frame(table(merged.baboon$major_cell_types))
colnames(df_baboon) <- c("major_cell_types", "count")

df_baboon$species <- "baboon"


## BAT
merged.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_annot.rds")


df_bat <- as.data.frame(table(merged.bat$major_cell_types))
colnames(df_bat) <- c("major_cell_types", "count")

df_bat$species <- "bat"


## HUMAN
merged.human <- readRDS("Data_EC_atlas/Seurat_obj/human/merged.human_annot.rds")


df_human <- as.data.frame(table(merged.human$major_cell_types))
colnames(df_human) <- c("major_cell_types", "count")

df_human$species <- "human"


majorcelltypes <- bind_rows(df_mouse, df_baboon, df_bat, df_human)


species_totals <- majorcelltypes %>%
  group_by(species) %>%
  summarise(total_cells = sum(count))

celltypes_totals <- majorcelltypes %>%
  group_by(major_cell_types) %>%
  summarise(total_cells = sum(count))

celltype_colors <- c(
  "EXN" = "#F8766D",
  "Glutamatergic" = "#F8766D",
  "IN" = "#A3A500",
  "GABAergic" = "#A3A500",
  "ASTR" = "#00BA38",
  "OPC" = "#00B0F6",
  "Oligo" = "#E76BF3",
  "Endothelial" = "#00BFC4",
  "VLMC/Astrocyte/Endothelial" = "#00BFC4",
  "Microglia" = "#CD9600",
  "Oligo/ASTR" = "#C77CFF"
)


ggplot(majorcelltypes, aes(x = species, y = count, fill = major_cell_types)) +
  geom_col(stat = "identity", width = 0.9) +
  geom_text(
    data = species_totals,
    mapping = aes(x = species, y = total_cells, label = total_cells),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 9  
  ) +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    title = "Cell Type Composition by Species",
    x = "Species",
    y = "Number of Cells",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),          
    axis.line = element_line(color = "black")) +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) + NoLegend()



#### Species neuronal cell types UMAPs ####
### BAT
merged.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedmerged.bat_neu_annot.rds")

DimPlot(merged.bat, reduction = "umap.rpca", label = F, label.size = 3, group.by = "neuronal_celltype_species", repel = T,
        cols = c("Sst" = "#033f63", "Pvalb" = "#98c1d9", "Pvalb_Vipr2" = "#3d5a80",
                 "Lamp5" = "#b5b682", "Lamp5_Lhx6" = "#28666e", "Vip/Sncg" = "#005A32",
                 "L2/3_IT_ENTL_Maml2" = "#ffcba9", "L2/3_IT_ENTL" = "#ff9b57", "L2_IT_ENTL" = "#ff791e", "L2/3_IT" = "#b34800",
                 "L2/3_IT_ENTL_Cux2" = "#C23A22", "L2/3_IT_Cald1" = "#E54C38", "L2/3_IT_PPP" = "#6B2113", "L2/3_IT_ENTm" = "#972E1A",
                 "L4-6" = "#BC746B", "L6b/CT_ENT" = "#fedc97", "L6b" = "#C7A317", "L2-6" = "#C4967D",
                 "L2/3_L6b" = "#C68E17",  "L4/5_IT" = "#ffb703")) + ggtitle("Fruit bat neurons") + NoLegend()

### PRIMATE
merged.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed_primate_annot.rds")

DimPlot(merged.primate, reduction = "umap.rpca", label = F, label.size = 3, group.by = "neuronal_celltype_species", repel = T,
        cols = c("Sst" = "#033f63", "Sst_Chodl" = "#1085CA", "Pvalb" = "#98c1d9", "Pvalb_Vipr2" = "#3d5a80",
                 "Lamp5" = "#b5b682", "Pax6_Adarb2" = "#7c9885", "Vip" = "#005A32", "Lamp5_Lhx6" = "#28666e",
                 "L2/3_IT_ENT" = "#ffcba9", "L2/3_IT_ENTL" = "#ff9b57", "L2/3_IT" = "#b34800", "L2/3_IT_PPP" = "#6B2113",
                 "L5/6_NP_CTX" = "#40130B", "L6_IT" = "#C68E17", "L6b/CT_ENT" = "#E8A317", "L5/6_NP" = "#FFF880", "L6b_ENT" = "#fedc97",
                 "L5/6_IT_TPE-ENT_Cux2" = "#ffb703", "L5/6_IT_TPE-ENT" = "#C4967D", "L5/6_IT" = "#FBCF3B")) + ggtitle("Primate neurons") + NoLegend()

### MOUSE
merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse.annot.rds")

DimPlot(merged.mouse, reduction = "umap.rpca", label = F, label.size = 3, group.by = "neuronal_celltype_species", repel = T,
        cols = c("Sst" = "#033f63", "Sst_Chodl" = "#1085CA", "Pvalb" = "#98c1d9", "Pvalb_Vipr2" = "#3d5a80",
                 "Lamp5" = "#b5b682", "Sncg" = "#7c9885", "Vip" = "#005A32", "Lamp5_Lhx6" = "#28666e",
                 "L2_IT_ENTm" = "#ffcba9", "L2/3_IT_ENTL" = "#ff9b57", "L2_IT_ENTL" = "#ff791e", "L2/3_IT_PPP" = "#6B2113", "L3_IT_ENT" = "#b34800",
                 "L6_IT_ENTL" = "#40130B", "L5/6_IT_TPE-ENT" = "#ffb703", "CR" = "#663c20", "L6b/CT_ENT" = "#fedc97", "L6b_CTX" = "#C7A317",
                 "NP_PPP" = "#FFF880", "L5_PPP" = "#fc0c00")) + ggtitle("Mouse neurons") + NoLegend()

