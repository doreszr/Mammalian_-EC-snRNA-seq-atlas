########################################################################################################################

### Quality control of samples and batch correction of species

########################################################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)

#### Violin plots ####
# Loading raw samples
human1 <- readRDS("Data_EC_atlas/Seurat_obj/human/human1raw.rds")
human2 <- readRDS("Data_EC_atlas/Seurat_obj/human/human2raw.rds")
human3 <- readRDS("Data_EC_atlas/Seurat_obj/human/human3raw.rds")

mouse1 <- readRDS("Data_EC_atlas/Seurat_obj/mouse/mouseraw1.rds")
mouse2 <- readRDS("Data_EC_atlas/Seurat_obj/mouse/mouseraw2.rds")
mouse3 <- readRDS("Data_EC_atlas/Seurat_obj/mouse/mouseraw3.rds")

baboon1 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamedbaboon1raw.rds")
baboon2 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamedbaboon2raw.rds")
baboon3 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamedbaboon3raw.rds")


bat3 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedbat3raw.rds")
bat4 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedbat4raw.rds")
bat5 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedbat5raw.rds")

# Find % mitochondrial reads in each cell
mouse1[["percent.mt"]] <- PercentageFeatureSet(mouse1, pattern = "^mt-")
mouse2[["percent.mt"]] <- PercentageFeatureSet(mouse2, pattern = "^mt-")
mouse3[["percent.mt"]] <- PercentageFeatureSet(mouse3, pattern = "^mt-")

human1[["percent.mt"]] <- PercentageFeatureSet(human1, pattern = "^MT-")
human2[["percent.mt"]] <- PercentageFeatureSet(human2, pattern = "^MT-")
human3[["percent.mt"]] <- PercentageFeatureSet(human3, pattern = "^MT-")

# Merge them into one object
merged_obj <- merge(human1, y = list(human2, human3, mouse1, mouse2, mouse3, baboon1, baboon2, baboon3, bat3, bat4, bat5),
                    add.cell.ids = c("human1", "human2", "human3", "mouse1", "mouse2", "mouse3", "baboon1", "baboon2", "baboon3", "bat3", "bat4", "bat5"))

# Re-join layers - next step
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])

# Now plot using VlnPlot
VlnPlot(merged_obj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), group.by = "orig.ident", raster=FALSE, pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### Samples UMAPs ####
### BABOON

baboon1 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/forplotrenamedbaboon1.rds")
baboon2 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/forplotrenamedbaboon2.rds")
baboon3 <- readRDS("Data_EC_atlas/Seurat_obj/baboon/forplotrenamedbaboon3.rds")

## Visualize singlets
DimPlot(baboon1)
DimPlot(baboon2)
DimPlot(baboon3)


### BAT

bat3 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/forplotrenamedbat3.rds")
bat4 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/forplotrenamedbat4.rds")
bat5 <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/forplotrenamedbat5.rds")


## Visualize singlets
DimPlot(bat3)
DimPlot(bat4)
DimPlot(bat5)

#### Batch correction UMAPs ####
merged.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_annot.rds")
DimPlot(merged.bat, reduction = "umap.rpca", group.by = "orig.ident") + ggtitle("Fruit bats")

merged.baboon <- readRDS("Data_EC_atlas/Seurat_obj/baboon/renamed.merged.baboon_annot.rds")
DimPlot(merged.baboon, reduction = "umap.rpca", group.by = "orig.ident") + ggtitle("Baboons")

merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")
DimPlot(merged.mouse, reduction = "umap.rpca", group.by = "orig.ident") + ggtitle("Mice")

merged.human <- readRDS("Data_EC_atlas/Seurat_obj/human/merged.human_annot.rds")
DimPlot(merged.human, reduction = "umap.rpca", group.by = "orig.ident") + ggtitle("Humans")


merged.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed_primate_annot.rds")
DimPlot(merged.primate, reduction = "umap.rpca", group.by = "species", cols = c("BABOON" = "#92BE36", "HUMAN" = "#927E1B")) + ggtitle("Primates") + NoLegend()

merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")
DimPlot(merged.mammal, reduction = "umap.rpca", group.by = "species") + ggtitle("Mammals")



