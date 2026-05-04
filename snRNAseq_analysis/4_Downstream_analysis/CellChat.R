########################################################################################################################

### CellChat analysis

########################################################################################################################


library(CellChat)
library(patchwork)
library(viridis)
library(Seurat)
library(ggplot2)
library(tidyverse)

# load object
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")


# subsetting to species
Idents(merged.mammal) <- merged.mammal$evolution

merged.mammal.primate <- subset(merged.mammal, ident = "PRIMATE")
merged.mammal.bat <- subset(merged.mammal, ident = "BAT")
merged.mammal.mouse <- subset(merged.mammal, ident = "MOUSE")

# now find shared clusters with more than 10 cells from each species
Idents(merged.mammal.primate) <- merged.mammal.primate$neuronal_celltype
Idents(merged.mammal.bat) <- merged.mammal.bat$neuronal_celltype
Idents(merged.mammal.mouse) <- merged.mammal.mouse$neuronal_celltype

table(merged.mammal.primate$neuronal_celltype)
table(merged.mammal.bat$neuronal_celltype)


primate.clusts <- names(which(table(merged.mammal.primate$neuronal_celltype) > 10))
bat.clusts <- names(which(table(merged.mammal.bat$neuronal_celltype) > 10))
mouse.clusts <- names(which(table(merged.mammal.mouse$neuronal_celltype) > 10))

shared.clusts <- Reduce(intersect, list(primate.clusts, bat.clusts, mouse.clusts))

# keep shared clusters only
primate.shared <- subset(merged.mammal.primate, 
                         subset = neuronal_celltype %in% shared.clusts)

bat.shared    <- subset(merged.mammal.bat, 
                        subset = neuronal_celltype %in% shared.clusts)

mouse.shared  <- subset(merged.mammal.mouse, 
                        subset = neuronal_celltype %in% shared.clusts)

primate.shared$neuronal_celltype <- droplevels(primate.shared$neuronal_celltype)
bat.shared$neuronal_celltype    <- droplevels(bat.shared$neuronal_celltype)
mouse.shared$neuronal_celltype  <- droplevels(mouse.shared$neuronal_celltype)

# do cell chat on each species separately
# PRIMATE
Idents(primate.shared) <- primate.shared$neuronal_celltype

data.input <- primate.shared[["SCT"]]$data # normalized data matrix

labels <- Idents(primate.shared)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat.primate <- createCellChat(object = primate.shared, group.by = "neuronal_celltype", assay = "SCT")


CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# set the used database in the object
cellchat.primate@DB <- CellChatDB.use
cellchat.primate@DB <- CellChatDB
cellchat.primate@DB <- as.list(CellChatDB.use)  # Convert and assign worked!

# subset the expression data of signaling genes for saving computation cost
cellchat.primate <- subsetData(cellchat.primate) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat.primate <- identifyOverExpressedGenes(cellchat.primate)
cellchat.primate <- identifyOverExpressedInteractions(cellchat.primate)


cellchat.primate <- computeCommunProb(cellchat.primate, type = "triMean")

cellchat.primate <- computeCommunProbPathway(cellchat.primate)

cellchat.primate <- aggregateNet(cellchat.primate)

# BAT
Idents(bat.shared) <- bat.shared$neuronal_celltype

data.input <- bat.shared[["SCT"]]$data # normalized data matrix
labels <- Idents(bat.shared)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat.bat <- createCellChat(object = bat.shared, group.by = "neuronal_celltype", assay = "SCT")


CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# set the used database in the object
cellchat.bat@DB <- CellChatDB.use
cellchat.bat@DB <- CellChatDB
cellchat.bat@DB <- as.list(CellChatDB.use)  # Convert and assign worked!

# subset the expression data of signaling genes for saving computation cost
cellchat.bat <- subsetData(cellchat.bat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat.bat <- identifyOverExpressedGenes(cellchat.bat)
cellchat.bat <- identifyOverExpressedInteractions(cellchat.bat)


cellchat.bat <- computeCommunProb(cellchat.bat, type = "triMean")

cellchat.bat <- computeCommunProbPathway(cellchat.bat)

cellchat.bat <- aggregateNet(cellchat.bat)


# MOUSE
Idents(mouse.shared) <- mouse.shared$neuronal_celltype

data.input <- mouse.shared[["SCT"]]$data # normalized data matrix
labels <- Idents(mouse.shared)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat.mouse <- createCellChat(object = mouse.shared, group.by = "neuronal_celltype", assay = "SCT")


CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# set the used database in the object
cellchat.mouse@DB <- CellChatDB.use
cellchat.mouse@DB <- CellChatDB
cellchat.mouse@DB <- as.list(CellChatDB.use)  # Convert and assign worked!

# subset the expression data of signaling genes for saving computation cost
cellchat.mouse <- subsetData(cellchat.mouse) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat.mouse <- identifyOverExpressedGenes(cellchat.mouse)
cellchat.mouse <- identifyOverExpressedInteractions(cellchat.mouse)


cellchat.mouse <- computeCommunProb(cellchat.mouse, type = "triMean")

cellchat.mouse <- computeCommunProbPathway(cellchat.mouse)

cellchat.mouse <- aggregateNet(cellchat.mouse)

# save cellchat objects
saveRDS(cellchat.mouse, file = "cellchatMOUSE.rds")
saveRDS(cellchat.primate, file = "cellchatPRIMATE.rds")
saveRDS(cellchat.bat, file = "cellchatBAT.rds")


# merging cell chat objects to check differential analysis
cellchat.mouse <- readRDS("cellchatMOUSE.rds")
cellchat.primate <- readRDS("cellchatPRIMATE.rds")
cellchat.bat <- readRDS("cellchatBAT.rds")

### MOUSE PRIMATE
object.list <- list(mouse = cellchat.mouse, primate = cellchat.primate)
cellchat.mp <- mergeCellChat(object.list, add.names = names(object.list))

# SEMA3
pathways.show <- c("SEMA3") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets


pdf("sema3circleMOUSEPRIMATE.pdf", width = 18, height = 9)
par(mfrow = c(1,2), xpd=TRUE, cex = 1.5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


# SLITRK
pathways.show <- c("SLITRK") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("slitrkcircleMOUSEPRIMATE.pdf", width = 18, height = 9)
par(mfrow = c(1,2), xpd=TRUE, cex = 1.5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()



### MOUSE BAT
object.list <- list(mouse = cellchat.mouse, bat = cellchat.bat)
cellchat.mb <- mergeCellChat(object.list, add.names = names(object.list))

# SEMA3
pathways.show <- c("SEMA3") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("sema3circleMOUSEBAT.pdf", width = 18, height = 9)
par(mfrow = c(1,2), xpd=TRUE, cex = 1.5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


# SLITRK
pathways.show <- c("SLITRK") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("slitrkcircleMOUSEBAT.pdf", width = 18, height = 9)
par(mfrow = c(1,2), xpd=TRUE, cex = 1.5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


### BAT PRIMATE
object.list <- list(bat = cellchat.bat, primate = cellchat.primate)
cellchat.bp <- mergeCellChat(object.list, add.names = names(object.list))


# SLITRK
pathways.show <- c("SLITRK") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("slitrkcircleBATPRIMATE.pdf", width = 18, height = 9)
par(mfrow = c(1,2), xpd=TRUE, cex = 1.5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# SEMA3
pathways.show <- c("SEMA3") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("sema3circleBATPRIMATE.pdf", width = 18, height = 9)
par(mfrow = c(1,2), xpd=TRUE, cex = 1.5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()
