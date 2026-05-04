## Script by Dorottya Ralbovszki
## Date: 11.03.2024
## Description: Produce an ortholog table between baboon and mouse where only the genes found in our cell matrices are present.
## On top of the one2one orthologs, one2one mapped many2many/one2many/many2one orthologs are also present.
## Output: one2one mapped ortholog table baboon vs mouse



# load the libraries
library(tidyverse)
library(Seurat)

# load the one2one ortholog data frame
orthombat <- read.csv("orthofinder/mappingorthologs/bat_vs_mouse/output/ortho_mbat_one.csv")

# check for duplicates
sum(duplicated(orthombat$gene_bat))

sum(duplicated(orthombat$gene_m))

orthombat[duplicated(orthombat$gene_m),]

orthombat[duplicated(orthombat$gene_bat),]

# save the duplicates
duplicates <- orthombat[duplicated(orthombat$gene_bat),]
dupli.many <- filter(orthombat,
                     gene_bat %in% duplicates$gene_bat)

orthombat.good <- filter(orthombat,
                       !gene_bat %in% duplicates$gene_bat)


bat3 <- readRDS("bat3raw.rds")
bat4 <- readRDS("bat4raw.rds")
bat5 <- readRDS("bat5raw.rds")
merged.bat <- readRDS("merged.bat.rds")


bat.genestest3 <- as.data.frame(rownames(bat3@assays[["RNA"]]))
bat.genestest4 <- as.data.frame(rownames(bat4@assays[["RNA"]]))
bat.genestest5 <- as.data.frame(rownames(bat5@assays[["RNA"]]))

colnames(bat.genestest3) <- "gene_bat"
colnames(bat.genestest4) <- "gene_bat"
colnames(bat.genestest5) <- "gene_bat"


bat.genestest <- full_join(bat.genestest3, bat.genestest4, bat.genestest5, by = "gene_bat")


mouse1 <- readRDS("mouseraw1.rds")
mouse2 <- readRDS("mouseraw2.rds")
mouse3 <- readRDS("mouseraw3.rds")
merged.mouse <- readRDS("merged_mouse.rds")
genes.mouse1 <- as.data.frame(rownames(mouse1@assays[["RNA"]]))
genes.mouse2 <- as.data.frame(rownames(mouse2@assays[["RNA"]]))
genes.mouse3 <- as.data.frame(rownames(mouse3@assays[["RNA"]]))
colnames(genes.mouse1) <- "gene_m"
colnames(genes.mouse2) <- "gene_m"
colnames(genes.mouse3) <- "gene_m"

all.mouse.genes <- full_join(genes.mouse1, genes.mouse2, genes.mouse3,  by = "gene_m")




one2onembat <- filter(orthombat.good,
                    gene_bat %in% bat.genestest$gene_bat)

many2manyortho <- read.csv("/orthofinder/mappingorthologs/bat_vs_mouse/output/ortho_mbat_many.csv")
many2manyortho <- bind_rows(dupli.many, many2manyortho)



# fix bat duplicates
many2manyortho <- many2manyortho[,2:4]
many2manyortho.good <- distinct(many2manyortho)

# filter out genes that are not present. in our data
many2manybat <- filter(many2manyortho.good,
                       gene_bat %in% bat.genestest$gene_bat)
many2manymouse <- filter(many2manybat,
                         gene_m %in% all.mouse.genes$gene_m)

# keep only ones that became one2one
duplimany <- many2manymouse[duplicated(many2manymouse$gene_bat),]
many2manymouse.good <- filter(many2manymouse,
                              !gene_bat%in% duplimany$gene_bat)
duplimany <- many2manymouse[duplicated(many2manymouse$gene_m),]
many2manymouse.good <- filter(many2manymouse.good,
                              !gene_m %in% duplimany$gene_m)


## check expression in mouse of one2many/many2many
sort_manybat <- many2manymouse[duplicated(many2manymouse$gene_bat),]
sort_manybat <- filter(many2manymouse,
                       gene_bat %in% sort_manybat$gene_bat)

# get raw gene expression
explist <- c()
for (row in 1:nrow(sort_manybat)) {
  gene_m <- sort_manybat[row, 3]
  # print(gene_b)
  exp <- AverageExpression(merged.mouse, assays = "SCT", layer = "data", features = gene_m)
  exp <- unname(exp$SCT[1,])
  if (is.null(exp)) {
    exp <- 0
    explist <- append(explist, exp)
  } else{
    explist <- append(explist, exp)
  }
  
}

# add gene expression to data frame
sort_manybat$gene_mexp <- explist

# pick the gene with the highest expression in our data
sort_manybat_new <- sort_manybat %>% arrange(desc(gene_mexp))
sort_manybat_newClean <- sort_manybat_new[!duplicated(sort_manybat_new$gene_bat),]



## check expression in baboon of one2many/many2many
sort_manymouse <- many2manybat[duplicated(many2manybat$gene_m),]
sort_manymouse <- filter(many2manybat,
                         gene_m %in% sort_manymouse$gene_m)

# get raw gene expression
# (had to modify loop because baboon5 sample is not part of the combined object genes
# that are only expressed in that sample are missing from combined)
explist <- c()
for (row in 1:nrow(sort_manymouse)) {
  gene_bat <- sort_manymouse[row, 1]
  # print(gene_b)
  exp <- AverageExpression(merged.bat, assays = "SCT", features = gene_bat)
  exp <- unname(exp$SCT[1,]) 
  if (is.null(exp)) {
    exp <- 0
    explist <- append(explist, exp)
  } else{
    explist <- append(explist, exp)
  }
  
  
}

# add gene expression to data frame
sort_manymouse$gene_bexp <- explist

# pick the gene with the highest expression in our data
sort_manymouse_new <- sort_manymouse %>% arrange(desc(gene_bexp))
sort_manymouse_newClean <- sort_manymouse_new[!duplicated(sort_manymouse_new$gene_m),]

# join the two data frames (baboon expression and human expression)
combimbmany <- full_join(sort_manymouse_newClean, sort_manybat_newClean, by = c("gene_bat", "geneid_m", "gene_m"), na_matches = "na")

# fix many2many orthologs
for (row in 1:nrow(combimbmany)) {
  gene_bat <- combimbmany[row, 1]
  gene_m <- combimbmany[row, 3]
  bexp <- combimbmany[row, 4]
  mexp <- combimbmany[row, 5]
  if(is.na(bexp) == 'TRUE'){
    exp <- AverageExpression(merged.bat, assays = "SCT", features = gene_bat)
    if (length(exp) == 0){
      combimbmany[row, 4] <- "missing"
    } else {
      combimbmany[row, 4] <- unname(exp$SCT[1,])
    } 
  } else {
    exp <- AverageExpression(merged.mouse, assays = "SCT", features = gene_m)
    if (length(exp) == 0){
      combimbmany[row, 5] <- "missing"
    } else {
      combimbmany[row, 5] <- unname(exp$SCT[1,])
    }
  }
}

# pick gene with highest expression
combimbmany_new <- combimbmany %>% arrange(desc(gene_bexp))
combimbmanyCleanm <- combimbmany_new[!duplicated(combimbmany_new$gene_m),]

combimbmany_new <- combimbmanyCleanm %>% arrange(desc(gene_mexp))
combimbmanyCleanb <- combimbmany_new[!duplicated(combimbmany_new$gene_bat),]
sum(duplicated(combimbmanyCleanb$gene_bat)) # 0

# remove columns we don't need
many2manybm <- combimbmanyCleanb[,1:3]
one2onembat <- one2onembat[2:4]

# join one2one and manys
allorthobm <- bind_rows(many2manybm, one2onembat, many2manymouse.good)

# check for duplicates
sum(duplicated(allorthobm$gene_bat))
sum(duplicated(allorthobm$gene_m))
allorthobm[duplicated(allorthobm$gene_bat),]
allorthobm[duplicated(allorthobm$gene_m),]

# LOC105297012 ENSMUSG00000047307 Pcdhb13
AverageExpression(merged.mouse, assays = "SCT", features = "Pcdhb13")
AverageExpression(merged.mouse, assays = "SCT", features = "Pcdhb7")

allorthobmtest <- filter(allorthobm,
                         geneid_m != "ENSMUSG00000047307")

# LOC105295323 ENSMUSG00000090053   Pakap
AverageExpression(merged.mouse, assays = "SCT", features = "Pakap")
AverageExpression(merged.mouse, assays = "SCT", features = "Pakap")
allorthobmtest <- filter(allorthobmtest,
                         geneid_m != "ENSMUSG00000090053")

# ERCC6 ENSMUSG00000054051 Ercc6
AverageExpression(merged.mouse, assays = "SCT", features = "Zfp14")
AverageExpression(merged.mouse, assays = "SCT", features = "Zfp82")
allorthobmtest <- filter(allorthobmtest,
                         geneid_m != "ENSMUSG00000098022")

# UBA52 ENSMUSG00000068240 Uba52rt
AverageExpression(merged.mouse, assays = "SCT", features = "Uba52rt") # not fount
AverageExpression(merged.mouse, assays = "SCT", features = "Uba52")
allorthobmtest <- filter(allorthobmtest,
                         geneid_m != "ENSMUSG00000068240")


# check duplicates
sum(duplicated(allorthobmtest$gene_bat)) # 0
sum(duplicated(allorthobmtest$gene_m)) # 0
allorthobmtest[duplicated(allorthobmtest$gene_m),]


# final data frame
allorthobm <- allorthobmtest
# save final data frame
write_csv2(allorthobm, "orthofinder/mappingorthologs/bat_vs_mouse/output/allorthobatm.csv")

