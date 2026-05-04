## Script by Dorottya Ralbovszki
## Date: 23.09.2024
## Description: Produce an ortholog table between human and mouse where only the genes found in our cell matrices are present.
## On top of the one2one orthologs, one2one mapped many2many/one2many/many2one orthologs are also present.
## Output: one2one mapped ortholog table baboon vs mouse



# load the libraries
library(tidyverse)
library(Seurat)

# load the one2one ortholog data frame
orthomh <- read.csv("orthofinder/mappingorthologs/human_vs_mouse/output/ortho_mh_one.csv")

# check for duplicates
sum(duplicated(orthomh$gene_h))

sum(duplicated(orthomh$gene_m))

orthomh[duplicated(orthomh$gene_m),]

orthomh[duplicated(orthomh$gene_h),]

# save the duplicates
duplicates <- orthomh[duplicated(orthomh$gene_m),]
dupli.many <- filter(orthomh,
                     gene_m %in% duplicates$gene_m)

orthomh.good <- filter(orthomh,
                       !gene_m %in% duplicates$gene_m)

# load the singe-cell matrices of baboon samples (combined as well for gene expression levels)
combined.human <- readRDS("merged.human.rds")

human1 <- readRDS("human1raw.rds")
human2 <- readRDS("human2raw.rds")
human3 <- readRDS("human3raw.rds")

# get the genes from the raw count assay
human.genes <- as.data.frame(rownames(combined.human@assays[["RNA"]]))
genes.human1 <- as.data.frame(rownames(human1@assays[["RNA"]]))
genes.human2 <- as.data.frame(rownames(human2@assays[["RNA"]]))
genes.human3 <- as.data.frame(rownames(human3@assays[["RNA"]]))
# set column name to match ortholog data frame
colnames(genes.human1) <- "gene_h"
colnames(genes.human2) <- "gene_h"
colnames(genes.human3) <- "gene_h"
colnames(human.genes) <- "gene_h"
# join the gene lists from all samples
test.human.genes <- full_join(genes.human1, genes.human2, genes.human3,  by = "gene_h")

# load the mouse single-cell matrices the same way
combined.mouse <- readRDS("merged_mouse.rds")
mouse.genes <- as.data.frame(rownames(combined.mouse@assays[["RNA"]]))
mouse1 <- readRDS("mouseraw1.rds")
mouse2 <- readRDS("mouseraw2.rds")
mouse3 <- readRDS("mouseraw3.rds")
genes.mouse1 <- as.data.frame(rownames(mouse1@assays[["RNA"]]))
genes.mouse2 <- as.data.frame(rownames(mouse2@assays[["RNA"]]))
genes.mouse3 <- as.data.frame(rownames(mouse3@assays[["RNA"]]))
colnames(genes.mouse1) <- "gene_m"
colnames(genes.mouse2) <- "gene_m"
colnames(genes.mouse3) <- "gene_m"
colnames(mouse.genes) <- "gene_m"
all.mouse.genes <- full_join(genes.mouse1, genes.mouse2, genes.mouse3,  by = "gene_m")



# keep only the one2one orthologs that are present in our baboon cell matrices
one2onemh <- filter(orthomh.good,
                    gene_h %in% test.human.genes$gene_h)

many2manyortho <- read.csv("orthofinder/mappingorthologs/human_vs_mouse/output/ortho_mh_many.csv")
many2manyortho <- bind_rows(dupli.many, many2manyortho)

# fix mouse duplicates
many2manyortho <- many2manyortho[,2:4]
many2manyortho.good <- distinct(many2manyortho)

# filter out genes that are not present. in our data
many2manyhuman <- filter(many2manyortho.good,
                       gene_h %in% test.human.genes$gene_h)
many2manymouse <- filter(many2manyhuman,
                         gene_m %in% all.mouse.genes$gene_m)

# keep only ones that became one2one
duplimany <- many2manymouse[duplicated(many2manymouse$gene_h),]
many2manymouse.good <- filter(many2manymouse,
                              !gene_h %in% duplimany$gene_h)
duplimany <- many2manymouse[duplicated(many2manymouse$gene_m),]
many2manymouse.good <- filter(many2manymouse.good,
                              !gene_m %in% duplimany$gene_m)


## check expression in mouse of one2many/many2many
sort_manyhuman <- many2manymouse[duplicated(many2manymouse$gene_h),]
sort_manyhuman <- filter(many2manymouse,
                       gene_h %in% sort_manyhuman$gene_h)

explist <- c()
for (row in 1:nrow(sort_manyhuman)) {
  gene_m <- sort_manyhuman[row, 3]
  # print(gene_b)
  exp <- AverageExpression(combined.mouse, assays = "SCT", layer = "data", features = gene_m)
  exp <- unname(exp$SCT[1,])
  if (is.null(exp)) {
    exp <- 0
    explist <- append(explist, exp)
  } else{
    explist <- append(explist, exp)
  }
  
}


# add gene expression to data frame
sort_manyhuman$gene_mexp <- explist

# pick the gene with the highest expression in our data
sort_manyhuman_new <- sort_manyhuman %>% arrange(desc(gene_mexp))
sort_manyhuman_newClean <- sort_manyhuman_new[!duplicated(sort_manyhuman_new$gene_h),]



## check expression in baboon of one2many/many2many
sort_manymouse <- many2manyhuman[duplicated(many2manyhuman$gene_m),]
sort_manymouse <- filter(many2manyhuman,
                         gene_m %in% sort_manymouse$gene_m)

# get raw gene expression
# (had to modify loop because baboon5 sample is not part of the combined object genes
# that are only expressed in that sample are missing from combined)
explist <- c()
for (row in 1:nrow(sort_manymouse)) {
  gene_h <- sort_manymouse[row, 1]
  # print(gene_b)
  exp <- AverageExpression(combined.human, assays = "SCT", features = gene_h)
  exp <- unname(exp$SCT[1,]) 
  if (is.null(exp)) {
    exp <- 0
    explist <- append(explist, exp)
  } else{
    explist <- append(explist, exp)
  }
  
  
}

# add gene expression to data frame
sort_manymouse$gene_hexp <- explist

# pick the gene with the highest expression in our data
sort_manymouse_new <- sort_manymouse %>% arrange(desc(gene_hexp))
sort_manymouse_newClean <- sort_manymouse_new[!duplicated(sort_manymouse_new$gene_m),]

# join the two data frames (baboon expression and human expression)
combimhmany <- full_join(sort_manymouse_newClean, sort_manyhuman_newClean, by = c("gene_h", "geneid_m", "gene_m"), na_matches = "na")

# fix many2many orthologs
for (row in 1:nrow(combimhmany)) {
  gene_h <- combimhmany[row, 1]
  gene_m <- combimhmany[row, 3]
  hexp <- combimhmany[row, 4]
  mexp <- combimhmany[row, 5]
  if(is.na(hexp) == 'TRUE'){
    exp <- AverageExpression(combined.human, assays = "SCT", features = gene_h)
    if (length(exp) == 0){
      combimhmany[row, 4] <- "missing"
    } else {
      combimhmany[row, 4] <- unname(exp$SCT[1,])
    } 
  } else {
    exp <- AverageExpression(combined.mouse, assays = "SCT", features = gene_m)
    if (length(exp) == 0){
      combimhmany[row, 5] <- "missing"
    } else {
      combimhmany[row, 5] <- unname(exp$SCT[1,])
    }
  }
}

# pick gene with highest expression
combimhmany_new <- combimhmany %>% arrange(desc(gene_hexp))
combimhmanyCleanm <- combimhmany_new[!duplicated(combimhmany_new$gene_m),]

combimhmany_new <- combimhmanyCleanm %>% arrange(desc(gene_mexp))
combimhmanyCleanh <- combimhmany_new[!duplicated(combimhmany_new$gene_h),]
sum(duplicated(combimhmanyCleanh$gene_h)) # 0

# remove columns we don't need
many2manyhm <- combimhmanyCleanh[,1:3]
one2onehm <- one2onemh[2:4]

# join one2one and manys
allorthohm <- bind_rows(many2manyhm, one2onehm, many2manymouse.good)

# check for duplicates
sum(duplicated(allorthohm$gene_h))
sum(duplicated(allorthohm$gene_m))
allorthohm[duplicated(allorthohm$gene_h),]
allorthohm[duplicated(allorthohm$gene_m),]

# PINX1 ENSMUSG00000021958  Pinx1
AverageExpression(combined.mouse, assays = "SCT", features = "Pinx1")
AverageExpression(combined.mouse, assays = "SCT", features = "Sox7")

allorthohmtest <- filter(allorthohm,
                         geneid_m != "ENSMUSG00000021958")



# check duplicates
sum(duplicated(allorthohmtest$gene_h)) # 0
sum(duplicated(allorthohmtest$gene_m)) # 1
allorthohmtest[duplicated(allorthohmtest$gene_m),]

# CXorf58 ENSMUSG00000118554 Fam90a1b
AverageExpression(combined.human, assays = "SCT", features = "CXorf58")
AverageExpression(combined.human, assays = "SCT", features = "FAM90A1")
allorthohmtest <- filter(allorthohmtest,
                         gene_h != "CXorf58")


# double check for duplicates
sum(duplicated(allorthohmtest$gene_h))
sum(duplicated(allorthohmtest$gene_m))

# final data frame
allorthohm <- allorthohmtest
# save final data frame
write_csv2(allorthohm, "orthofinder/mappingorthologs/human_vs_mouse/output/allorthohm.csv")

