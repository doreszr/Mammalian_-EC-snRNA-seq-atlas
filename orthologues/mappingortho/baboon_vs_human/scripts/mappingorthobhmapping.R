## Script by Dorottya Ralbovszki
## Date: 19.02.2024
## Description: Produce an ortholog table between baboon and human where only the genes found in our cell matrices are present.
## On top of the one2one orthologs, one2one mapped many2many/one2many/many2one orthologs are also present.
## Output: one2one mapped ortholog table baboon vs human


# load the libraries
library(tidyverse)
library(Seurat)

# load the one2one ortholog data frame
orthohb <- read.csv("orthofinder/mappingorthologs/baboon_vs_human/output/ortho_hb_one.csv")

# check for duplicates
sum(duplicated(orthohb$gene_b))

sum(duplicated(orthohb$gene_h))

orthohb[duplicated(orthohb$gene_h),]

orthohb[duplicated(orthohb$gene_b),]

# save the duplicates
duplicates <- orthohb[duplicated(orthohb$gene_b),]
dupli.many <- filter(orthohb,
                     gene_b %in% duplicates$gene_b)

# clean up the one2one data frame from those duplicates
orthohb.good <- filter(orthohb,
                       !gene_b %in% duplicates$gene_b)

# load the singe-cell matrices of baboon samples (combined as well for gene expression levels)
combined.baboon <- readRDS("merged.baboon.rds")
bab.genes <- as.data.frame(rownames(combined.baboon@assays[["RNA"]]))
baboon1 <- readRDS("rawbaboon1.rds")
baboon2 <- readRDS("rawbaboon2.rds")
baboon3 <- readRDS("rawbaboon3.rds")
baboon5 <- readRDS("rawbaboon5.rds")

# get the genes from the raw count assay
bab.genestest1 <- as.data.frame(rownames(baboon1@assays[["RNA"]]))
bab.genestest2 <- as.data.frame(rownames(baboon2@assays[["RNA"]]))
bab.genestest3 <- as.data.frame(rownames(baboon3@assays[["RNA"]]))
bab.genestest5 <- as.data.frame(rownames(baboon5@assays[["RNA"]]))
# set column name to match ortholog data frame
colnames(bab.genestest1) <- "gene_b"
colnames(bab.genestest2) <- "gene_b"
colnames(bab.genestest3) <- "gene_b"
colnames(bab.genestest5) <- "gene_b"
colnames(bab.genes) <- "gene_b"
# join the gene lists from all samples
bab.genestest <- full_join(bab.genestest1, bab.genestest2,  by = "gene_b")
bab.genestest <- full_join(bab.genestest, bab.genestest3, by = "gene_b")
bab.genestest <- full_join(bab.genestest, bab.genestest5, by = "gene_b")

# load the human single-cell matrices the same way
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



# keep only the one2one orthologs that are present in our baboon cell matrices
one2onehb <- filter(orthohb.good,
                    gene_b %in% bab.genestest$gene_b)


# load the many2many/many2one/one2many orthologs
many2manyortho <- read.csv("orthofinder/mappingorthologs/baboon_vs_human/output/ortho_hb_many.csv")
many2manyortho <- bind_rows(dupli.many, many2manyortho)

# explore for discrepancies
weirddupli <- filter(many2manyortho,
                     gene_b == gene_h)
weirddupli.gene <- weirddupli$gene_h
sum(duplicated(weirddupli.gene))
weirddupli.gene[duplicated(weirddupli.gene)]

many2manyortho.dupli <- filter(many2manyortho,
                         gene_h %in% weirddupli.gene)

leak.genes <- filter(many2manyortho.dupli,
                    !gene_b == gene_h)

# fix baboon duplicates
many2manyortho <- many2manyortho[,2:4]
many2manyortho.good <- distinct(many2manyortho)

# filter out genes that are not present in our cell matrices (human and baboon)
many2manybab <- filter(many2manyortho.good,
                       gene_b %in% bab.genestest$gene_b)
many2manyhuman <- filter(many2manybab,
                         gene_h %in% test.human.genes$gene_h)

# keep only ones that became one2one
duplimany <- many2manyhuman[duplicated(many2manyhuman$gene_b),]
many2manyhuman.good <- filter(many2manyhuman,
                              !gene_b %in% duplimany$gene_b)
duplimany <- many2manyhuman[duplicated(many2manyhuman$gene_h),]
many2manyhuman.good <- filter(many2manyhuman.good,
                              !gene_h %in% duplimany$gene_h)


## check expression in human of one2many/many2many/many2one
sort_manybab <- many2manyhuman[duplicated(many2manyhuman$gene_b),]
sort_manybab <- filter(many2manyhuman,
                       gene_b %in% sort_manybab$gene_b)

# get raw gene expression
explist <- c()
for (row in 1:nrow(sort_manybab)) {
  gene_h <- sort_manybab[row, 3]
  # print(gene_b)
  exp <- AverageExpression(combined.human, assays = "SCT", layer = "data", features = gene_h)
  exp <- unname(exp$SCT[1,])
  if (is.null(exp)) {
    exp <- 0
    explist <- append(explist, exp)
  } else{
    explist <- append(explist, exp)
  }
  
}

# add gene expression to data frame
sort_manybab$gene_hexp <- explist

# pick the gene with the highest expression in our data
sort_manybab_new <- sort_manybab %>% arrange(desc(gene_hexp))
sort_manybab_newClean <- sort_manybab_new[!duplicated(sort_manybab_new$gene_b),]



## check expression in baboon of one2many/many2many/many2one
sort_manyhuman <- many2manybab[duplicated(many2manybab$gene_h),]
sort_manyhuman <- filter(many2manybab,
                         gene_h %in% sort_manyhuman$gene_h)

# get raw gene expression
# (had to modify loop because baboon5 sample is not part of the combined object genes
# that are only expressed in that sample are missing from combined)
explist <- c()
for (row in 1:nrow(sort_manyhuman)) {
  gene_bab <- sort_manyhuman[row, 1]
  # print(gene_b)
  exp <- AverageExpression(combined.baboon, assays = "SCT", features = gene_bab)
  exp <- unname(exp$SCT[1,]) 
  if (is.null(exp)) {
    exp <- 0
    explist <- append(explist, exp)
  } else{
    explist <- append(explist, exp)
  }
  
  
}

# add gene expression to data frame
sort_manyhuman$gene_bexp <- explist

# pick the gene with the highest expression in our data
sort_manyhuman_new <- sort_manyhuman %>% arrange(desc(gene_bexp))
sort_manyhuman_newClean <- sort_manyhuman_new[!duplicated(sort_manyhuman_new$gene_h),]

# join the two data frames (baboon expression and human expression)
combihbmany <- full_join(sort_manyhuman_newClean, sort_manybab_newClean, by = c("gene_b", "geneid_h", "gene_h"), na_matches = "na")

# fix many2many orthologs
for (row in 1:nrow(combihbmany)) {
  gene_bab <- combihbmany[row, 1]
  gene_h <- combihbmany[row, 3]
  bexp <- combihbmany[row, 4]
  hexp <- combihbmany[row, 5]
  if(is.na(bexp) == 'TRUE'){
    exp <- AverageExpression(combined.baboon, assays = "SCT", features = gene_bab)
    if (length(exp) == 0){
      combihbmany[row, 4] <- "missing"
    } else {
      combihbmany[row, 4] <- unname(exp$SCT[1,])
    } 
  } else {
    exp <- AverageExpression(combined.human, assays = "SCT", features = gene_h)
    if (length(exp) == 0){
      combihbmany[row, 5] <- "missing"
    } else {
      combihbmany[row, 5] <- unname(exp$SCT[1,])
    }
  }
}


# pick gene with highest expression
combihbmany_new <- combihbmany %>% arrange(desc(gene_bexp))
combihbmanyCleanh <- combihbmany_new[!duplicated(combihbmany_new$gene_h),]

combihbmany_new <- combihbmanyCleanh %>% arrange(desc(gene_hexp))
combihbmanyCleanb <- combihbmany_new[!duplicated(combihbmany_new$gene_b),]
sum(duplicated(combihbmanyCleanb$gene_b)) # 0

# remove columns we don't need
many2manybh <- combihbmanyCleanb[,1:3]
one2onebh <- one2onehb[2:4]

# join one2one and manys
allorthobh <- bind_rows(many2manybh, one2onebh, many2manyhuman.good)

# check for duplicates
sum(duplicated(allorthobh$gene_b))
sum(duplicated(allorthobh$gene_h))
allorthobh[duplicated(allorthobh$gene_b),]
allorthobh[duplicated(allorthobh$gene_h),]

# LOC101000556 ENSG00000153201     RANBP2
AverageExpression(combined.human, assays = "SCT", features = "RGPD5")
AverageExpression(combined.human, assays = "SCT", features = "RANBP2")
allorthobhtest <- filter(allorthobh,
                         geneid_h != "ENSG00000015568")

# CALCA ENSG00000110680      CALCA
AverageExpression(combined.human, assays = "SCT", features = "CALCA")
AverageExpression(combined.human, assays = "SCT", features = "CALCB")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000175868")

# LOC101019904 ENSG00000270757 HSPE1-MOB4
AverageExpression(combined.human, assays = "SCT", features = "HSPE1-MOB4")
AverageExpression(combined.human, assays = "SCT", features = "HSPE1")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000270757")

# LOC101003803 ENSG00000177462      OR2T8
AverageExpression(combined.human, assays = "SCT", features = "OR2T8")
AverageExpression(combined.human, assays = "SCT", features = "OR2AK2")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000177462")

# LOC101018472 ENSG00000178685     PARP10
AverageExpression(combined.human, assays = "SCT", features = "PARP10")
AverageExpression(combined.human, assays = "SCT", features = "PLEC")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000178685")

# LOC101022334 ENSG00000131042     LILRB2
AverageExpression(combined.human, assays = "SCT", features = "LILRB2")
AverageExpression(combined.human, assays = "SCT", features = "LILRB4")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000131042")


# TRIM39 ENSG00000204599     TRIM39
AverageExpression(combined.human, assays = "SCT", features = "TRIM39")
AverageExpression(combined.human, assays = "SCT", features = "TRIM39-RPP21")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000248167")


# LOC101005875 ENSG00000025708       TYMP
AverageExpression(combined.human, assays = "SCT", features = "TYMP")
AverageExpression(combined.human, assays = "SCT", features = "SCO2")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000025708")


#  WASHC1 ENSG00000181404     WASHC1
AverageExpression(combined.human, assays = "SCT", features = "WASHC1")
AverageExpression(combined.human, assays = "SCT", features = "WASH6P")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000182484")


#  LOC100997559 ENSG00000269404       SPIB
AverageExpression(combined.human, assays = "SCT", features = "SPIB")
AverageExpression(combined.human, assays = "SCT", features = "POLD1")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000269404")


#  LOC101016189 ENSG00000186314    PRELID2
AverageExpression(combined.human, assays = "SCT", features = "PRELID2")
AverageExpression(combined.human, assays = "SCT", features = "GRXCR2")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000204928")



#  LOC101002450 ENSG00000188629     ZNF177
AverageExpression(combined.human, assays = "SCT", features = "ZNF177")
AverageExpression(combined.human, assays = "SCT", features = "ZNF559")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000188629")


# check duplicates
sum(duplicated(allorthobhtest$gene_b))
sum(duplicated(allorthobhtest$gene_h))
allorthobhtest[duplicated(allorthobhtest$gene_h),]


# fix duplicates
# LOC101025445 ENSG00000015479  MATR3
AverageExpression(combined.baboon, assays = "SCT", features = "MATR3")
AverageExpression(combined.baboon, assays = "SCT", features = "LOC101025445")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000015479")

# PINX1 ENSG00000254093  PINX1
AverageExpression(combined.baboon, assays = "SCT", features = "PINX1")
AverageExpression(combined.baboon, assays = "SCT", features = "SOX7")
allorthobhtest <- filter(allorthobhtest,
                         geneid_h != "ENSG00000258724")



# double chec for duplicates
sum(duplicated(allorthobhtest$gene_b))
sum(duplicated(allorthobhtest$gene_h))

# final data frame
allorthobh <- allorthobhtest
# save final data frame
write_csv2(allorthobh, "orthofinder/mappingorthologs/baboon_vs_human/output/allorthobh.csv")





