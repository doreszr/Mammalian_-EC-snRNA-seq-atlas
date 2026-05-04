################################################################################

# Pearson correlation from DOI: 10.1126/science.abo7257

################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)


## Load data set and set up path

inputdir <- "your_path"

merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

Idents(merged.mammal) <- merged.mammal$neuronal_celltype
merged.mammal <- RenameIdents(merged.mammal, "Vip_primate" = "IN", "Sncg_Pax6" = "IN", "Lamp5" = "IN",  "Pvalb_Vipr2" = "IN", "Pvalb" = "IN", "Sst" = "IN", "Sst_Chodl_mouse" = "IN",
                              "L2_IT_ENTL" = "EXN", "L2/3_IT_ENTL" = "EXN", "L2/3_IT" = "EXN", "L2/3_IT_Reln" = "EXN", "L3_IT_ENT" = "EXN",
                              "L2/3_IT_PPP_mouse" = "EXN","L2/3_IT_ENTL_bat" = "EXN", "L2/3_IT_ENTL_Maml2_bat" = "EXN", "L2-6" = "EXN", "L2-6_bat" = "EXN",
                              "L2-6_primate" = "EXN", "L4-6" = "EXN", "L5/6_NP" = "EXN", "L5/6_IT_TPE-ENT" = "EXN", "L6b" = "EXN", "L6b/CT_ENT" = "EXN", "L6_IT_ENTL_mouse" = "EXN",
                              "L6b/CT_ENT_mouse" = "EXN", "L6_IT_primate" = "EXN", "CR" = "EXN")
merged.mammal$invsexn <- Idents(merged.mammal)

for (ctp in c("EXN", "IN")){
  subseu <- subset(merged.mammal, invsexn == ctp)
  
  hvg.list <- SplitObject(subseu, split.by = "orig.ident") %>%
    lapply(., function(x) FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)) %>%
    lapply(., function(x) VariableFeatures(x))
  saveRDS(hvg.list, file = paste0(inputdir, "2_HVG.byCTP.", ctp, ".rds"))
}




## Select the top HVG for each major cell type and get the union
sp_list <- list(Human = c("human1", "human2", "human3"),
                Baboon = c("baboon1", "baboon2", "baboon3"), 
                Mouse = c("mouse1", "mouse2", "mouse3"),
                Bat = c("bat3", "bat4", "bat5"))

SelectHVG <- function(hvg.list, nfeatures = 2000) {
  
  var.features <- unname(obj = unlist(x = hvg.list)) %>%
    table() %>%
    sort(., decreasing = TRUE)
  
  tie.val <- var.features[min(nfeatures, length(x = var.features))]
  features <- names(x = var.features[which(x = var.features > tie.val)])
  if (length(x = features) > 0) {
    feature.ranks <- sapply(X = features, FUN = function(x) {
      ranks <- sapply(X = hvg.list, FUN = function(vf) {
        if (x %in% vf) {
          return(which(x = x == vf))
        }
        return(NULL)
      })
      median(x = unlist(x = ranks))
    })
    features <- names(x = sort(x = feature.ranks))
  }
  features.tie <- var.features[which(x = var.features == tie.val)]
  tie.ranks <- sapply(X = names(x = features.tie), FUN = function(x) {
    ranks <- sapply(X = hvg.list, FUN = function(vf) {
      if (x %in% vf) {
        return(which(x = x == vf))
      }
      return(NULL)
    })
    median(x = unlist(x = ranks))
  })
  features <- c(features, names(x = head(x = sort(x = tie.ranks), 
                                         nfeatures - length(x = features))))
  return(features)
}

flist <- list()
for (nfeatures in c(1000, 1250, 1500)){
  final.hvg <- lapply(c("EXN", "IN"), function(ctp) {
    hvg.list <- readRDS(file = paste0(inputdir, "2_HVG.byCTP.", ctp, ".rds"))
    
    
    new.hvg <- lapply(names(sp_list), function(sp) 
      SelectHVG(hvg.list = hvg.list[sp_list[[sp]]], nfeatures = nfeatures)
    ) %>% 
      unlist() %>% unique() %>% unname()
    return(new.hvg)
  }) %>% 
    unlist() %>% unique() %>% unname()
  
  flist[[as.character(nfeatures)]] <- final.hvg
}
saveRDS(flist, file = paste0(inputdir, "2_HVG.byCTP.combine.rds"))



## save <= 100 cells per cluster
merged.mammal@meta.data$spcls <- paste0(as.character(merged.mammal@meta.data$species), "|", as.character(merged.mammal@meta.data$neuronal_celltype))

# check whihch clusters contain cells from all four species
metadata <- merged.mammal@meta.data
length(levels(metadata$neuronal_celltype))
for (i in levels(metadata$neuronal_celltype)) {
  print(i)
  # Focus on one specific cluster
  cluster_of_interest <- as.character(i)
  metadata_i  <- metadata[metadata$neuronal_celltype == cluster_of_interest, ]
  
  # Create a summary table of cell counts per cell types
  distribution <- table(metadata_i$species)
  
  # Convert to a data frame for plotting
  distribution_df <- as.data.frame(distribution)
  colnames(distribution_df) <- c("Cell_type", "Cell_Count")
  
  
  # Normalize counts to percentages for each cluster
  distribution_df <- distribution_df %>%
    mutate(Percentage = Cell_Count / sum(Cell_Count) * 100)
  # Filter out low %
  distribution_df <- distribution_df %>%
    filter(Percentage > 2)
  
  # Add a column with cluster info for plotting
  distribution_df <- distribution_df %>%
    mutate(Cluster = as.character(i))
  
  # Calculate number of cell in that cluster
  total_sum <- sum(distribution_df$Cell_Count)
  
  # Calculate percentage of cell in that cluster after >5% filter
  total_sum_perc <- sum(distribution_df$Percentage)
  
  # Round the sum to one decimal place
  rounded_sum <- round(total_sum_perc, 1)
  
  # Plot the stacked bar plot with percentages displayed
  p <- ggplot(distribution_df, aes(x = Cluster, y = Percentage, fill = Cell_type)) +
    geom_bar(stat = "identity", position = "stack") + # Stacked bar plot
    geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
              position = position_stack(vjust = 0.5), size = 3) + # Add percentage labels
    labs(title = paste("Distribution of Cell Types in Cluster", cluster_of_interest, ", no of Cells:", total_sum, ", Sum of %", rounded_sum),
         x = "Cluster",
         y = "Percentage",
         fill = "Cell Type") +
    theme_minimal()
  
  print(p)
}



clusters <- c( "Vip_primate", "CR","L3_IT_ENT","L2/3_IT_Reln","L2/3_IT_ENTL_Maml2_bat", "L2-6_bat","L2/3_IT_ENTL_bat",
               "L5_PPP_mouse","Sst_Chodl_mouse","L4-6","L2-6","L6_IT_primate","L2/3_IT_ENTL","L6b/CT_ENT","L5/6_IT_TPE-ENT",
               "Sst", "L2-6_primate","Pvalb","L2_IT_ENTL","L2/3_IT","Lamp5","L6b/CT_ENT_mouse","Sncg_Pax6","L6b",
               "L5/6_NP","L6_IT_ENTL_mouse","Pvalb_Vipr2","L2/3_IT_PPP_mouse")

# exclude clusters where not all species are present
clusters <- setdiff(clusters, c("L2/3_IT_PPP_mouse", "L6_IT_ENTL_mouse", "L6b/CT_ENT_mouse", "L2_IT_ENTL","L2-6", "L5_PPP_mouse", "L2-6_bat", "L2/3_IT_ENTL_Maml2_bat", "Sst_Chodl_mouse", "Vip_primate")) 
clusters <- clusters[c(11, 13, 19, 15, 16, 14, 3, 8, 5, 4, 2, 6, 10, 18, 17, 9, 7, 12, 1)]



all_sps <- c("HUMAN", "BABOON", "MOUSE", "BAT")
all_cls <- rep(all_sps, each = length(clusters)) %>%
  paste0(., "|", rep(clusters, times = 4))


set.seed(20210804)
cells <- lapply(all_cls, function(x) {
  subc <- colnames(merged.mammal)[merged.mammal@meta.data$spcls == x]
  if (length(subc) >= 100){
    subc <- sample(subc, 100)
  }
  return(subc)
}) %>%
  unlist()


seu <- merged.mammal[, cells]
Idents(seu) <- "spcls"
avgs <- log(AggregateExpression(seu)$SCT + 1)
saveRDS(avgs, file = paste0(inputdir, "2_Subset.avg.rds"))


avg <- readRDS(file = paste0(inputdir, "2_Subset.avg.rds"))
avg_sup <- readRDS(file = paste0(inputdir, "Subset.avg.sup.rds"))
avg <- cbind(avg, avg_sup[rownames(avg), , drop = FALSE])

hvg <- readRDS(file = paste0(inputdir, "2_HVG.byCTP.combine.rds"))[["1500"]]


## Pair-wise correlation of clusters
sp_pairs <- c("HUMAN|BABOON", "HUMAN|MOUSE", "HUMAN|BAT", "BABOON|MOUSE", "BABOON|BAT", "MOUSE|BAT")

clusters <- gsub("_", "-", clusters)
dense_avg <- as.matrix(avg)

cor_res <- lapply(sp_pairs, function(pair) {
  sps <- strsplit(pair, "|", fixed = TRUE)[[1]]
  sp1 <- sps[1]; sp2 <- sps[2]
  
  cls1 <- paste0(sp1, "|", clusters);
  cls2 <- paste0(sp2, "|", clusters)
  
  aa <- cor(dense_avg[hvg, cls1], dense_avg[hvg, cls2], method = "pearson") 
  res <- diag(aa)%>%
    setNames(., clusters)
  return(res)
}) %>%
  setNames(., sp_pairs) %>%
  as.data.frame(., check.names = FALSE)
div_mat <- (1-cor_res) %>% as.matrix()




pair_order <- rev(c("HUMAN|BABOON", "HUMAN|MOUSE", "BABOON|MOUSE", "HUMAN|BAT", "BABOON|BAT", "MOUSE|BAT"))


col_use <- viridis(10)[c(1,2,5,6,7,8,10)]

gg_heat <- function(data, div_col = "div", cols) {
  p <- ggplot(data, aes_string(x = "cluster", y = "pair", fill = div_col)) +
    geom_tile(width = 1, height = 1, size = 0.1, color = "black") +
    scale_fill_gradientn(colors = cols) +
    theme_classic() +
    RotatedAxis() + 
    coord_fixed() +
    theme(legend.position = "right", axis.text.x = element_text(size = rel(0.5)), axis.line = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.text.y = element_blank())
  return(p)
}


## Get the raw divergence scores
raw_data <- div_mat %>%
  reshape2::melt() %>%
  setNames(., c("cluster", "pair", "div")) %>%
  mutate(cluster = factor(as.character(cluster), levels = clusters)) %>%
  mutate(pair = factor(as.character(pair), levels = pair_order))


raw_data$div <- MinMax(raw_data$div, min = 0, max = quantile(raw_data$div, 0.98))
p_raw <- gg_heat(data = raw_data, div_col = "div", cols = col_use) 




## Get the scaled divergence scores
scale_data <- div_mat %>%
  scale() %>%
  reshape2::melt() %>%
  setNames(., c("cluster", "pair", "div")) %>%
  mutate(cluster = factor(as.character(cluster), levels = clusters)) %>%
  mutate(pair = factor(as.character(pair), levels = pair_order))

p_scale <- gg_heat(data = scale_data, div_col = "div", cols = col_use) 




p_raw <- p_raw
p_scale <- p_scale
pdf(paste0(inputdir, "2_DIV_1minusCor_Subset_combine.pdf"), width = 10, height = 4)
plot_grid(p_raw, p_scale, nrow = 2, ncol = 1, align = "v") %>% print()
dev.off()


