################################################################################

# Module scores

################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)

#### Loading genes ####

### Learning GO terms

# GO:0007611 learning and memory
learning.memory <- read_tsv("GO_0007611.tsv")

learning.memory <- as.data.frame(learning.memory)

learning.memory <- learning.memory %>% distinct(SYMBOL)
learning.memory.genes <- learning.memory$SYMBOL

## GO:0008542 visual learning
visual.learning <- read_tsv("GO_0008542.tsv")

visual.learning <- as.data.frame(visual.learning)

visual.learning <- visual.learning %>% distinct(SYMBOL)
visual.learning.genes <- visual.learning$SYMBOL

## GO:0008355 olfactory learning
olfactory <- c("Gucy2d", "Grin1", "Mup20", "Adcy3", "Tpbg")

## GO:0042297 vocal learning
vocal <- c("Stra6", "Nrxn2", "Cntnap2", "Htt", "Foxp2", "Kiaa0319")

## GO:0031223 auditory behavior
auditory.behavior <- c("Stra6", "Nrxn2", "Cntnap2", "Htt", "Slc1a3", "Foxp2", "Drd2", "Neurog1", "Abl2", "Kiaa0319", "Slitrk6", "Tifab")


list_GO<- list(learning.memory.genes, visual.learning.genes, olfactory, auditory.behavior)
names(list_GO) <- c("memory", "visual", "olfactory", "auditory")




### Disease associated genes

# WP_ALZHEIMERS_DISEASE mouse gene set
AD.genes <- c("Adam10","Adam17","Apaf1","Apbb1","Apoe","App","Atp2a1","Atp2a2","Bad","Bid","Cacna1c","Cacna1d","Cacna1s","Calm2","Capn1","Capn2","Casp12","Casp3",
              "Casp7","Casp8","Casp9","Cdk5","Cdk5r1","Eif2ak3","Fas","Gnaq","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Hsd17b10","Ide","Il1b","Itpr1","Itpr2","Itpr3",
              "Lpl","Lrp1","Mme","Mapt","Nos1","Plcb1","Plcb2","Plcb3","Plcb4","Ppp3ca","Ppp3cb","Ppp3cc","Ppp3r1","Ppp3r2","Psen1","Psen2","Ryr3","Snca","Tnf","Tnfrsf1a",
              "Trp53","Aph1a","Atf6","Nae1","Bace1","Mapk1","Mapk3","Atp2a3","Cacna1f","Chp1","Gsk3b","Ncstn","Psenen","Chp2","Calml3","Ern1","Calm4")
AD.list <- list(AD.genes)

# WP_PARKINSONS_DISEASE mouse gene set
Parkinson.genes <- c("Mir873a","Mir1224","Mir873b","Mir30f","Eprs1","Apaf1","Casp2","Casp3","Casp6","Casp7","Casp9","Ccne1","Ccne2","Cycs","Slc6a3","Ddc","Ube2j2","Gpr37",
                     "Septin5","Mapk11","Atxn2","Snca","Th","Ubb","Ube2l3","Uba1","Ube2g2","Uchl1","Syt11","Mapk13","Mapk14","Mapk12","Mir127","Mir136","Mir212","Mir26b","Mirlet7g",
                     "Prkn","Ube2j1","Ube2l6","Park7","Htra2","Lrrk2","Ube2g1","Sncaip","Pink1","Mir30e","Mir338","Mir34b","Mir370","Mir409","Mir485","Mir503","Mir19a","Mir10a","Mir375",
                     "Mir433","Mir16-2","Uba7")


list_disease<- list(AD.genes, Parkinson.genes)
names(list_disease) <- c("AD", "PD")

#### INTEGRATED ####
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

merged.mammal <- AddModuleScore(
  object = merged.mammal,
  features = list_GO,
  name = "learning"
)
colnames(merged.mammal@meta.data)[grep("learning", colnames(merged.mammal@meta.data))] <- names(list_GO)


merged.mammal <- AddModuleScore(
  object = merged.mammal,
  features = list_disease,
  name = "disease"
)
colnames(merged.mammal@meta.data)[grep("disease", colnames(merged.mammal@meta.data))] <- names(list_disease)

# UMAPs
FeaturePlot(merged.mammal, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_GO), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()



FeaturePlot(merged.mammal, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_disease), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()

# violoin plots
merged.mammal$neuronal_celltype <- factor(merged.mammal$neuronal_celltype,
                                          levels = c("Vip_primate", "Sncg_Pax6", "Lamp5", "Pvalb_Vipr2", "Pvalb", "Sst", "Sst_Chodl_mouse",
                                                     "L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat",
                                                     "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
                                                     "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse",
                                                     "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
                                                     "L2-6_bat", "L2-6_primate", "L2-6", "CR"))

features <- c("memory", "visual", "olfactory", "auditory", "AD", "PD")

plots <- lapply(features, function(f) {
  VlnPlot(
    merged.mammal,
    features = f,
    group.by = "neuronal_celltype",
    pt.size = 0
  ) + NoLegend() +
    ggtitle(f) + geom_hline(yintercept = 0, color = "black") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

for (i in seq_along(plots)) {
  ggsave(
    filename = paste0("plots/modulescores_violin_", features[i], ".png"),
    plot = plots[[i]],
    dpi = 300,
    width = 8,
    height = 5
  )
}

### Statistical testing
metadata <- merged.mammal@meta.data
celltypes <- unique(merged.mammal$neuronal_celltype)

all_results <- list()

for (i in features) {
  
  df_sample <- metadata %>%
    group_by(sample, neuronal_celltype) %>%
    summarise(mean_score = mean(.data[[i]]), .groups = "drop")
  
  results <- data.frame()
  
  for (ct in celltypes) {
    
    in_group <- df_sample$mean_score[df_sample$neuronal_celltype == ct]
    out_group <- df_sample$mean_score[df_sample$neuronal_celltype != ct]
    
    if (length(unique(df_sample$sample[df_sample$neuronal_celltype == ct])) < 2) next
    
    test <- wilcox.test(in_group, out_group)
    
    median_in <- median(in_group)
    median_out <- median(out_group)
    
    results <- rbind(results, data.frame(
      feature = i,
      celltype = ct,
      p_value = test$p.value,
      median_in = median_in,
      median_out = median_out,
      effect_size = median_in - median_out,
      n_in_samples = length(in_group),
      n_out_samples = length(out_group)
    ))
  }
  
  results$adj_p <- p.adjust(results$p_value, method = "BH")
  
  all_results[[i]] <- results
}

final_results <- do.call(rbind, all_results)

# show effect size on heatmap
heat_df <- final_results %>%
  select(feature, celltype, effect_size) %>%
  pivot_wider(
    names_from = celltype,
    values_from = effect_size
  )

heat_long <- heat_df %>%
  pivot_longer(
    cols = -feature,
    names_to = "celltype",
    values_to = "effect_size"
  )
heat_mat <- heat_df %>%
  column_to_rownames("feature") %>%
  as.matrix()

heat_long <- heat_long %>%
  left_join(final_results %>% select(feature, celltype, adj_p),
            by = c("feature", "celltype"))

cell_order <- c("Vip_primate", "Sncg_Pax6", "Lamp5", "Pvalb_Vipr2", "Pvalb", "Sst", "Sst_Chodl_mouse",
                "L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat",
                "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
                "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse",
                "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
                "L2-6_bat", "L2-6_primate", "L2-6", "CR")

feature_order <- c("memory", "visual", "auditory", "olfactory", "AD", "PD")


heat_long$celltype <- factor(heat_long$celltype, levels = cell_order)
heat_long$feature <- factor(heat_long$feature, levels = feature_order)

ggplot(heat_long, aes(x = celltype, y = feature, fill = effect_size)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) + geom_text(aes(label = ifelse(adj_p < 0.1, "*", ""))) + ggtitle("Module Score Enrichment Across Cell Types")



#### PRIMATE ####
merged.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed_primate_annot.rds")

merged.primate <- AddModuleScore(
  object = merged.primate,
  features = list_GO,
  name = "learning"
)
colnames(merged.primate@meta.data)[grep("learning", colnames(merged.primate@meta.data))] <- names(list_GO)


merged.primate <- AddModuleScore(
  object = merged.primate,
  features = list_disease,
  name = "disease"
)
colnames(merged.primate@meta.data)[grep("disease", colnames(merged.primate@meta.data))] <- names(list_disease)


# UMAPs
FeaturePlot(merged.primate, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_GO), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()

FeaturePlot(merged.primate, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_disease), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()


# violoin plots
plots <- lapply(features, function(f) {
  VlnPlot(
    merged.primate,
    features = f,
    group.by = "neuronal_celltype_species",
    pt.size = 0
  ) + NoLegend() +
    ggtitle(f) + geom_hline(yintercept = 0, color = "black") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

for (i in seq_along(plots)) {
  ggsave(
    filename = paste0("plots/modulescores_PRIMATEviolin_", features[i], ".png"),
    plot = plots[[i]],
    dpi = 300,
    width = 8,
    height = 5
  )
}

### Statistical testing
metadata <- merged.primate@meta.data
celltypes <- unique(merged.primate$neuronal_celltype_species)

all_results <- list()

for (i in features) {
  
  df_sample <- metadata %>%
    group_by(sample, neuronal_celltype_species) %>%
    summarise(mean_score = mean(.data[[i]]), .groups = "drop")
  
  results <- data.frame()
  
  for (ct in celltypes) {
    
    in_group <- df_sample$mean_score[df_sample$neuronal_celltype_species == ct]
    out_group <- df_sample$mean_score[df_sample$neuronal_celltype_species != ct]
    
    if (length(unique(df_sample$sample[df_sample$neuronal_celltype_species == ct])) < 2) next
    
    test <- wilcox.test(in_group, out_group)
    
    median_in <- median(in_group)
    median_out <- median(out_group)
    
    results <- rbind(results, data.frame(
      feature = i,
      celltype = ct,
      p_value = test$p.value,
      median_in = median_in,
      median_out = median_out,
      effect_size = median_in - median_out,
      n_in_samples = length(in_group),
      n_out_samples = length(out_group)
    ))
  }
  
  results$adj_p <- p.adjust(results$p_value, method = "BH")
  
  all_results[[i]] <- results
}

final_results <- do.call(rbind, all_results)

# show effect size on heatmap
heat_df <- final_results %>%
  select(feature, celltype, effect_size) %>%
  pivot_wider(
    names_from = celltype,
    values_from = effect_size
  )

heat_long <- heat_df %>%
  pivot_longer(
    cols = -feature,
    names_to = "celltype",
    values_to = "effect_size"
  )
heat_mat <- heat_df %>%
  column_to_rownames("feature") %>%
  as.matrix()

heat_long <- heat_long %>%
  left_join(final_results %>% select(feature, celltype, adj_p),
            by = c("feature", "celltype"))


feature_order <- c("memory", "visual", "auditory", "olfactory", "AD", "PD")


heat_long$feature <- factor(heat_long$feature, levels = feature_order)

ggplot(heat_long, aes(x = celltype, y = feature, fill = effect_size)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) + geom_text(aes(label = ifelse(adj_p < 0.1, "*", "")))


### BAT ####
merged.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_annot.rds")

merged.bat <- AddModuleScore(
  object = merged.bat,
  features = list_GO,
  name = "learning"
)
colnames(merged.bat@meta.data)[grep("learning", colnames(merged.bat@meta.data))] <- names(list_GO)



merged.bat <- AddModuleScore(
  object = merged.bat,
  features = list_disease,
  name = "disease"
)
colnames(merged.bat@meta.data)[grep("disease", colnames(merged.bat@meta.data))] <- names(list_disease)

# UMAPs
FeaturePlot(merged.bat, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_GO), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()

FeaturePlot(merged.bat, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_disease), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()


# violoin plots
plots <- lapply(features, function(f) {
  VlnPlot(
    merged.bat,
    features = f,
    group.by = "neuronal_celltype_species",
    pt.size = 0
  ) + NoLegend() +
    ggtitle(f) + geom_hline(yintercept = 0, color = "black") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

for (i in seq_along(plots)) {
  ggsave(
    filename = paste0("plots/modulescores_BATviolin_", features[i], ".png"),
    plot = plots[[i]],
    dpi = 300,
    width = 8,
    height = 5
  )
}

### Statistical testing
metadata <- merged.bat@meta.data
celltypes <- unique(merged.bat$neuronal_celltype_species)

all_results <- list()

for (i in features) {
  
  df_sample <- metadata %>%
    group_by(sample, neuronal_celltype_species) %>%
    summarise(mean_score = mean(.data[[i]]), .groups = "drop")
  
  results <- data.frame()
  
  for (ct in celltypes) {
    
    in_group <- df_sample$mean_score[df_sample$neuronal_celltype_species == ct]
    out_group <- df_sample$mean_score[df_sample$neuronal_celltype_species != ct]
    
    if (length(unique(df_sample$sample[df_sample$neuronal_celltype_species == ct])) < 2) next
    
    test <- wilcox.test(in_group, out_group)
    
    median_in <- median(in_group)
    median_out <- median(out_group)
    
    results <- rbind(results, data.frame(
      feature = i,
      celltype = ct,
      p_value = test$p.value,
      median_in = median_in,
      median_out = median_out,
      effect_size = median_in - median_out,
      n_in_samples = length(in_group),
      n_out_samples = length(out_group)
    ))
  }
  
  results$adj_p <- p.adjust(results$p_value, method = "BH")
  
  all_results[[i]] <- results
}

final_results <- do.call(rbind, all_results)

# show effect size on heatmap
heat_df <- final_results %>%
  select(feature, celltype, effect_size) %>%
  pivot_wider(
    names_from = celltype,
    values_from = effect_size
  )

heat_long <- heat_df %>%
  pivot_longer(
    cols = -feature,
    names_to = "celltype",
    values_to = "effect_size"
  )
heat_mat <- heat_df %>%
  column_to_rownames("feature") %>%
  as.matrix()

heat_long <- heat_long %>%
  left_join(final_results %>% select(feature, celltype, adj_p),
            by = c("feature", "celltype"))


feature_order <- c("memory", "visual", "auditory", "olfactory", "AD", "PD")


heat_long$feature <- factor(heat_long$feature, levels = feature_order)

ggplot(heat_long, aes(x = celltype, y = feature, fill = effect_size)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) + geom_text(aes(label = ifelse(adj_p < 0.1, "*", "")))



#### MOUSE ####
merged.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouse_annot.rds")

merged.mouse <- AddModuleScore(
  object = merged.mouse,
  features = list_GO,
  name = "learning"
)
colnames(merged.mouse@meta.data)[grep("learning", colnames(merged.mouse@meta.data))] <- names(list_GO)


merged.mouse <- AddModuleScore(
  object = merged.mouse,
  features = list_disease,
  name = "disease"
)
colnames(merged.mouse@meta.data)[grep("disease", colnames(merged.mouse@meta.data))] <- names(list_disease)

# UMAPs
FeaturePlot(merged.mouse, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_GO), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()


FeaturePlot(merged.mouse, reduction = "umap.rpca", pt.size = 0.1,
            features = names(list_disease), order = TRUE) &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()


# violoin plots
plots <- lapply(features, function(f) {
  VlnPlot(
    merged.mouse,
    features = f,
    group.by = "neuronal_celltype_species",
    pt.size = 0
  ) + NoLegend() +
    ggtitle(f) + geom_hline(yintercept = 0, color = "black") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
})

for (i in seq_along(plots)) {
  ggsave(
    filename = paste0("plots/modulescores_MOUSEviolin_", features[i], ".png"),
    plot = plots[[i]],
    dpi = 300,
    width = 8,
    height = 5
  )
}

### Statistical testing
metadata <- merged.mouse@meta.data
celltypes <- unique(merged.mouse$neuronal_celltype_species)

all_results <- list()

for (i in features) {
  
  df_sample <- metadata %>%
    group_by(sample, neuronal_celltype_species) %>%
    summarise(mean_score = mean(.data[[i]]), .groups = "drop")
  
  results <- data.frame()
  
  for (ct in celltypes) {
    
    in_group <- df_sample$mean_score[df_sample$neuronal_celltype_species == ct]
    out_group <- df_sample$mean_score[df_sample$neuronal_celltype_species != ct]
    
    if (length(unique(df_sample$sample[df_sample$neuronal_celltype_species == ct])) < 2) next
    
    test <- wilcox.test(in_group, out_group)
    
    median_in <- median(in_group)
    median_out <- median(out_group)
    
    results <- rbind(results, data.frame(
      feature = i,
      celltype = ct,
      p_value = test$p.value,
      median_in = median_in,
      median_out = median_out,
      effect_size = median_in - median_out,
      n_in_samples = length(in_group),
      n_out_samples = length(out_group)
    ))
  }
  
  results$adj_p <- p.adjust(results$p_value, method = "BH")
  
  all_results[[i]] <- results
}

final_results <- do.call(rbind, all_results)

# show effect size on heatmap
heat_df <- final_results %>%
  select(feature, celltype, effect_size) %>%
  pivot_wider(
    names_from = celltype,
    values_from = effect_size
  )

heat_long <- heat_df %>%
  pivot_longer(
    cols = -feature,
    names_to = "celltype",
    values_to = "effect_size"
  )
heat_mat <- heat_df %>%
  column_to_rownames("feature") %>%
  as.matrix()

heat_long <- heat_long %>%
  left_join(final_results %>% select(feature, celltype, adj_p),
            by = c("feature", "celltype"))


feature_order <- c("memory", "visual", "auditory", "olfactory", "AD", "PD")


heat_long$feature <- factor(heat_long$feature, levels = feature_order)

ggplot(heat_long, aes(x = celltype, y = feature, fill = effect_size)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) + geom_text(aes(label = ifelse(adj_p < 0.1, "*", "")))


