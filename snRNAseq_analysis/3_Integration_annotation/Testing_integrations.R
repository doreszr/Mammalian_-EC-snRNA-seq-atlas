########################################################################################################################

### Different integration methods

########################################################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(scIntegrationMetrics)

#### Integration methods UMAPs ####
# Load primate Seurat object

merged.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/merged.primate_int.rds")

DimPlot(merged.primate, reduction = "umap.harmony", group.by = "sample") + ggtitle("harmony_primate")

DimPlot(merged.primate, reduction = "umap.cca", group.by = "sample") + ggtitle("CCA_primate")

DimPlot(merged.primate, reduction = "umap.rpca", group.by = "sample") + ggtitle("RPCA_primate")

# SCVI
merged.primate.scvi <- readRDS("Data_EC_atlas/Seurat_obj/primate/merged.primate.scvi.rds")

DimPlot(merged.primate, reduction = "umap.scvi", group.by = "sample") + ggtitle("scvi_primate")

#### Benchmarking ####

# calculate metrics
metrics.pca <- getIntegrationMetrics(merged.primate, meta.label = "major_cell_types", meta.batch = "sample", iLISI_perplexity = 20)

metrics.rpca <- getIntegrationMetrics(merged.primate, meta.label = "major_cell_types", meta.batch = "sample", method.reduction = "integrated.rpca", iLISI_perplexity = 20)

metrics.cca <- getIntegrationMetrics(merged.primate, meta.label = "major_cell_types", meta.batch = "sample", method.reduction = "integrated.cca", iLISI_perplexity = 20)

metrics.harmony <- getIntegrationMetrics(merged.primate, meta.label = "major_cell_types", meta.batch = "sample", method.reduction = "harmony", iLISI_perplexity = 20)

metrics.scvi <- getIntegrationMetrics(merged.primate.scvi, meta.label = "cell_types", meta.batch = "sample", method.reduction = "integrated.scvi", iLISI_perplexity = 20)

# concatonate metrics
metrics_to_df <- function(metrics_list, method_name) {
  enframe(metrics_list, name = "metric", value = "value") %>%
    unnest(value) %>%   # handles vectors
    mutate(method = method_name)
}

df_metrics <- bind_rows(
  metrics_to_df(metrics.pca, "PCA"),
  metrics_to_df(metrics.cca, "CCA"),
  metrics_to_df(metrics.rpca, "RPCA"),
  metrics_to_df(metrics.harmony, "Harmony"),
  metrics_to_df(metrics.scvi, "scVI")
)

df_metrics$method <- factor(df_metrics$method, levels = c("PCA", "CCA", "RPCA", "Harmony", "scVI"))

# plot metrics
ggplot(df_metrics, aes(x = method, y = value)) +
  geom_line(group = 1, color = "black") +
  geom_point(aes(color = method), size = 3) +
  facet_wrap(~metric, scales = "free_y") +
  theme_minimal() + theme(
    axis.text.x = element_text(angle = 45, hjust = 1))


df_wide <- df_metrics %>%
  select(method, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

ggplot(df_wide, aes(x = CiLISI, y = iLISI, color = method)) +
  geom_point(size = 4) +
  geom_text(aes(label = method), vjust = -1, size = 4) +
  theme_minimal(base_size = 14) + theme(legend.position = "none")
