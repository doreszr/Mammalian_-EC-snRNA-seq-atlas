################################################################################

# Visualizing cell types in cross-species integrated data set

################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggpubr)

#### Integrated neuronal cell types ####
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")

DimPlot(merged.mammal, reduction = "umap.rpca", label = F, label.size = 2, group.by = "neuronal_celltype", pt.size = 0.05,
        cols = c("Sst" = "#033f63", "Sst_Chodl_mouse" = "#1085CA", "Pvalb" = "#98c1d9", "Pvalb_Vipr2" = "#3d5a80",
                 "Lamp5" = "#b5b682", "Sncg_Pax6" = "#7c9885", "Vip_primate" = "#005A32",
                 "L2/3_IT_Reln" = "#ffcba9", "L2/3_IT_ENTL" = "#ff9b57", "L2_IT_ENTL" = "#ff791e", "L2/3_IT" = "#b34800",
                 "L2/3_IT_ENTL_bat" = "#C23A22", "L2/3_IT_ENTL_Maml2_bat" = "#E54C38", "L6_IT_ENTL_mouse" = "#40130B",
                 "L2/3_IT_PPP_mouse" = "#6B2113", "L3_IT_ENT" = "#972E1A",
                 "L4-6" = "#FBCF3B", "L6_IT_primate" = "#C68E17", "L6b/CT_ENT" = "#fedc97", "L5/6_IT_TPE-ENT" = "#ffb703", "L5/6_NP" = "#FFF880", "L6b" = "#C7A317",
                 "L6b/CT_ENT_mouse" = "#E8A317", "L5_PPP_mouse" = "#d24e01",
                 "L2-6_primate" = "#C4967D", "L2-6_bat" = "#BC746B", "L2-6" = "#fc0c00", "CR" = "#663c20")) + ggtitle("Mammalian neurons") + NoLegend()


#### Number of nuclei in all neuronal cell types ####
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")
metadata <- merged.mammal@meta.data

celltype_count_summary <- metadata %>%
  group_by(neuronal_celltype) %>%
  summarise(count = n()) %>%
  ungroup()

celltype_count_summary$neuronal_celltype <- factor(celltype_count_summary$neuronal_celltype, levels = c("Vip_primate", "Sncg_Pax6", "Lamp5", "Pvalb_Vipr2", "Pvalb", "Sst", "Sst_Chodl_mouse",
                                                                                                        "L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat",
                                                                                                        "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
                                                                                                        "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse",
                                                                                                        "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
                                                                                                        "L2-6_bat", "L2-6_primate", "L2-6", "CR"))



# Plot the data
ggplot(celltype_count_summary, aes(x = neuronal_celltype, y = count)) +
  geom_bar(stat = "identity") + 
  labs(x = "Subtype", y = "Number", title = "No of nuclei") +
  theme_minimal() +
  theme(
    panel.grid = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  coord_flip() +
  theme(panel.grid = element_blank())



#### Species distribution in all neuronal cell types ####
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")
metadata <- merged.mammal@meta.data
# Calculate proportions
# Step 1: Calculate cell type proportions per sample
proportions_per_sample <- metadata %>%
  group_by(sample, neuronal_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))


# Step 2: Add species info to each sample
sample_species <- metadata %>% distinct(sample, species)
proportions_per_sample <- proportions_per_sample %>%
  left_join(sample_species, by = "sample")

# Step 3: Calculate average proportions per species
average_proportions <- proportions_per_sample %>%
  group_by(species, neuronal_celltype) %>%
  summarise(mean_proportion = mean(proportion),
            sd_proportion = sd(proportion),
            .groups = "drop")

# Step 4: Calculate average proportions per species per subtypes
average_proportions_subtypes <- average_proportions %>%
  group_by(neuronal_celltype) %>%
  mutate(normalized_proportion = mean_proportion / sum(mean_proportion)) %>%
  ungroup()

average_proportions_subtypes$species <- factor(average_proportions_subtypes$species, levels = c("HUMAN", "BABOON", "MOUSE", "BAT"))
average_proportions_subtypes$neuronal_celltype <- factor(average_proportions_subtypes$neuronal_celltype, levels = c("Vip_primate", "Sncg_Pax6", "Lamp5", "Pvalb_Vipr2", "Pvalb", "Sst", "Sst_Chodl_mouse",
                                                                                                                    "L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat",
                                                                                                                    "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
                                                                                                                    "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse",
                                                                                                                    "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
                                                                                                                    "L2-6_bat", "L2-6_primate", "L2-6", "CR"))

# Plot the data using ggplot2
ggplot(average_proportions_subtypes, aes(x = neuronal_celltype, y = normalized_proportion, fill = species)) +
  geom_bar(stat = "identity", position = "fill") + 
  labs(x = "Subtype", y = "Percentage", title = "% of nuclei") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Rotate x-axis labels
  scale_y_continuous(labels = scales::percent) +
  coord_flip() + # Display as percentages
  scale_fill_manual(values = c(
    "BABOON" = "#92BE36",
    "HUMAN" = "#927E1B",
    "MOUSE" = "#49338B",
    "BAT" = "#A12E71"
  )) +
  theme(panel.grid = element_blank())

#### Cluster markers dot plot ####
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")


Idents(merged.mammal) <- merged.mammal$neuronal_celltype
desired_order.all <- c("Vip_primate", "Sncg_Pax6", "Lamp5", "Pvalb_Vipr2", "Pvalb", "Sst", "Sst_Chodl_mouse",
                       "L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat",
                       "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
                       "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse",
                       "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
                       "L2-6_bat", "L2-6_primate", "L2-6", "CR")
# Set identities in the specified order
merged.mammal@active.ident <- factor(Idents(merged.mammal), levels = desired_order.all)

interesting.markers <- c("Slc17a7", "Cux2", "Pcdh8", "Lef1", "Chn2", "Plch1", "Fign",  "Reln", "Grik1", "Maml2", "Rorb", "Foxp2", "Tle4", "Bcl11b", "Gad1", "Adarb2", "Lhx6", "Chodl")
dotplot <- DotPlot(merged.mammal, features = interesting.markers) +
  theme_minimal() +
  scale_x_discrete(position = "top") +
  scale_color_gradient2(
    low = "blue",
    mid = "grey",
    high = "red",
    midpoint = 0
  ) +
  theme(
    axis.title.x = element_blank(),
    
    axis.text.x.top = element_text(
      angle = 45,
      hjust = 0,
      vjust = 0,
      face = "italic",
      size = 15,
      margin = margin(b = 2)
    ),
    axis.ticks.x.top = element_line(),
    axis.line.x.top = element_line(),

    axis.text.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank()
  )


