################################################################################

# Inhibitory neuron subtypes

################################################################################


library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)


#### Species distribution IN cell types ####
merged.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammal_annot.rds")
metadata <- merged.mammal@meta.data

# IN species-wise
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


IN_celltypes <- c("Sst_Chodl_mouse", "Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5_Lhx6", "Lamp5", "Sncg_Pax6", "Vip_primate")

filtered_df <- average_proportions %>%
  filter(neuronal_celltype %in% IN_celltypes)
filtered_df$species <- factor(filtered_df$species, levels = c("HUMAN", "BABOON", "MOUSE", "BAT"))
filtered_df$neuronal_celltype <- factor(filtered_df$neuronal_celltype, levels = c("Sst_Chodl_mouse", "Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5_Lhx6", "Lamp5", "Sncg_Pax6", "Vip_primate"))

ggplot(filtered_df, aes(x = neuronal_celltype, y = mean_proportion, fill = species)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Inhibitory Subtype Proportions by Species",
    x = "Subtype",
    y = "Mean Proportion",
    fill = "Species"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "HUMAN" = "#917f2f",
    "BABOON" = "#92bf3e",
    "MOUSE" = "#49348c",
    "BAT" = "#a22d71"
  )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# IN evolution-wise
# Calculate proportions
# Step 1: Calculate cell type proportions per sample
proportions_per_sample <- metadata %>%
  group_by(sample, neuronal_celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(proportion = count / sum(count))


# Step 2: Add evolution info to each sample
sample_evolution <- metadata %>% distinct(sample, evolution)
proportions_per_sample <- proportions_per_sample %>%
  left_join(sample_evolution, by = "sample")

# Step 3: Calculate average proportions per species
average_proportions <- proportions_per_sample %>%
  group_by(evolution, neuronal_celltype) %>%
  summarise(mean_proportion = mean(proportion),
            sd_proportion = sd(proportion),
            .groups = "drop")


IN_celltypes <- c("Sst_Chodl_mouse", "Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5_Lhx6", "Lamp5", "Sncg_Pax6", "Vip_primate")

filtered_df <- average_proportions %>%
  filter(neuronal_celltype %in% IN_celltypes)
filtered_df$evolution <- factor(filtered_df$evolution, levels = c("PRIMATE", "MOUSE", "BAT"))
filtered_df$neuronal_celltype <- factor(filtered_df$neuronal_celltype, levels = c("Sst_Chodl_mouse", "Sst", "Pvalb", "Pvalb_Vipr2", "Lamp5_Lhx6", "Lamp5", "Sncg_Pax6", "Vip_primate"))

ggplot(filtered_df, aes(x = neuronal_celltype, y = mean_proportion, fill = evolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Inhibitory Subtype Proportions by Species",
    x = "Subtype",
    y = "Mean Proportion",
    fill = "Species"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "PRIMATE" = "#d5e59a",
    "MOUSE" = "#49348c",
    "BAT" = "#a22d71"
  )) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        legend.text = element_text(size = 20),
        axis.title.y = element_text(size = 22),
        title = element_text(size = 22),
        plot.margin = margin(l = 20))


#### IN subtypes markers dot plot ####
in.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalIN_annot.rds")



inhibitory.markers <- c("Adarb2", "Nr2f2", "Prox1", "Vip", "Calb2", "Reln", "Kcng1", "Sncg", "Pax6", "Lamp5",
                        "Lhx6", "Sox6", "Slit2", "Pvalb", "Trps1", "Vipr2", "Sst", "Npy", "Chodl")

Idents(in.mammal) <- in.mammal$evo_neuronal_celltype
# Specify the order of your identities
desired_order.in2 <- c("MOUSE_Sst_Chodl_mouse", "PRIMATE_Sst", "MOUSE_Sst", "BAT_Sst",
                       "PRIMATE_Pvalb", "MOUSE_Pvalb", "BAT_Pvalb",  
                       "PRIMATE_Pvalb_Vipr2", "MOUSE_Pvalb_Vipr2", "BAT_Pvalb_Vipr2",  
                       "PRIMATE_Lamp5", "MOUSE_Lamp5", "BAT_Lamp5",
                       "PRIMATE_Sncg_Pax6", "MOUSE_Sncg_Pax6", "BAT_Sncg_Pax6",
                       "PRIMATE_Vip_primate", "MOUSE_Vip_primate", "BAT_Vip_primate")
# Set identities in the specified order
in.mammal@active.ident <- factor(Idents(in.mammal), levels = desired_order.in2)


dotplot <- DotPlot(in.mammal, features = inhibitory.markers)  +  coord_flip() + theme_minimal() + RotatedAxis() 
dotplot + scale_color_gradientn(colors = c("blue", "grey", "red")) + theme(axis.text.y = element_text(face = "italic", size = 20))
