################################################################################

# Sankey plots

################################################################################

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(networkD3)
library(dplyr)
library(htmlwidgets)
library(webshot)

### Cortical layers Sankey plots ####
### EXN subset
exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")

### Create layers metadata
Idents(exn.mammal) <- exn.mammal$neuronal_celltype
superficial.layers <- c("L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat",  "L2/3_IT_ENTL_Maml2_bat",
                        "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT_mouse")


deep.layers <- c("L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP",  "L6b", "L6_IT_ENT_mouse", 
                 "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse")


clusts <- c("L2_IT_ENTL", "L2/3_IT_ENTL", "L2/3_IT_ENTL_bat", "L2/3_IT_ENTL_Maml2_bat", "L2/3_IT_PPP_mouse", "L2/3_IT", "L2/3_IT_Reln", "L3_IT_ENT",
            "L4-6", "L5/6_IT_TPE-ENT",  "L5/6_NP", "L6b", "L6_IT_ENTL_mouse", "L6b/CT_ENT_mouse", "L6b/CT_ENT", "L6_IT_primate", "L5_PPP_mouse",
            "L2-6_bat", "L2-6_primate", "L2-6")

new_idents <- setNames(rep("mix", length(clusts)), clusts) # Default all to "mix"
new_idents[clusts %in% superficial.layers] <- "superficial"
new_idents[clusts %in% deep.layers] <- "deep"

exn.mammal <- RenameIdents(exn.mammal, new_idents)
table(Idents(exn.mammal))

exn.mammal$layers <- Idents(exn.mammal)

### Load EXN species subsets
exn.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_EXN_annot.rds")
exn.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed.primateEXN_annot.rds")
exn.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouseEXN_annot.rds")
# set right ident
Idents(exn.primate) <- exn.primate$neuronal_celltype_species
Idents(exn.bat) <- exn.bat$neuronal_celltype_species
Idents(exn.mouse) <- exn.mouse$neuronal_celltype_species

### Cortical layers PRIMATE-MAMMAL EXN neurons
# Seurat objects
clusters.primate <- data.frame(Cell = colnames(exn.primate), Cluster = Idents(exn.primate))
clusters.mammal <- data.frame(Cell = colnames(exn.mammal), Cluster = Idents(exn.mammal))


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("p_", sankey_data.primate$Cluster_primate),
    paste0("mm_", sankey_data.primate$Cluster_mammal)
    
  ))
)


links.primate <- sankey_data.primate %>%
  mutate(
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    SourceID = match(paste0("p_", Cluster_primate), nodes$name) - 1
    
  ) %>%
  select(SourceID, TargetID, Count)

filtered_links <- links.primate[links.primate$Count >= 50, ]

sankey.raw <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey.raw, "sankey_raw_layers_primate.html", selfcontained = TRUE)

# convert to pdf

webshot("sankey_raw_layers_primate.html", "sankey_raw_layers_primate.pdf", zoom = 3)


# Remove labels by setting them to empty strings
nodes$name <- ""

sankey <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey, "sankey_layers_primate.html", selfcontained = TRUE)

# convert to pdf
webshot("sankey_layers_primate.html", "sankey_layers_primate.pdf")


### Cortical layers BAT-MAMMAL EXN neurons
# Seurat objects
clusters.bat <- data.frame(Cell = colnames(exn.bat), Cluster = Idents(exn.bat))
clusters.mammal <- data.frame(Cell = colnames(exn.mammal), Cluster = Idents(exn.mammal))

# Find cells shared between obj1 and obj3
shared_cells.bat <- intersect(clusters.bat$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
bat_clusters <- clusters.bat%>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_bat = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.bat <- bat_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.bat <- transition_data.bat %>%
  group_by(Cluster_bat, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.bat <- sankey_data.bat %>%
  mutate(
    Cluster_bat_label = paste0("b_", Cluster_bat),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("b_", sankey_data.bat$Cluster_bat),
    paste0("mm_", sankey_data.bat$Cluster_mammal)
  ))
)

links.bat <- sankey_data.bat %>%
  mutate(
    SourceID = match(paste0("b_", Cluster_bat), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

filtered_links <- links.bat[links.bat$Count >= 50, ]

sankey.raw <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey.raw, "sankey_raw_layers_bat.html", selfcontained = TRUE)

# convert to pdf

webshot("sankey_raw_layers_bat.html", "sankey_raw_layers_bat.pdf", zoom = 3)


# Remove labels by setting them to empty strings
nodes$name <- ""

sankey <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey, "sankey_layers_bat.html", selfcontained = TRUE)

# convert to pdf
webshot("sankey_layers_bat.html", "sankey_layers_bat.pdf")


### Cortical layers MOUSE-PRIMATE EXN neurons

# set right ident
Idents(exn.mouse) <- exn.mouse$layers
Idents(exn.primate) <- exn.primate$layers
Idents(exn.mammal) <- exn.mammal$layers

# Seurat objects
clusters.mouse <- data.frame(Cell = colnames(exn.mouse), Cluster = Idents(exn.mouse))
clusters.primate <- data.frame(Cell = colnames(exn.primate), Cluster = Idents(exn.primate))
clusters.mammal <- data.frame(Cell = colnames(exn.mammal), Cluster = Idents(exn.mammal))


# Find cells shared between obj1 and obj3
shared_cells.mouse <- intersect(clusters.mouse$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
mouse_clusters <- clusters.mouse%>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mouse = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.mouse <- mouse_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.mouse <- transition_data.mouse %>%
  group_by(Cluster_mouse, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.mouse <- sankey_data.mouse %>%
  mutate(
    Cluster_mouse_label = paste0("m_", Cluster_mouse),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("m_", sankey_data.mouse$Cluster_mouse),
    paste0("mm_", sankey_data.mouse$Cluster_mammal),
    paste0("p_", sankey_data.primate$Cluster_primate)
  ))
)

links.mouse <- sankey_data.mouse %>%
  mutate(
    SourceID = match(paste0("m_", Cluster_mouse), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

links.primate <- sankey_data.primate %>%
  mutate(
    SourceID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    TargetID = match(paste0("p_", Cluster_primate), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)


links_combined <- bind_rows(links.mouse, links.primate)


# Add a 'group' column to the nodes data frame:
nodes$group <- as.factor(c("superficial", "deep", "CR", "superficial", "deep", "mix", "CR", "superficial", "deep" ))

nodes$label <- as.character(c("","","","","","","","",""))

# Remove labels by setting them to empty strings
nodes$name <- ""

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["superficial","deep", "CR", "mix",]) .range(["#e5e3df", "#231f20", "#ffffff", "#5b5b59"])'

sankey <- sankeyNetwork(Links = links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells",
                        colourScale=my_color, NodeGroup="group")

saveWidget(sankey, "sankey_layers_mouse_primate.html", selfcontained = TRUE)

# convert to pdf
webshot("sankey_layers_mouse_primate.html", "sankey_layers_mouse_primate.pdf")


#### L2-6_primate cluster Sankey plots ####
exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")

# set the right ident
Idents(exn.mammal) <- exn.mammal$neuronal_celltype

# subset right layer
deep.sub.mammal <- subset(exn.mammal, idents = "L2-6")

### Load EXN species subsets
exn.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamedmerged.bat_neu_EXN_annot.rds")
exn.primate <- readRDS("Data_EC_atlas/Seurat_obj/primaterenamed.primateEXN_annot.rds")
exn.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouseEXN_annot.rds")
# set right ident
Idents(exn.primate) <- exn.primate$neuronal_celltype_species
Idents(exn.bat) <- exn.bat$neuronal_celltype_species
Idents(exn.mouse) <- exn.mouse$neuronal_celltype_species

### BAT-PRIMATE

# Seurat objects
clusters.bat <- data.frame(Cell = colnames(exn.bat), Cluster = Idents(exn.bat))
clusters.primate <- data.frame(Cell = colnames(exn.primate), Cluster = Idents(exn.primate))
clusters.mammal <- data.frame(Cell = colnames(deep.sub.mammal), Cluster = Idents(deep.sub.mammal))


# Find cells shared between obj1 and obj3
shared_cells.bat <- intersect(clusters.bat$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
bat_clusters <- clusters.bat%>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_bat = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.bat <- bat_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.bat <- transition_data.bat %>%
  group_by(Cluster_bat, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.bat <- sankey_data.bat %>%
  mutate(
    Cluster_bat_label = paste0("b_", Cluster_bat),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("b_", sankey_data.bat$Cluster_bat),
    paste0("mm_", sankey_data.primate$Cluster_mammal),
    paste0("p_", sankey_data.primate$Cluster_primate)
  ))
)

links.bat <- sankey_data.bat %>%
  mutate(
    SourceID = match(paste0("b_", Cluster_bat), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

links.primate <- sankey_data.primate %>%
  mutate(
    SourceID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    TargetID = match(paste0("p_", Cluster_primate), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)


links_combined <- bind_rows(links.bat, links.primate)

filtered_links_combined <- links_combined[links_combined$Count >= 10, ]


sankey.raw <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey.raw, "plots/sankey_raw_L26_bat_primate110725.html", selfcontained = TRUE)

# convert to pdf

webshot("plots/sankey_raw_L26_bat_primate.html", "plots/sankey_raw_L26_bat_primate.pdf", zoom = 3)

# Remove labels by setting them to empty strings
nodes$name <- ""
sankey <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey, "plots/sankey_L26_bat_primate.html", selfcontained = TRUE)

# convert to pdf

webshot("plots/sankey_L26_bat_primate.html", "plots/sankey_L26_bat_primate.pdf", zoom = 3)



### MOUSE-PRIMATE


# Example for three Seurat objects
clusters.mouse <- data.frame(Cell = colnames(exn.mouse), Cluster = Idents(exn.mouse))
clusters.primate <- data.frame(Cell = colnames(exn.primate), Cluster = Idents(exn.primate))
clusters.mammal <- data.frame(Cell = colnames(deep.sub.mammal), Cluster = Idents(deep.sub.mammal))



# Find cells shared between obj1 and obj3
shared_cells.mouse <- intersect(clusters.mouse$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
mouse_clusters <- clusters.mouse%>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mouse = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.mouse <- mouse_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.mouse <- transition_data.mouse %>%
  group_by(Cluster_mouse, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.mouse <- sankey_data.mouse %>%
  mutate(
    Cluster_mouse_label = paste0("m_", Cluster_mouse),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("m_", sankey_data.mouse$Cluster_mouse),
    paste0("mm_", sankey_data.mouse$Cluster_mammal),
    paste0("p_", sankey_data.primate$Cluster_primate)
  ))
)

links.mouse <- sankey_data.mouse %>%
  mutate(
    SourceID = match(paste0("m_", Cluster_mouse), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

links.primate <- sankey_data.primate %>%
  mutate(
    SourceID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    TargetID = match(paste0("p_", Cluster_primate), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)


links_combined <- bind_rows(links.mouse, links.primate)

filtered_links_combined <- links_combined[links_combined$Count >= 10, ]

sankey.raw <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey.raw, "sankey_raw_L26primate_mouse_primate.html", selfcontained = TRUE)

# convert to pdf
webshot("sankey_raw_L26primate_mouse_primate.html", "sankey_raw_L26primate_mouse_primate.pdf")

# Remove labels by setting them to empty strings
nodes$name <- ""
sankey <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey, "sankey_L26primate_mouse_primate.html", selfcontained = TRUE)

# convert to pdf
webshot("sankey_L26primate_mouse_primate.html", "sankey_L26primate_mouse_primate.pdf")



#### MIxed layer clusters Sankey plots #####

### BAT PRIMATE

exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")
Idents(exn.mammal) <- exn.mammal$neuronal_celltype
deep.sub.mammal <- subset(exn.mammal, idents = c("L2-6_bat", "L2-6", "L2-6_primate"))

exn.bat <- readRDS("Data_EC_atlas/Seurat_obj/fruitbat/renamed.merged.bat_neu_EXN_annot.rds")
Idents(exn.bat) <- exn.bat$neuronal_celltype_species

exn.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed.primateEXN_annot.rds")
Idents(exn.primate) <- exn.primate$neuronal_celltype_species


# Example for three Seurat objects
clusters.bat <- data.frame(Cell = colnames(exn.bat), Cluster = Idents(exn.bat))
clusters.primate <- data.frame(Cell = colnames(exn.primate), Cluster = Idents(exn.primate))
clusters.mammal <- data.frame(Cell = colnames(deep.sub.mammal), Cluster = Idents(deep.sub.mammal))



# Find cells shared between obj1 and obj3
shared_cells.bat <- intersect(clusters.bat$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
bat_clusters <- clusters.bat%>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_bat = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.bat <- bat_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.bat <- transition_data.bat %>%
  group_by(Cluster_bat, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.bat <- sankey_data.bat %>%
  mutate(
    Cluster_bat_label = paste0("b_", Cluster_bat),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("b_", sankey_data.bat$Cluster_bat),
    paste0("mm_", sankey_data.primate$Cluster_mammal),
    paste0("p_", sankey_data.primate$Cluster_primate)
  ))
)

links.bat <- sankey_data.bat %>%
  mutate(
    SourceID = match(paste0("b_", Cluster_bat), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

links.primate <- sankey_data.primate %>%
  mutate(
    SourceID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    TargetID = match(paste0("p_", Cluster_primate), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)


links_combined <- bind_rows(links.bat, links.primate)

filtered_links_combined <- links_combined[links_combined$Count >= 10, ]

sankey.raw <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey.raw, "sankey_raw_mixed_bat_primate.html", selfcontained = TRUE)

# convert to pdf

webshot("sankey_raw_mixed_bat_primate.html", "sankey_raw_mixed_bat_primate.pdf", zoom = 3)


# Remove labels by setting them to empty strings
nodes$name <- ""
sankey <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey, "sankey_mixed_bat_primate.html", selfcontained = TRUE)

# convert to pdf

webshot("/sankey_mixed_bat_primate.html", "sankey_mixed_bat_primate.pdf", zoom = 3)




### MOUSE PRIMATE

exn.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalEXN_annot.rds")
Idents(exn.mammal) <- exn.mammal$neuronal_celltype
deep.sub.mammal <- subset(exn.mammal, idents = c("L2-6_bat", "L2-6", "L2-6_primate"))

exn.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouseEXN_annot.rds")
Idents(exn.mouse) <- exn.mouse$neuronal_celltype_species

exn.primate <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed.primateEXN_annot.rds")
Idents(exn.primate) <- exn.primate$neuronal_celltype_species



# Example for three Seurat objects
clusters.mouse <- data.frame(Cell = colnames(exn.mouse), Cluster = Idents(exn.mouse))
clusters.primate <- data.frame(Cell = colnames(exn.primate), Cluster = Idents(exn.primate))
clusters.mammal <- data.frame(Cell = colnames(deep.sub.mammal), Cluster = Idents(deep.sub.mammal))



# Find cells shared between obj1 and obj3
shared_cells.mouse <- intersect(clusters.mouse$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
mouse_clusters <- clusters.mouse%>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mouse = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.mouse <- mouse_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.mouse <- transition_data.mouse %>%
  group_by(Cluster_mouse, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.mouse <- sankey_data.mouse %>%
  mutate(
    Cluster_mouse_label = paste0("m_", Cluster_mouse),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("m_", sankey_data.mouse$Cluster_mouse),
    paste0("mm_", sankey_data.mouse$Cluster_mammal),
    paste0("p_", sankey_data.primate$Cluster_primate)
  ))
)

links.mouse <- sankey_data.mouse %>%
  mutate(
    SourceID = match(paste0("m_", Cluster_mouse), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

links.primate <- sankey_data.primate %>%
  mutate(
    SourceID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    TargetID = match(paste0("p_", Cluster_primate), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)


links_combined <- bind_rows(links.mouse, links.primate)

filtered_links_combined <- links_combined[links_combined$Count >= 10, ]

sankey.raw <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey.raw, "sankey_raw_mixed_mouse_primate.html", selfcontained = TRUE)

# convert to pdf

webshot("sankey_raw_mixed_mouse_primate.html", "sankey_raw_mixed_mouse_primate.pdf", zoom = 3)


# Remove labels by setting them to empty strings
nodes$name <- ""
sankey <- sankeyNetwork(Links = filtered_links_combined, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells")

saveWidget(sankey, "sankey_mixed_mouse_primate.html", selfcontained = TRUE)

# convert to pdf

webshot("sankey_mixed_mouse_primate.html", "sankey_mixed_mouse_primate.pdf", zoom = 3)






#### IN Sankey plots ####

# load IN subsets
in.mammal <- readRDS("Data_EC_atlas/Seurat_obj/mammals/mammalIN_annot.rds")
in.primates <- readRDS("Data_EC_atlas/Seurat_obj/primate/renamed.primateIN_annot.rds")
in.mouse <- readRDS("Data_EC_atlas/Seurat_obj/mouse/merged.mouseIN_annot.rds")
in.bat <- readRDS("Data_EC_atlas/Seurat_obj/renamed.merged.bat_neu_IN_annot.rds")


### MOUSE-PRIMATE

# set right ident
Idents(in.mouse) <- in.mouse$neuronal_celltype_species
Idents(in.primates) <- in.primates$neuronal_celltype_species
Idents(in.mammal) <- in.mammal$neuronal_celltype

# Seurat objects
clusters.mouse <- data.frame(Cell = colnames(in.mouse), Cluster = Idents(in.mouse))
clusters.primate <- data.frame(Cell = colnames(in.primates), Cluster = Idents(in.primates))
clusters.mammal <- data.frame(Cell = colnames(in.mammal), Cluster = Idents(in.mammal))


# Find cells shared between obj1 and obj3
shared_cells.mouse <- intersect(clusters.mouse$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
mouse_clusters <- clusters.mouse%>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mouse = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.mouse) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.mouse <- mouse_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.mouse <- transition_data.mouse %>%
  group_by(Cluster_mouse, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.mouse <- sankey_data.mouse %>%
  mutate(
    Cluster_mouse_label = paste0("m_", Cluster_mouse),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("m_", sankey_data.mouse$Cluster_mouse),
    paste0("mm_", sankey_data.mouse$Cluster_mammal),
    paste0("p_", sankey_data.primate$Cluster_primate)
  ))
)

links.mouse <- sankey_data.mouse %>%
  mutate(
    SourceID = match(paste0("m_", Cluster_mouse), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

links.primate <- sankey_data.primate %>%
  mutate(
    SourceID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    TargetID = match(paste0("p_", Cluster_primate), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)


links_combined <- bind_rows(links.mouse, links.primate)

filtered_links <- links_combined[links_combined$Count >= 10, ]

sankey.raw <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")



# Add a 'group' column to the nodes data frame:
nodes$group <- as.factor(c("Lamp5", "Sst_Chodl", "Lamp5_Lhx6", "Pvalb_Vipr2", "Vip", "Sncg","Pvalb", "Sst", "Sst", "Pvalb", "Lamp5", "Sncg_Pax6", "Sst_Chodl_mouse", "Pvalb_Vipr2", "Vip", "Vip_primate", "Sst", "Pvalb", "Lamp5", "Pax6_Adarb2", "Lamp5_Lhx6", "Pvalb_Vipr2", "Sst_Chodl"))

nodes$label <- as.character(rep("", length(nodes$group)))

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["Lamp5", "Sst_Chodl", "Lamp5_Lhx6", "Pvalb_Vipr2", "Vip", "Vip_primate", "Sncg", "Pvalb", "Sst", "Sncg_Pax6", "Sst_Chodl_mouse", "Pax6_Adarb2"]) .range(["#b5b682", "#1085CA", "#28666e", "#3d5a80", "#005A32", "#005A32", "#7c9885", "#98c1d9", "#033f63", "#7c9885", "#1085CA", "#7c9885"])' 

nodes$color <- ifelse(nodes$group == "Lamp5", "#b5b682",
                      ifelse(nodes$group == "Pvalb_Vipr2", "#3d5a80",
                             ifelse(nodes$group == "Pvalb", "#98c1d9",
                                    ifelse(nodes$group %in% c("Sst_Chodl", "Sst_Chodl_mouse"), "#1085CA",
                                           ifelse(nodes$group == "Lamp5_Lhx6", "#28666e",
                                                  ifelse(nodes$group == "Sst", "#033f63",
                                                         ifelse(nodes$group %in% c("Sncg", "Sncg_Pax6", "Pax6_Adarb2"), "#7c9885", "#fedc97",
                                                                ifelse(nodes$group %in% c("Vip", "Vip_primate"), "#005A32"))))))))

# Remove labels by setting them to empty strings
nodes$name <- ""

# Use color column instead of colourScale
sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
              Value = "Count", NodeID = "name", units = "cells",
              NodeGroup = "color", colourScale = my_color)

sankey <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells",
                        NodeGroup = "color", colourScale = my_color)

saveWidget(sankey, "sankey_IN_mouse_primate.html", selfcontained = TRUE)

# convert to pdf
webshot("sankey_IN_mouse_primate.html", "sankey_IN_mosue_primate.pdf")


### bat-primate

# set right ident
Idents(in.bat) <- in.bat$neuronal_celltype_species
Idents(in.primates) <- in.primates$neuronal_celltype_species
Idents(in.mammal) <- in.mammal$neuronal_celltype

# Example for three Seurat objects
clusters.bat <- data.frame(Cell = colnames(in.bat), Cluster = Idents(in.bat))
clusters.primate <- data.frame(Cell = colnames(in.primates), Cluster = Idents(in.primates))
clusters.mammal <- data.frame(Cell = colnames(in.mammal), Cluster = Idents(in.mammal))



# Find cells shared between obj1 and obj3
shared_cells.bat <- intersect(clusters.bat$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
bat_clusters <- clusters.bat%>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_bat = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.bat) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.bat <- bat_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.bat <- transition_data.bat %>%
  group_by(Cluster_bat, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.bat <- sankey_data.bat %>%
  mutate(
    Cluster_bat_label = paste0("b_", Cluster_bat),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


# Find cells shared between obj1 and obj3
shared_cells.primate <- intersect(clusters.primate$Cell, clusters.mammal$Cell)

# Subset cluster information for shared cells
primate_clusters <- clusters.primate%>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_primate = Cluster)

mammal_clusters <- clusters.mammal %>%
  filter(Cell %in% shared_cells.primate) %>%
  select(Cell, Cluster) %>%
  rename(Cluster_mammal = Cluster)

# Merge the cluster data for shared cells
transition_data.primate <- primate_clusters %>%
  inner_join(mammal_clusters, by = "Cell")

# Count the number of transitions between clusters
sankey_data.primate <- transition_data.primate %>%
  group_by(Cluster_primate, Cluster_mammal) %>%
  summarise(Count = n()) %>%
  ungroup()

sankey_data.primate <- sankey_data.primate %>%
  mutate(
    Cluster_primate_label = paste0("p_", Cluster_primate),
    Cluster_mammal_label = paste0("mm_", Cluster_mammal)
  )


nodes <- data.frame(
  name = unique(c(
    paste0("b_", sankey_data.bat$Cluster_bat),
    paste0("mm_", sankey_data.bat$Cluster_mammal),
    paste0("p_", sankey_data.primate$Cluster_primate)
  ))
)

links.bat <- sankey_data.bat %>%
  mutate(
    SourceID = match(paste0("b_", Cluster_bat), nodes$name) - 1,
    TargetID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)

links.primate <- sankey_data.primate %>%
  mutate(
    SourceID = match(paste0("mm_", Cluster_mammal), nodes$name) - 1,
    TargetID = match(paste0("p_", Cluster_primate), nodes$name) - 1
  ) %>%
  select(SourceID, TargetID, Count)


links_combined <- bind_rows(links.bat, links.primate)

filtered_links <- links_combined[links_combined$Count >= 10, ]

sankey.raw <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                            Value = "Count", NodeID = "name", units = "cells")

# Add a 'group' column to the nodes data frame:
nodes$group <- as.factor(c("Sst", "Sncg/Vip", "Lamp5", "Lamp5_Lhx6", "Pvalb_Vipr2", "Pvalb", "Sst", "Pvalb", "Vip_primate", "Lamp5", "Sncg_Pax6", "Pvalb_Vipr2", "Vip", "Sst", "Pvalb", "Lamp5", "Pax6_Adarb2", "Lamp5_Lhx6", "Pvalb_Vipr2", "Sst_Chodl"))

nodes$label <- as.character(rep("", length(nodes$group)))

# Give a color for each group:
my_color <- 'd3.scaleOrdinal() .domain(["Sst", "Sncg/Vip", "Lamp5", "Lamp5_Lhx6", "Pvalb_Vipr2", "Pvalb", "Vip", "Vip_primate", "Sncg_Pax6", "Pax6_Adarb2", "Sst_Chodl"]) .range(["#033f63", "#7c9885", "#b5b682", "#28666e", "#3d5a80", "#98c1d9", "#005A32", "#005A32", "#7c9885","#7c9885", "#1085CA"])'

# Remove labels by setting them to empty strings
nodes$name <- ""

sankey <- sankeyNetwork(Links = filtered_links, Nodes = nodes, Source = "SourceID", Target = "TargetID", 
                        Value = "Count", NodeID = "name", units = "cells",
                        colourScale=my_color, NodeGroup="group")

saveWidget(sankey, "sankey_IN_bat_primate.html", selfcontained = TRUE)

# convert to pdf
webshot("sankey_IN_bat_primate.html", "sankey_IN_bat_primate.pdf")

