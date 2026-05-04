# Generate reference data set for cell typing using Allen mouse snRNA-seq atlas

## Data
Data was downloaded from [Allen Brain Map](https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x)

Information about the data from the web site:

"This data set includes single-cell transcriptomes from multiple cortical areas and the hippocampal formation, including 1.1M total cells. Samples were collected from dissections of brain regions from ~8 week-old male and female mice, from pan-neuronal transgenic lines."

## Creating Seurat object from EC
The data set was immediately subsetted to contain only cells from EC and was then read into Seurat.

A neuron only EC reference data set was also generated using the same input.

Input: ```/home/projects/misc/evoec/production/mouseAllen/input/expression_matrix.hdf5```
```/home/projects/misc/evoec/production/mouseAllen/input/metadata_mouse.csv```

Output:
EC cells ```/home/projects/misc/evoec/production/mouseAllen/output/mouse_atlasv5.rds```

EC neurons ```/home/projects/misc/evoec/production/mouseAllen/output/mouse_atlasv5_neu.rds```

```lodadmouseallen.R``` and ```loadmouseallen_neu.R``` were run like:

```bash
#!/bin/bash

  

source /home/local/opt/miniconda3/etc/profile.d/conda.sh

  

Rscript loadmouseallen.R
```

