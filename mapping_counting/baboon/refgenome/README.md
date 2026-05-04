# Preparing reference genome for mapping baboon snRNA-seq data

## Filtering baboon gtf for protein-coding genes
Because CellRanger *mkgtf* command cannot handle RefSeq reference gtf files

Reference genome used as input:
GCF_008728515.1


### Get gene IDs of protein-coding genes from gtf file
Get the gene IDs of genes that are protein-coding from a gtf file from RefSeq

Input: ```GCF_008728515.1_Panubis1.0_genomic.gtf```

Output:
```baboon/refgenome/output/geneids.txt```

Run it like: 

```bash
python get_gene_ids.py genomic.gtf idlist.txt
```

### Filter gtf based on ID list created in previous step
Filter gtf file from RefSeq so output gtf only contains protein-coding genes. Using gene ID list of protein-coding genes from previous step to filter. Special case for mitochondrial genes - fixing those lines as well.

Input: ```GCF_008728515.1_Panubis1.0_genomic.gtf```
```baboon/refgenome/output/geneids.txt```

Output:
```baboon/refgenome/output/baboon.gtf```


Run it like: 

```bash
python filtering_for_cellranger.py genomic.gtf idlist output.gtf
```
