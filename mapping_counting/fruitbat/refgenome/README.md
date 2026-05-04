# Preparing reference genome for mapping fruit bat snRNA-seq data

Reference genome used as input:
GCF_000151845.1


## Preparing reference fasta for mapping
This reference genome only contains genomic sequence so the mitochondrial part had to be added according to the following steps:

1. Downloaded mitochondrial fasta from https://www.ncbi.nlm.nih.gov/nuccore/NC_026542.1?report=fasta into a text file without the header
2. transform it to single line (```cat NC_026542.1.fna |  tr -d '\n' > single_line_mito.txt```)
3. Run the single_to_multifasta.py on that single_line_mito.txt so line length matches genomic fasta line length
4. Add header from https://www.ncbi.nlm.nih.gov/nuccore/NC_026542.1?report=fasta (>NC_026542.1 Pteropus vampyrus mitochondrion, complete genome)
5. Concatonate to genomic fasta file via cat (```bash cat *fna > combined.fasta```) when in a folder genomic.fna and mito.fna so mito.fna is added at the end of genomic

Input: ```GCF_000151845.1_Pvam_2.0_genomic.fna```

```NC_026542.1.fna```

Output:
from step 2 ```fruitbat/refgenome/fasta/output/single_line_mito_bat.txt```

from step 3 ```fruitbat/refgenome/fasta/output/fixednoheader.txt```

from step 4 ```fruitbat/refgenome/fasta/output/mito_fixed.fna.txt```

from step 5 ```fruitbat/refgenome/fasta/output/combined_mito_bat.fna.gz```

## Preparing reference GTF for mapping
Because reference GTF only contains genomic genes

Created the mitochondrial GTF manually getting the information from https://www.ncbi.nlm.nih.gov/nuccore/NC_026542.1?report=genbank
[mitogenes gtf](fruitbat/refdata_vampyrus_refseq/mitogenes.gtf) and using the baboon GTF lines (changing gene ID, coordinates ect). 

Output: ```fruitbat/refgenome/gtf/output/cleanmitogenes.gtf```

## Filtering fruit bat gtf for protein-coding genes
Because CellRanger *mkgtf* command cannot handle RefSeq reference gtf files
### Get gene IDs of protein-coding genes from gtf file
Get the gene IDs of genes that are protein-coding from a gtf file from RefSeq

Input: ```GCF_000151845.1_Pvam_2.0_genomic.gtf```

Output: ```fruitbat/refgenome/gtf/output/idlist.txt```

Run it like: 

```bash
python get_gene_ids.py genomic.gtf idlist.txt
```

### Filter gtf based on ID list created in previous step
Filter gtf file from RefSeq so output gtf only contains protein-coding genes. Using gene ID list of protein-coding genes from previous step to filter.

Input: ```GCF_000151845.1_Pvam_2.0_genomic.gtf```

```fruitbat/refgenome/gtf/output/idlist.txt```

Output:
```fruitbat/refgenome/gtf/output/genomicfixed.gtf```

Run it like: 

```bash
python filtering_for_cellranger.py genomic.gtf idlist output.gtf
```

### Combine mitochondrial and genomic gtf files
Concatonate the same was as fasta files before via cat.

Output:
```fruitbat/refgenome/gtf/output/fruitbat.gtf```