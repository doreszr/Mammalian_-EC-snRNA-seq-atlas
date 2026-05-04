
# Mapping and counting baboon snRNA-seq data
## Data
4 baboon (Papio hamadryas) EC samples. More info about sample preparation at the raw data folders.

## Mapping and counting with CellRanger
Mapping and counting the reads from 10X experiment baboon EC samples with 10x Genomics Cell Ranger 7.1.0.
Output is standard CellRanger output
### Run CellRanger mkref
Indexing step: creating index for mapping from filtered reference GTF file and reference FASTA

Input:
```GCF_008728515.1_Panubis1.0_genomic.fna```
```baboon.gtf```

Output:
```baboon/cellranger/indexing/Papio_anubis_genome```

### Run CellRanger count
Mapping and counting step: alignment of reads and counting transcripts

Input: 
the paired-end fastq files found at ```baboon/sample ID```

the generated index ```baboon/cellranger/indexing/Papio_anubis_genome```

Output:
CellRanger output for each sample found at ```baboon/cellranger/mapping_counting/sample ID/outs``` 
