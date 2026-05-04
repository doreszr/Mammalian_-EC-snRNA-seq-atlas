# Mapping and counting fruit bat snRNA-seq data
## Data
4 fruit bat (Rousettus aegyptiacus) EC samples. More info about sample preparation at the raw data folders.

## Mapping and counting with cellranger
Mapping and counting the reads from 10X experiment fruit bat EC samples with 10x Genomics Cell Ranger 7.1.0.
Output is standard CellRanger output
### Run CellRanger mkref
Indexing step: creating index for mapping from filtered reference GTF file and reference FASTA

Input: ```fruitbat/refgenome/fasta/output/combined_mito_bat.fna.gz```
```fruitbat/refgenome/gtf/output/fruitbat.gtf```

Output:  ```fruitbat/cellranger/indexing/P_vapmpyrus_genome```

### Run cellRanger count
Mapping and counting step: alignment of reads and counting transcripts

Input: 
the paired-end fastq files found at ```fruitbat/sample ID```

the generated index ```fruitbat/cellranger/indexing/P_vapmpyrus_genome```

Output: 
CellRanger output for each sample found at ```fruitbat/cellranger/mapping/sample ID/outs```

