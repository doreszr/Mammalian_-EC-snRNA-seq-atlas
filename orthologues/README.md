# Identify orthologoues genes with OrthoFinder
Standard workflow was followed from [OrthoFinder tutorial](https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html).

## Download proteomes from Ensembl
For the following species:
- Mus musculus (mouse)
- Rattus norvegicus (rat)
- Danio rerio (zebra fish)
- Homo sapiens (human)
- Pan paniscus (bonobo)
- (chimpanzee)
- Gorilla gorilla (gorilla)
- Macaca mulatta (macaque)
- Papio anubis (baboon)
- Pteropus vampyrus (fruit bat)
- Salmo salar (salmon)

## Process proteomes
### Filter out mitochondrial proteins from previous genomic assembly - salmon
Because reference genome faa only contains genomic proteins but previous version of the reference proteome contains the mitochondrial genes GCF_000233375.1_ICSASG_v2.

Input: ```salmon/protein_salmonv2.faa```

Output: ```orthologues/process_proteomes/mito/salmon/output/salmon_mito.faa```

Run it like: 

```bash
python orthologues/process_proteomes/mito/salmon/scripts/salmon_mitoprot.py protein_salmonv2.faa salmon_mito.faa
```

### Create mitochondrial proteome for fruit bat
Download the protein sequences from RefSeq one by one (13 proteins). Same mitochondrion assembly was used as for mitochondrial GTF ([NC_026542.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_026542.1))

Output: ```orthologues/process_proteomes/mito/fruit_bat/output/mito_bat.faa```

### Run primary_transcript.py on proteomes
OrthoFinder provied a script to filter proteome files to extract just the longest transcript variant per gene and run OrthoFinder on these files.

Input: ```reference_proteomes/species/fa```

Output: 
vertebrates ```orthologues/process_proteomes/primary_transcripts/output```

primates ```orthologues/process_proteomes/primary_transcripts/output/primate```

Run it like: 

```bash
for f in *fa ; do python orthologues/process_proteomes/primary_transcripts/scripts/primary_transcript.py $f ; done
```

### Combine genomic and mitochondrial proteomes
Pteropus vampyrus and Salmo salar proteomes only contained genomic proteins so mitochondrial proteome was added.

Concatonate to genomic faa file via cat (```bash cat *faa > combined.faa```)

Input:
salmon ```orthologues/process_proteomes/primary_transcripts/salmon.faa```
```orthologues/process_proteomes/mito/salmon/output/salmon_mito.faa```

fruit bat ```orthologues/process_proteomes/primary_transcripts/bat.faa```
```orthologues/process_proteomes/mito/fruit_bat/output/mito_bat.faa```

Output:
salmon ```orthologues/process_proteomes/combined/salmon/output/combined_salmon.faa```

fruit bat ```orthologues/process_proteomes/combined/fruit_bat/output/combined_bat.faa```

### Filter out old protein name duplicates
Filter out genes that were kept from previous assembly with old gene id so only keep gene ids that are present in the GTF file.

Input: ```orthologues/process_proteomes/primary_transcripts/output/human.faa```
```orthologues/process_proteomes/primary_transcripts/output/mouse.faa```
```orthologues/process_proteomes/primary_transcripts/output/danio.faa```

Output:  ```orthologues/process_proteomes/filter/output/human_filt.fa```
```orthologues/process_proteomes/filter/output/mouse_filt.fa```
```orthologues/process_proteomes/filter/output/danio_filt.fa```

Run it like: 

```bash
python orthologues/process_proteomes/filter/scripts/filter_faa_files.py human.fa filtered.fa
```

## Run OrthoFinder
To find the orthologues genes between baboon and human, there was a primate run seperately from vertebrate species spanning from fish to human.
### Primates
Species:
- Homo sapiens
- Pan paniscus
- Pan troglodytes
- Gorilla gorilla
- Macaca mulatta
- Papio anubis

Input:
human ```orthologues/process_proteomes/filter/output/human_filt.fa```
rest ```orthologues/process_proteomes/primary_transcripts/output/primate```

Output: ```orthologues/OrthoFinder/primates/Results```

### Vertebrate species
Species:
- Homo sapiens
- Papio anubis
- Mus musculus
- Rattus norvegicus
- Danio rerio
- Salmo salar
- Pteropus vampyrus

Input:
human, mouse, zebra fish ```orthologues/process_proteomes/filter/output/```
salmon ```orthologues/process_proteomes/combined/salmon/output/combined_salmon.faa```
fruit bat ```orthologues/process_proteomes/combined/fruit_bat/output/combined_bat.faa```
baboon, rat ```orthologues/process_proteomes/primary_transcripts/output```

Output: ```orthologues/OrthoFinder/vertebrates/Results_Feb06```

Run it like: 

```bash
orthofinder -f primary_transcripts/primates/
```

## Mapping orthologous genes between target species

### Human - Baboon
#### Build one2one many ortholog data frames
From OrthoFinder output which uses protein IDs for baboon and Ensembl gene IDs for human. Fetch matching RefSeq gene IDs from GTF files.

Input: ```baboon/refgenome/output/0705baboon.gtf```
```/orthologues/OrthoFinder/primates/Results_Feb06/Orthologues/Orthologues_human_filt/human_filt__v__baboon.tsv```
```reference_genomes/human/Homo_sapiens.GRCh38.111.gtf```

Output: ```orthologues/mappingortho/baboon_vs_human/output/ortho_hb_one.csv```
```orthologues/mappingortho/baboon_vs_human/output/ortho_hb_many.csv```
```orthologues/mappingortho/baboon_vs_human/output/df_baboon.csv```
```orthologues/mappingortho/baboon_vs_human/output/df_human.csv```
```orthologues/mappingortho/baboon_vs_human/output/orthohbmany.csv```

Run it like: 
```bash
python orthologues/mappingortho/baboon_vs_human/scripts/get_prot_ids_hb.py baboon.gtf human_filt__v__baboon.tsv Homo_sapiens.GRCh38.111.gtf
```

#### Find the orthologues for genes that are found in our data
Translate baboon genes to their human orthologues. One2one are straight forward. Manys: when only one of the manys were found in our data (human or baboon) it was picked as one2one, but when multiple were found the one with a higher expression rate was picked. Idea from BENGAL pipeline.

Output: ```orthologues/mappingortho/baboon_vs_human/output/allorthobh.csv```

Run the R script ```orthologues/mappingortho/baboon_vs_human/scripts/mappingorthobhmapping.R``` 

### Mouse - Human
#### Build one2one many ortholog data frames
From OrthoFinder output which uses Ensembl gene IDs . Fetch matching RefSeq gene IDs from GTF files.

Input: ```reference_genomes/human/Homo_sapiens.GRCh38.111.gtf```
```orthologues/OrthoFinder/vertebrates/Results/Orthologues/Orthologues_human_filt/mouse_filt__v__human_filt.tsv```
```reference_genomes/mouse/Mus_musculus.GRCm39.111.gtf```

Output: ```orthologues/mappingortho/human_vs_mouse/output/ortho_mh_one.csv```
```orthologues/mappingortho/human_vs_mouse/output/ortho_mh_many.csv```
```orthologues/mappingortho/human_vs_mouse/output/df_mouse.csv```
```orthologues/mappingortho/human_vs_mouse/output/df_human.csv```
```orthologues/mappingortho/human_vs_mouse/output/orthomhmany.csv```

Run it like: 

```bash
python orthologues/mappingortho/human_vs_mouse/scripts/get_prot_ids_mh.py Homo_sapiens.GRCh38.111.gtf mouse_filt__v__human_filt.tsv Mus_musculus.GRCm39.111.gtf
```

#### Find the orthologues for genes that are found in our data
Translate human genes to their mouse orthologues. One2one are straight forward. Manys: when only one of the manys were found in our data (human or mouse) it was picked as one2one, but when multiple were found the one with a higher expression rate was picked. Idea from BENGAL pipeline.

Output: ```orthologues/mappingortho/human_vs_mouse/output/allorthohm.csv```

primate integration ```orthologues/mappingortho/human_vs_mouse/output/allorthohm_primate.csv```

Run the R script ```orthologues/mappingortho/human_vs_mouse/scripts/mappingorthomousehuman.R``` 

or for primate integration ```orthologues/mappingortho/human_vs_mouse/scripts/mappingorthomousehuman_primate.R```

### Mouse - Fruit bat
#### Build one2one many ortholog data frames
From OrthoFinder output which uses protein IDs for bat and Ensembl gene IDs for mouse. Fetch matching RefSeq gene IDs from GTF files.

Input: ```fruitbat/refgenome/gtf/output/fruitbat.gtf```
```Results_Feb06/Orthologues/Orthologues_mouse_filt/mouse_filt__v__bat.tsv```
```reference_genomes/mouse/Mus_musculus.GRCm39.111.gtf```

Output: ```orthologues/mappingortho/bat_vs_mouse/output/ortho_mb_one.csv```
```orthologues/mappingortho/bat_vs_mouse/output/ortho_mb_many.csv```
```orthologues/mappingortho/bat_vs_mouse/output/df_mouse.csv```
```orthologues/mappingortho/bat_vs_mouse/output/df_bat.csv```
```orthologues/mappingortho/bat_vs_mouse/output/orthombmany.csv```

Run it like: 

```bash
python orthologues/mappingortho/bat_vs_mouse/scripts/get_prot_ids_mb.py fruitbat.gtf mouse_filt__v__bat.tsv Mus_musculus.GRCm39.111.gtf
```

#### Find the orthologues for genes that are found in our data
Translate bat genes to their mouse orthologues. One2one are straight forward. Manys: when only one of the manys were found in our data (bat or mouse) it was picked as one2one, but when multiple were found the one with a higher expression rate was picked. Idea from BENGAL pipeline.

Run the R script ```orthologues/mappingortho/bat_vs_mouse/scripts/mappingorthomousebat.R``` 
