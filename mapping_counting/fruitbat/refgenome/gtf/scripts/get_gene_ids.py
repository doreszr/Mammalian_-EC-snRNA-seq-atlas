#!/usr/bin/python3
### Author: Dorottya
### Date: 230524

### get the gene ids of genes that are protein coding from fruit bat gtf GCF_000151845.1_Pvam_2.0_genomic.gtf (mito part is missing and will be added later)
### run it like: python get_gene_ids_bat.py gtf idlist


import sys

# setting up list to save the ids into
genelist = []
# setting up counter to be able to quickly check if parser catches all the protein coding genes
genecount = 0

# open input gtf and ouput file
with open(sys.argv[1], 'r') as gtf, open(sys.argv[2], 'w') as fout:
    # loop through gtf file
    for lanes in gtf:
        lane = lanes.rstrip()
        # skip header
        if not lane.startswith('#'):
            # split line into columns
            column = lane.split('\t')
            
            # chatching gene lines only
            if column[2] == 'gene':
                # extracting gene id
                info = column[8].split(';')
                # accomodate of all line structures (different number of columns in info)
                gene_type = info[-2]
                gene_typemore = info[-3]
                gene_typemoremore = info[-4]
                gene_typemoremoremore = info[-5]

            # print(gene_typemoremore)
               
                # saving protein coding gene ids into list
                if gene_type == ' gene_biotype \"protein_coding\"':
                    idl = info[1].split(':')
                    # print(idl)
                    id = idl[1]
                    id = id.rstrip('\"')
                    genelist.append(id)
                    # keeping the gene count updated
                    genecount += 1
                    
                elif gene_typemore == ' gene_biotype \"protein_coding\"':
                    idl = info[1].split(':')
                    id = idl[1]
                    id = id.rstrip('\"')
                    genelist.append(id)
                    # keeping the gene count updated
                    genecount += 1
                
                elif gene_typemoremore == ' gene_biotype \"protein_coding\"':
                    idl = info[1].split(':')
                    id = idl[1]
                    id = id.rstrip('\"')
                    genelist.append(id)
                    # keeping the gene count updated
                    genecount += 1
                
                else:
                    continue
                    

    # print gene list into one column
    print(*genelist, sep = '\n', file = fout)
    # print final count for monitoring
    print(genecount)