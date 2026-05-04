#!/usr/bin/python3
### Author: Dorottya
### Date: 061223

### get the gene IDs of genes that are protein coding from a gtf file from RefSeq
### run it like: python get_gene_ids.py gtf idlist


import sys

# setting up list to save the ids into
genelist = []
# setting up counter to be able to quickly check if parser catches all the protein coding genes
genecount = 0

# open input gtf and ouput list file
with open(sys.argv[1], 'r') as gtf, open(sys.argv[2], 'w') as fout:
    # loop through gtf file
    for lanes in gtf:
        lane = lanes.rstrip()
        # skip header
        if not lane.startswith('#'):
            # print(lane)
            column = lane.split('\t')
            # print(column[2])
            
            # chatching gene lines only
            if column[2] == 'gene':
                # extracting gene id
                info = column[8].split(';')
                # print(info[-2])
                # accomodate of all line structures (different number of columns in info)
                gene_type = info[-2]
                gene_typemore = info[-3]
                gene_typemoremore = info[-4]
                gene_typemoremoremore = info[-5]
                # print(gene_type)
                # saving protein coding gene ids into list
                if gene_type == ' gene_biotype \"protein_coding\"':
                    name = info[0].split('\"')
                    idl = info[2].split(':')
                    id = idl[1]
                    id = id.rstrip('\"')
                    genelist.append(id)
                    genecount += 1
                elif gene_typemore == ' gene_biotype \"protein_coding\"':
                    name = info[0].split('\"')
                    idl = info[2].split(':')
                    id = idl[1]
                    id = id.rstrip('\"')
                    genelist.append(id)
                    genecount += 1
                elif gene_typemoremore == ' gene_biotype \"protein_coding\"':
                    name = info[0].split('\"')
                    idl = info[2].split(':')
                    id = idl[1]
                    id = id.rstrip('\"')
                    genelist.append(id)
                    genecount += 1
                elif gene_typemoremoremore == ' gene_biotype \"protein_coding\"':
                    name = info[0].split('\"')
                    idl = info[2].split(':')
                    id = idl[1]
                    id = id.rstrip('\"')
                    genelist.append(id)
                    genecount += 1
                else:
                    continue

    # print gene list into one column
    print(*genelist, sep = '\n', file = fout)
    # print final count for monitoring
    print(genecount)