#!/usr/bin/python3
### Author: Dorottya
### Date: 030524

### filter fruit bat genomic gtf file from RefSeq GCF_000151845.1_Pvam_2.0_genomic.gtf (mito part is missing and will be added later)
### so output only contains protein coding genes
### run it like: python filtering_gtf.py gtf.gtf idlist output.gtf
### using gene ID list of protein coding genes to filter


import sys
# create empty list to store gene ids
genelist = []
counter = 0



# open input and output files: gtf gene id list and filtered gtf as output
with open(sys.argv[1], 'r') as gtf, open(sys.argv[2], 'r') as geneid,  open(sys.argv[3], 'w') as fout:
    # looping through the gene id list to remove \n character and add them to a list
    for genes in geneid:
        genes = genes.rstrip()
        genelist.append(genes)

    # looping through the gtf file
    for lanes in gtf:
        lane = lanes.rstrip()
        # print header lines
        if lane.startswith('#'):
            # skip closing line
            if lane == '###':
                continue
            else:
                print(lane, file = fout)
        else:
            
            # splitting lines so we can extract gene ids
            column = lane.split('\t')
            # chatching gene lines
            if column[2] == 'gene':
                info = column[8].split(';')
                idl = info[1].split(':')
                # print(idl)
                
                id = idl[1]
                id = id.rstrip('\"')
                # print(id)
                
                # grabbing protein-coding genes
                if id in genelist:
                    counter += 1
                    print(lane, file = fout)
                else:
                    continue
            
            else:
                # check other position so no protien coding genes are skipped
                info = column[8].split(';')
                idl = info[2].split(':')
                
                
                id = idl[1]
                id = id.rstrip('\"')
                
                # grabbing protein-coding genes
                if id in genelist:
                    counter += 1
                    print(lane, file = fout)
                else:
                    continue
    # to see number of protein coding genes        
    print(counter)