#!/usr/bin/python3
### Author: Dorottya
### Date: 10122023

### filter gtf file from RefSeq GCF_008728515.1 so output gtf only contains protein coding genes
### run it like: python filtering_for_cellranger.py gtf.gtf idlist output.gtf
### using gene ID list of protein coding genes to filter
### special case for mt genes - fixing those lines as well

import sys
# create empty list to store gene ids
genelist = []

# duplicate genes
dupligenes = ["100137232", "100137263", "100137470", "100337524"]

# open input files: gtf and gene id list
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
            print(lane, file = fout)
        else:
            # splitting lines so we can extract gene ids
            column = lane.split('\t')
            info = column[8].split(';')
            name = info[0].split('\"')
            idl = info[2].split(':')
            # print(len(idl))
            
            # regular genes
            if len(idl) == 2:
                id = idl[1]
                id = id.rstrip('\"')
                # print(id)

                # grabbing protein-coding genes
                if id in genelist:
                    # skip duplicate gene ids
                    if id in dupligenes:
                        continue
                    else:
                         # printing the lines that are protien coding genes
                        print(lane, file = fout)
            # irregular genes (missing gene id from non-gene lines in gtf)
            elif len(idl) < 2:
                # extracting gene name
                column = lane.split('\t')
                info = column[8].split(';')
                name = info[0].split('\"')
                # print(name)

                # catching mt genes and printing out fixed lines
                if name[1] == 'CYTB':
                    info[1] =  ' transcript_id \"unassigned_transcript_634\"; db_xref \"GeneID:14444641\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ATP6':
                    info[1] =  ' transcript_id \"unassigned_transcript_621\"; db_xref \"GeneID:14444634\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ATP8':
                    info[1] =  ' transcript_id \"unassigned_transcript_620\"; db_xref \"GeneID:14444633\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ND1':
                    info[1] =  ' transcript_id \"unassigned_transcript_605\"; db_xref \"GeneID:14444642\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ND2':
                    info[1] =  ' transcript_id \"unassigned_transcript_609\"; db_xref \"GeneID:14444630\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ND3':
                    info[1] =  ' transcript_id \"unassigned_transcript_624\"; db_xref \"GeneID:14444636\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ND4':
                    info[1] =  ' transcript_id \"unassigned_transcript_627\"; db_xref \"GeneID:14444638\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ND4L':
                    info[1] =  ' transcript_id \"unassigned_transcript_626\"; db_xref \"GeneID:14444637\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ND5':
                    info[1] =  ' transcript_id \"unassigned_transcript_631\"; db_xref \"GeneID:14444639\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'ND6':
                    info[1] =  ' transcript_id \"unassigned_transcript_632\"; db_xref \"GeneID:14444640\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'COX1':
                    info[1] =  ' transcript_id \"unassigned_transcript_615\"; db_xref \"GeneID:14444631\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'COX2':
                    info[1] =  ' transcript_id \"unassigned_transcript_618\"; db_xref \"GeneID:14444632\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                elif name[1] == 'COX3':
                    info[1] =  ' transcript_id \"unassigned_transcript_622\"; db_xref \"GeneID:14444635\"'
                    # print(info)
                    info = ';'.join(info)
                    column[8] = info
                    line = '\t'.join(column)
                    print(line, file = fout)
                # skipping non-protein-coding genes
                else:
                    continue
            