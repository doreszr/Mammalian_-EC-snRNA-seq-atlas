#!/usr/bin/python3

### Author: Dorottya
### Date: 06/02-24

### Ensemble Salmo salar assembly does not contain mitochondrial sequence
### but previous version of the Ref Seq reference genome contains the mitochondrial genes GCF_000233375.1_ICSASG_v2
### This script filters out the mitochondrial proteins from proteome.

import sys

# set up empty dictionary
prot_dict = dict()
# set up empty list for saving protein sequence
seq = []

# open input and output files: previous assembly faa as input and filtered faa as output
with open(sys.argv[1], 'r') as mitoprot, open(sys.argv[2], 'w') as fout:
    for line in mitoprot:
        line = line.rstrip()
        # catch header
        if line.startswith('>'):
            # check if there is anything saved in seq
            if len(seq) > 0:
                # save previously saved protein sequence
                prot_dict[key] = seq
                # empty the list again
                seq = []
                # save the new header
                key = line
            else:
                # save the new header
                key = line
        else:
            # save protein sequence
            seq.append(line)
    
    # iterate through the dictionary
    for key, value in prot_dict.items():
        # catch the mitochondrial proteins
        if 'mitochondrion' in key:
            # print those proteins with header and sequence into output file
            print(key, file = fout)
            print(*value, sep = '\n', file = fout)
                

