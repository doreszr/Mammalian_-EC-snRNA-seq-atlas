#!/usr/bin/python3

### Author: Dorottya
### Date: 06/02-24

### filter out genes that were kept from previous assembly with old gene id so only keep gene ids that are present in the gtf file
### run it like: python filter_faa_files.py x.fa filtered.fa

import sys
import pandas as pd

fa_depo = dict()
genes = []
line_depo = []
counter = 0


with open(sys.argv[1], 'r') as gtf, open(sys.argv[2], 'r') as fa, open(sys.argv[3], 'w') as fout, open(sys.argv[4], 'w') as out:
    for line in fa:
        if line.startswith('>'):
            counter += 1
            line = line.rstrip()
            if len(line_depo) == 0:
               key = line
            else:
                fa_depo[key] = line_depo
                line_depo = []
                key = line
        else:
            line = line.rstrip()
            line_depo.append(line)
    fa_depo[key] = line_depo
    # print(fa_depo)

    for lane in gtf:
        lane = lane.rstrip()
        if lane.startswith('#'):
            continue
        else:
            column = lane.split('\t')
            # print(column[2])
            if column[2] == 'gene':
                info = column[8].split(';')
                name = info[0].split('\"')
                # print(name[1])
                genes.append(name[1])

    for key, value in fa_depo.items():
        geneid = key.lstrip('>')
        geneid = geneid.split('.')
        if geneid[0] in genes:
            seq = '\n'.join(value)
            print(key + '\n' + seq, file = fout)
        else:
            print(geneid[0], file = out)
            
    print(counter)
    print(len(fa_depo))
        

