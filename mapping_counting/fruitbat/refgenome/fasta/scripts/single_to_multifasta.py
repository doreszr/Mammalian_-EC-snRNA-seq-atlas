#!/usr/bin/python3
### Author: Dorottya
### Date: 2023-05-12

### This program will convert a single-line fasta file into multi-line fasta file.
### run it like: python single_to_multifasta.py gtf idlist output.gtf

##################
import sys

# setting up empty counter for counting nucleotides
count = 0

# setting up empty lists for saving each line
linelist = []
# and saving the multi-line fasta
alllines = []

# setting up counter to catch the last line not reaching max line lengh
i = 0

# opening input and output files
with open(sys.argv[1], 'r') as fin, open(sys.argv[2], 'w') as fout:
    # opening the line
    for line in fin:
        # checking if input fasta is truly single-line
        # print('test')
        line=line.rstrip()
        # iterating through the nucleotides in the line
        for element in line:
            # keeping a count of nucleotides
            i += 1
            # when we reach end of line
            if count == 80:
                # add line to allline list
                alllines.append(''.join(linelist))
                # setting back counter
                count = 1
                # emptying line list
                linelist = []
                # saving nucleotides for next line
                linelist.append(element)
            # catching the last nucleotide
            elif i == len(line):
                # saving nucleotides for next line
                linelist.append(element)
                # saving that line to the allline list
                alllines.append(''.join(linelist))
            # saving nucleotides in a line
            else:
                # keeping a count of nucleotides
                count += 1
                # saving nucleotides
                linelist.append(element)
    # lines are joined by new line chr and are printed into output file
    print('\n'.join(alllines) + '\n', file = fout)
