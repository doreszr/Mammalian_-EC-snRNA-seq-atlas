#!/usr/bin/python3

### Author: Dorottya
### Date: 23/09-24

### extracting many2many and one2many from OrthoFinder species_vs_species tsv output 
### run it like: python get_prot_ids.py gtfspecieshuman tsv gtfspeciesmouse

import sys
import pandas as pd

# setting up empty dictionaries to store

# protein ids human
human_prot = dict()

# one2one orthologs
orthomh1 = dict()

# mouse protein ids
mouse_prot = dict()


dict_list = []

# setting up data frame to store human and mouse gene ids
orthomhmany = pd.DataFrame(columns=[ 'geneid_h', 'geneid_m'])

# open human gtf, orthofinder ortholog file, mouse gtf
with open(sys.argv[1], 'r') as gtfh, open(sys.argv[2], 'r') as orthos, open(sys.argv[3], 'r') as gtfm:
    # creating a RefSeq gene ID gene name dictionary
    for lines in gtfh:
        lines = lines.rstrip()
        # skip header lines
        if lines.startswith('#'):
            continue
        # split into columns    
        else:
            column = lines.split('\t')
            # extract gene ID and name
            if column[2] == 'gene':
                info = column[8].split(';')
                gene_id = info[0].split('\"')
                gene_name = info[2].split('\"')
                # print(gene_nameq)
                # accommodating unconsistent formatting
                if 'gene_source' in gene_name[0]:
                    # print(info)
                    key = gene_id[1]
                    human_prot[key] = gene_id[1]
                else:
                    key = gene_id[1]
                    human_prot[key] = gene_name[1]

    # extracting one2one and many orthologs    
    for lines in orthos:
        line = lines.rstrip()
        # catching the orthologs
        if line.startswith('OG'):
            # checking if one2one ortholog
            if ',' in line:
                line = line.split('\t')
                # figuring out one2many/many2many/many2one mouse-human
                # many2..
                if ',' in line[1]:
                    mouse = line[1].split(',')
                    if ',' in line[2]:
                        # many2many mouse-human
                        human = line[2].split(',')
                        # saving each mouse gene id in that line
                        for genes in mouse:
                            # extracting mouse gene ids
                            genes = genes.lstrip()
                            genes = genes.split('.')
                            for prot in human:
                                # extracting human gene ids
                                prot = prot.lstrip()
                                prot = prot.split('.')
                                # populating dictionary
                                row_dict = {'geneid_h': prot[0], 'geneid_m': genes[0]}
                                dict_list.append(row_dict)
                    # many2one mouse-human 
                    else:
                        for genes in mouse:
                            genes = genes.lstrip()
                            genes = genes.split('.')
                            prot = line[2].split('.')
                            row_dict = {'geneid_h': prot[0], 'geneid_m': genes[0]}
                            dict_list.append(row_dict)
                # one2many mouse-human
                else:
                    human = line[2].split(',')
                    for prot in human:
                            genes = line[1].split('.')
                            prot = prot.lstrip()
                            prot = prot.split('.')
                            row_dict = {'geneid_h': prot[0], 'geneid_m': genes[0]}
                            dict_list.append(row_dict)       
            # one2one ortholog (no comma in line)
            else:
                # print(line)
                columns = line.split('\t')
                mouse = columns[1].split('.')
                # print(human)
                # print(human[0])
                human = columns[2].split('.')
                key = human[0]
                orthomh1[key] = mouse[0]
    # create data frame from dictionary 
    orthomhmany = pd.DataFrame.from_dict(dict_list)

    # print(orthohb)


    # extracting gene ids from mouse gtf for matching protein ids
    for lines in gtfm:
        lines = lines.rstrip()
        # skipping header
        if lines.startswith('#'):
            continue
        else:
            column = lines.split('\t')
            if column[2] == 'gene':
                info = column[8].split(';')
                gene_id = info[0].split('\"')
                gene_name = info[2].split('\"')
                # print(gene_nameq)
                # accommodating unconsistent formatting
                if 'gene_source' in gene_name[0]:
                    # print(info)
                    key = gene_id[1]
                    mouse_prot[key] = gene_id[1]
                else:
                    key = gene_id[1]
                    mouse_prot[key] = gene_name[1]
    # print(human_prot)



    # print(id_prot)

    # building the data frames
    df_human = pd.DataFrame(human_prot.items(), columns=['geneid_h','gene_h'])
    df_mouse = pd.DataFrame(mouse_prot.items(), columns=['geneid_m','gene_m'])
    df_ortho1 = pd.DataFrame(orthomh1.items(), columns=['geneid_h','geneid_m'])
    combi1one = pd.merge(df_human, df_ortho1, how="right", on="geneid_h")
    combi2one = pd.merge(combi1one, df_mouse, how="left", on="geneid_m")
    combi1many = pd.merge(df_human, orthomhmany, how="right", on="geneid_h", validate="one_to_many")
    combi2many = pd.merge(combi1many, df_mouse, how="left", on="geneid_m", validate="many_to_one")
    print(combi1many)
    print(df_ortho1)
    print(combi2many)
    combi2one.to_csv('ortho_mh_one.csv', index=False)
    combi2many.to_csv('ortho_mh_many.csv', index=False)
    df_human.to_csv('df_human.csv', index=False)
    df_mouse.to_csv('df_mouse.csv', index=False)
    orthomhmany.to_csv('orthomhmany.csv', index=False)