#!/usr/bin/python3

### Author: Dorottya
### Date: 07/02-24

### extracting many2many and one2many from OrthoFinder species_vs_species tsv output 
### run it like: python get_prot_ids.py gtfspeciesbaboon tsv gtfspecieshuman

import sys
import pandas as pd

# setting up empty dictionaries to store

# protein ids baboon
id_prot = dict()

# one2one orthologs
orthohb1 = dict()

# human protein ids
human_prot = dict()


dict_list = []

# setting up data frame to store baboon protein ids and human gene ids
orthohbmany = pd.DataFrame(columns=[ 'protid_b', 'geneid_h'])

# open baboon gtf, orthofinder ortholog file, human gtf
with open(sys.argv[1], 'r') as gtf, open(sys.argv[2], 'r') as orthos, open(sys.argv[3], 'r') as gtfh:
    for lanes in gtf:
        lane = lanes.rstrip()
        # skip header lines
        if lane.startswith('#'):
            continue
        # split into columns
        else:
            column = lane.split('\t')
            # print(column[2])
            # catch the CDS lines (containing protein ids)
            if column[2] == 'CDS':
                # splitting last column (containing information)
                info = column[8].split(';')
                name = info[0].split('\"')
                # print(name)
                # getting the protein id
                protidl = info[6].split('\"')
                # print(protidl)

                # accommodating for different info collumns
                if 'protein_id' not in protidl[0]:
                    # print(protidl)
                    if 'protein_id' not in protidl[0]:
                        protidl = info[7].split('\"')
                        if 'protein_id' not in protidl[0]:
                            protidl = info[8].split('\"')
                            if 'protein_id' not in protidl[0]:
                                protidl = info[9].split('\"')
                                if 'protein_id' not in protidl[0]:
                                    protidl = info[10].split('\"')
                                    if 'protein_id' not in protidl[0]:
                                        protidl = info[11].split('\"')
                                        # print(protidl)
                                        protid = protidl[1]
                                        key = protid
                                        id_prot[key] = name[1]
                                    else:
                                        # print(protidl)
                                        protid = protidl[1]
                                        key = protid
                                        id_prot[key] = name[1] 
                                else:
                                    # print(protidl)
                                    protid = protidl[1]
                                    key = protid
                                    id_prot[key] = name[1] 
                            else:
                                # print(protidl)
                                protid = protidl[1]
                                key = protid
                                id_prot[key] = name[1] 
                        else:
                             # print(protidl)
                            protid = protidl[1]
                            key = protid
                            id_prot[key] = name[1]
                    else:
                             # print(protidl)
                            protid = protidl[1]
                            key = protid
                            id_prot[key] = name[1]
                else:
                    protid = protidl[1]
                    key = protid
                    id_prot[key] = name[1]

    # extracting one2one and many orthologs    
    for lines in orthos:
        line = lines.rstrip()
        # catching the orthologs
        if line.startswith('OG'):
            # checking if one2one ortholog
            if ',' in line:
                line = line.split('\t')
                # figuring out one2many/many2many/many2one human-baboon
                # many2..
                if ',' in line[1]:
                    human = line[1].split(',')
                    # many2many human-baboon
                    if ',' in line[2]:
                        baboon = line[2].split(',')
                        # saving each human gene id in that line
                        for genes in human:
                            # extracting human gene ids
                            genes = genes.lstrip()
                            genes = genes.split('.')
                            # saving baboon protein id to each human gene id in the line
                            for prot in baboon:
                                # extracting baboon protein ids
                                prot = prot.lstrip()
                                # populating dictionary
                                row_dict = {'protid_b': prot, 'geneid_h': genes[0]}
                                dict_list.append(row_dict)
                    # many2one human-baboon            
                    else:
                        for genes in human:
                            genes = genes.lstrip()
                            genes = genes.split('.')
                            row_dict = {'protid_b': line[2], 'geneid_h': genes[0]}
                            dict_list.append(row_dict)
                # one2many human-baboon            
                else:
                    baboon = line[2].split(',')
                    for prot in baboon:
                            genes = line[1].split('.')
                            prot = prot.lstrip()
                            row_dict = {'protid_b': prot, 'geneid_h': genes[0]}
                            dict_list.append(row_dict)
            # one2one ortholog (no comma in line)                       
            else:
                # print(line)
                columns = line.split('\t')
                human = columns[1].split('.')
                # print(human)
                # print(human[0])
                key = columns[2]
                orthohb1[key] = human[0]
    # create data frame from dictionary            
    orthohbmany = pd.DataFrame.from_dict(dict_list)

    # print(orthohb)

    # extracting gene ids from baboon gtf for matching protein ids
    for lines in gtfh:
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
                # print(gene_name)
                # accommodating unconsistent formatting
                if 'gene_source' in gene_name[0]:
                    # print(info)
                    key = gene_id[1]
                    human_prot[key] = gene_id[1]
                else:
                    key = gene_id[1]
                    human_prot[key] = gene_name[1]
    # print(human_prot)



    # print(id_prot)

    
    df_baboon = pd.DataFrame(id_prot.items(), columns=['protid_b','gene_b'])
    df_human = pd.DataFrame(human_prot.items(), columns=['geneid_h','gene_h'])
    df_ortho1 = pd.DataFrame(orthohb1.items(), columns=['protid_b','geneid_h'])
    combi1one = pd.merge(df_baboon, df_ortho1, how="right", on="protid_b")
    combi2one = pd.merge(combi1one, df_human, how="left", on="geneid_h")
    combi1many = pd.merge(df_baboon, orthohbmany, how="right", on="protid_b", validate="one_to_many")
    combi2many = pd.merge(combi1many, df_human, how="left", on="geneid_h", validate="many_to_one")
    print(combi1many)
    print(df_ortho1)
    print(combi2many)
    print(combi2one)
    combi2one.to_csv('ortho_hb_one.csv', index=False)
    combi2many.to_csv('ortho_hb_many.csv', index=False)
    df_baboon.to_csv('df_baboon.csv', index=False)
    df_human.to_csv('df_human.csv', index=False)
    orthohbmany.to_csv('orthohbmany.csv', index=False)

    

