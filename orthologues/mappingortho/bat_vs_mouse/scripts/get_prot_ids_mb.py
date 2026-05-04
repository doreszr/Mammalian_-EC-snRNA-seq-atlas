#!/usr/bin/python3

### Author: Dorottya
### Date: 28/05-24

### extracting many2many and one2many from OrthoFinder species_vs_species tsv output 
### run it like: python get_prot_ids.py gtfspeciesbat mouse_vs_bat.tsv gtfspeciesmouse

import sys
import pandas as pd

# setting up empty dictionaries to store

# protein ids bat
id_prot = dict()

# one2one orthologs
orthombat1 = dict()

# mouse protein ids
mouse_prot = dict()


dict_list = []

# setting up data frame to store bat protein ids and mouse gene ids
orthombatmany = pd.DataFrame(columns=[ 'protid_bat', 'geneid_m'])

# open bat gtf, orthofinder ortholog file, mouse gtf
with open(sys.argv[1], 'r') as gtf, open(sys.argv[2], 'r') as orthos, open(sys.argv[3], 'r') as gtfm:
    # extracting baboon protein IDs and gene names from gtf file
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
                # continue
                if ',' in line:
                    line = line.split('\t')
                    # figuring out one2many/many2many/many2one mouse-bat
                    # many2..
                    if ',' in line[1]:
                        mouse = line[1].split(',')
                        # many2many mouse-bat
                        if ',' in line[2]:
                            bat = line[2].split(',')
                            # saving each mouse gene id in that line
                            for genes in mouse:
                                # extracting mouse gene ids
                                genes = genes.lstrip()
                                genes = genes.split('.')
                                # saving bat protein id to each mouse gene id in the line
                                for prot in bat:
                                    # extracting bat protein ids
                                    prot = prot.lstrip()
                                    # populating dictionary
                                    row_dict = {'protid_bat': prot, 'geneid_m': genes[0]}
                                    dict_list.append(row_dict)
                        # many2one mouse-bat             
                        else:
                            for genes in mouse:
                                genes = genes.lstrip()
                                genes = genes.split('.')
                                row_dict = {'protid_bat': line[2], 'geneid_m': genes[0]}
                                dict_list.append(row_dict)
                    # one2many mouse-bat            
                    else:
                        bat = line[2].split(',')
                        for prot in bat:
                                genes = line[1].split('.')
                                prot = prot.lstrip()
                                row_dict = {'protid_bat': prot, 'geneid_m': genes[0]}
                                dict_list.append(row_dict)
            # one2one ortholog (no comma in line)                           
            else:
                # print(line)
                columns = line.split('\t')
                mouse = columns[1].split('.')
                # print(human)
                # print(human[0])
                key = columns[2]
                orthombat1[key] = mouse[0]
    # create data frame from dictionary            
    orthombatmany = pd.DataFrame.from_dict(dict_list)

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
    df_bat = pd.DataFrame(id_prot.items(), columns=['protid_bat','gene_bat'])
    df_mouse = pd.DataFrame(mouse_prot.items(), columns=['geneid_m','gene_m'])
    df_ortho1 = pd.DataFrame(orthombat1.items(), columns=['protid_bat','geneid_m'])
    combi1one = pd.merge(df_bat, df_ortho1, how="right", on="protid_bat")
    combi2one = pd.merge(combi1one, df_mouse, how="left", on="geneid_m")
    combi1many = pd.merge(df_bat, orthombatmany, how="right", on="protid_bat", validate="one_to_many")
    combi2many = pd.merge(combi1many, df_mouse, how="left", on="geneid_m", validate="many_to_one")
    print(combi1many)
    print(df_ortho1)
    print(combi2many)
    combi2one.to_csv('ortho_mbat_one.csv', index=False)
    combi2many.to_csv('ortho_mbat_many.csv', index=False)
    df_bat.to_csv('df_bat.csv', index=False)
    df_mouse.to_csv('df_mouse.csv', index=False)
    orthombatmany.to_csv('orthombatmany.csv', index=False)