#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 14:33:29 2025

@author: bradenkeiser
"""

import os, sys
import pandas as pd
import numpy as np
from Bio import SeqIO
import subprocess

'''
This script makes a db of all phages using both genomad and phastest phage regions.
An all-against-all BLASTn search is then performed to begin matching phages. Sequences
are taken from the database by querying the blast_hits. Various transformations are 
used to isolate equivalent non-equivalent phages identified across the programs. 
    These include: 
    1. whether or not the query and target are the same strain and contig
        (phastest and genomad received the same contig information at input - 
         the contig numberings will be the same)
    2. whether or not the lengths of the query and target match represent at least
        95% of the length of the smallest phage in the comparison
    3. Whether each query-target pair is an even-odd (phastest-genomad) pair
        --> we dont want to count phastest-phastest pairs until later
    
'''



def make_blast_db(input_fasta, db_type, output_db):
    """
    Creates a BLAST database from a FASTA file.

    Parameters:
    - input_fasta (str): Path to the input FASTA file
    - db_type (str): Type of sequences ("prot" for proteins, "nucl" for nucleotides)
    - output_db (str): Name of the output database

    Returns:
    - None
    """

    # Ensure the input file exists
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"FASTA file '{input_fasta}' not found.")

    # Command to run makeblastdb
    cmd = [
        "makeblastdb",
        "-in", input_fasta,
        "-dbtype", db_type,
        "-out", output_db
    ]

    # Run the command
    try:
        subprocess.run(cmd, check=True)
        print(f"BLAST database '{output_db}' created successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
    except FileNotFoundError:
        print("Error: makeblastdb not found. Ensure BLAST+ is installed and in your PATH.")


blast_type = 'n'

dir_in = '' # ideally, where dir_genomad and dir_phastest are located

dir_genomad = '' # preferably within the input directory
dir_phastest  = '' #preferably within the input directory

file_genomad = 'genomad_phageome.fasta' # full phage genomes, 1 seq per genome
file_phastest = 'phastest_phageome.fasta' # full phage genomes, 1 seq per genome

path_all_genomad = os.path.join(dir_in, dir_genomad, file_genomad)

path_all_phastest = os.path.join(dir_in, dir_phastest, file_phastest)

phages_genomad = SeqIO.parse(path_all_genomad, 'fasta')

phages_phastest= SeqIO.parse(path_all_phastest, 'fasta')
'''PART 1: MAKE A DATABASE OF ALL THE PHAGES '''

# Merge the two FASTA files into a single one for BLAST database creation
merged_fasta = os.path.join(dir_in, "merged_phageome.fasta")

with open(merged_fasta, "w") as outfile:
    for fasta_file in [path_all_phastest, path_all_genomad]:
        with open(fasta_file, "r") as infile:
            outfile.write(infile.read())

# Create the BLAST database
output_db = os.path.join(dir_in, "phastest_genomad_db/phas-genomad_phage_db")
make_blast_db(merged_fasta, 'nucl', output_db)


'''PART 2: BLAST THE PHASTEST PHAGES AGAINST THE TOTAL DB'''

'''We want to blast all against all'''

blast_output = os.path.join(dir_in, 'phgc_phages_blasthits.blast.out')
if blast_type == 'n':
    print('Working on nucleotide BLAST')
    blast_command = [
        "blastn",
        "-query", merged_fasta,
        "-db", output_db,
        "-evalue", "0.000001",
        "-outfmt", "6",
        "-out", blast_output
    ]
if not os.path.exists(blast_output):
    try:
        subprocess.run(blast_command, check = True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running BLAST: {e}")
        sys.exit(1)
else:
    print('blast already performed, output file present, skipping...')
    
    
phas_phages = SeqIO.parse(path_all_phastest, 'fasta')
gnmd_phages = SeqIO.parse(path_all_genomad, 'fasta')
phas_len_dict = {}
phage_seq_dict = {} # get a dictionary of program_genome : sequence
for phage in phas_phages:
    key = phage.name
    seq = phage.seq
    description=phage.description.split(';')[0]
    length = description.split('=')[1]
    phas_len_dict[key] = length
    phage_seq_dict[key] = seq
gnmd_len_dict = {}
for phage in gnmd_phages:
    key = phage.name
    seq = phage.seq
    description=phage.description.split(';')[0]
    length = description.split('=')[1]
    gnmd_len_dict[key] = length
    phage_seq_dict[key] = seq

phage_len_dict = gnmd_len_dict | phas_len_dict

print('the total number of genomad phages here is: {}'.format(len(gnmd_len_dict)))
print('the total number of phastest phages here is: {}'.format(len(phas_len_dict)))


columns = ["query", "blast_hit", "seqid", "length", "mismatch", "gapopen", 
           "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

blast_output_df = pd.read_table(blast_output, header = None)
blast_output_df.columns = columns

#get the length of the query
blast_output_df['query_len'] = blast_output_df['query'].map(phage_len_dict).astype(int)
blast_output_df['hit_len'] = blast_output_df['blast_hit'].map(phage_len_dict).astype(int)

blast_output_df['same_strain'] = blast_output_df.apply(
    lambda row: 1 if '_'.join(row['query'].split('_')[1:5]) == '_'.join(
        row['blast_hit'].split('_')[1:5]) else 0, axis=1)

blast_output_df['phage_program'] = blast_output_df.apply(
    lambda row: 'both' if row['query'].split('_')[0] != row['blast_hit'].split('_')[0]
    else ('phastest_only' if row['query'].split('_')[0] == 'PH' else 'genomad_only'), axis = 1
    )
blast_output_df['matched_length_plausible'] = blast_output_df.apply(
    lambda row: 1 if row['length'] >= 0.95*min(row['query_len'], row['hit_len']) 
    else 0, axis =1)
blast_output_df['longer_phage'] = blast_output_df.apply(lambda row: row['query'] if row['query_len'] > row['hit_len'] else row['blast_hit'], axis = 1)
# it does not matter if the query or the blast_hit are selected; it matters later if that's gc or ph
blast_output_df['longer_seq'] = blast_output_df['longer_phage'].map(phage_seq_dict)
#Now, start filtering the data by these new, comparative results


'''PART THREE: START TRANSFORMING BY THE NEW FEATURES'''

#ensure the matched lengths of sequences are within 10% genome size to each other
blast_minlen = blast_output_df[blast_output_df['matched_length_plausible'] == 1]


blast_samestrain = blast_minlen[blast_minlen['same_strain'] == 1] # this necessitates that we use paired strain_contig viruses

blast_hitlength = blast_samestrain[
    (blast_samestrain['query_len'] >= blast_samestrain['hit_len'] * 0.80) | 
    (blast_samestrain['hit_len'] >= blast_samestrain['query_len'] * 0.80)
    ]

#drop duplicates of the paired columns
blast_hits_paired = blast_hitlength.drop_duplicates(['query', 'blast_hit'])
blast_hits_paired = blast_hits_paired.drop_duplicates(['blast_hit', 'query'])

blast_hits_paired['ph-gc'] = blast_hits_paired.apply(
    lambda row: 1 if (row['query'].split('_')[0] == 'PH' and row['blast_hit'].split('_')[0] == 'GC') or
                      (row['query'].split('_')[0] == 'GC' and row['blast_hit'].split('_')[0] == 'PH') 
                      else 0, axis=1
)

blast_hits_both = blast_hits_paired[blast_hits_paired['ph-gc'] == 1]
# Step 1: Sort pairs to ensure (A, B) and (B, A) are stored as (A, B)
unique_pairs = set()  # Initialize a set to track unique pairs
both_dict = {}  # Initialize the dictionary to store key-value pairs

# Iterate through query and blast_hit columns
for q, h in zip(blast_hits_both['query'], blast_hits_both['blast_hit']):
    key = tuple(([q, h]))  # Sort and store as tuple to ensure consistency
    # Only add the pair if it's not already in the set
    if key not in unique_pairs: #and antikey not in unique_pairs:
        unique_pairs.add(key)  # Add the pair to the set
        both_dict[key[0]] = key[1]  # Store only one direction in the dictionary
        
#use gc since we know it is larger
gc_list = blast_hits_paired[blast_hits_paired['query'].str.startswith('GC')]
ph_list = blast_hits_paired[blast_hits_paired['blast_hit'].str.startswith('PH')]


genomad_frame = pd.DataFrame({'genomad':gc_list['query'], 
                              'phastest':gc_list['query'].map(both_dict)})
phastest_frame = pd.DataFrame({'phastest':ph_list['query'], 
                               'genomad':ph_list['query'].map(both_dict)}) 

# isolate only genomad results, only phastest results
only_genomad = genomad_frame[genomad_frame['phastest'].isna() == True]
only_phastest = phastest_frame[phastest_frame['genomad'].isna() == True]

#isolate those phages found across programs
both_programs_genomad = genomad_frame[genomad_frame['phastest'].isna() == False]
both_programs_phastest = phastest_frame[phastest_frame['genomad'].isna() == False]

both_programs = pd.concat([both_programs_genomad], ignore_index=True)
both_programs = both_programs.drop_duplicates('genomad')

longer_phage_dict_prep = blast_hits_paired.drop_duplicates('query', ignore_index = True)
longer_phage_dict_prep = blast_hits_paired.drop_duplicates('blast_hit', ignore_index = True)

blast_hits_paired_longer_phage_deduplicate = longer_phage_dict_prep.drop_duplicates(
    'longer_phage', 
    ignore_index = True
    )

longer_phage_dict = blast_hits_paired_longer_phage_deduplicate.set_index('query')['longer_phage'].to_dict()


'''PART FOUR: COMBINE THE LISTS '''

blast_all = pd.concat([both_programs, only_genomad, only_phastest], ignore_index = True)
blast_all['longer_phage'] = blast_all['phastest'].map(longer_phage_dict); blast_all.head()
blast_all['longer_phage'] = blast_all['longer_phage'].fillna(0)
blast_all['longer_phage'] = blast_all.apply(lambda row: row['genomad'] if row['longer_phage'] == 0
                                            else row['longer_phage'], 
                                            axis = 1)

merged_df = pd.merge(blast_all, blast_hits_paired_longer_phage_deduplicate, 
                     on='longer_phage', how = 'left'); len(merged_df)


final_csv = merged_df[['longer_phage', 'longer_seq', 'phage_program']]
final_csv['chosen_program_phage'] = final_csv['longer_phage'].apply(
    lambda x: 'phastest' if x.split('_')[0] == 'PH' else 'genomad'
    )
final_csv['seq_length'] = final_csv['longer_seq'].apply(lambda x: len(x))
new_columns = ['longer_phage', 'sequence', 'phage_program', 'chosen_program_phage', 'seq_length']
final_csv.columns = new_columns
final_csv = final_csv[['longer_phage', 'seq_length','phage_program', 'chosen_program_phage', 'sequence']]
final_csv['longer_phage'] = final_csv['longer_phage'].apply(lambda x: '_'.join(x.split('_')[1:]))

#prepare the output files accordingly
output_file_prefix = '' # choose an appropriate name to represent the fasta and csv
output_fasta = os.path.join(dir_in, f'{output_file_prefix}.fasta')
output_excel= os.path.join(dir_in, f'{output_file_prefix}.csv')

#open the final fasta output file to write filtered data lists
with open(output_fasta,'w') as output_handle:
    seq_records = []
    for index, record in final_csv.iterrows():
        sequence = record['sequence']
        name = record['longer_phage']
        phage_program = record['phage_program'] if record['phage_program'] != 'both' else 0
        header  = '_'.join(list(name.split('_')) + [phage_program]) if (phage_program) != 0 else name
        seq_record = SeqIO.SeqRecord(sequence, id=header, 
                                     description=f"seq_len={len(sequence)};program={record['chosen_program_phage']}"
                                     )
        seq_records.append(seq_record)
    SeqIO.write(seq_records, output_fasta, 'fasta')

#write the excel-compatible csv file for easy viewing by all
final_csv.to_csv(output_excel, sep = ',', header = 0)
