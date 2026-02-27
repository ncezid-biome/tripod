#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import pandas as pd
import utilities
import settings
import os
import run_tripod_analysis as rta
from Bio import SeqRecord
from Bio.SeqRecord import SeqRecord
import re

'''
running scripts below 
set threshold to 38 for isolate HMAS, and 23 for stool HMAS
python3 get_FP_sequences_by_samples.py 
    -i /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/isolate_HMAS/MN_M05688_240618-422066253/MN_M05688_240618_v1.2.1_20240906_163532 
    -o test 
    -r /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/tripod/dev/MN_M05688_240322_v1.0.0_20241007_104321/reference.fasta 
    -p ../MN_M05688_240618-422066253_sample_isolates_IS_mapping.csv 
    -t isolate

python3 get_FP_sequences_by_samples.py 
    -i /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/stool_HMAS/MN_M05688_240322-413314107/MN_M05688_240322-413314107_v1.2.1_20240908_131907 
    -o test 
    -r /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/tripod/dev/MN_M05688_240322_v1.0.0_20241007_104321/reference.fasta 
    -p ../MN_M05688_240618-422066253_sample_isolates_mapping.csv 
    -t stool
'''

sample_list = ["CIMS-MN-055", "CIMS-MN-063", "CIMS-MN-051", "CIMS-MN-069", "CIMS-MN-064", "CIMS-MN-067", "CIMS-MN-060", "CIMS-MN-053"]

# Total 110 FP primers (primers in 14 out of 14 'good' STOOL samples)
FP_primer_list = [
    "OG0002993primerGroup0",
    "OG0000406primerGroup1",
    "OG0002066primerGroup7",
    "OG0002615primerGroup4",
    "OG0003163primerGroup9",
    "OG0002423primerGroup4",
    "OG0002293primerGroup9",
    "OG0002892primerGroup0",
    "OG0001427primerGroup5",
    "OG0001618primerGroup0",
    "OG0002038primerGroup9",
    "OG0003002primerGroup7",
    "OG0002026primerGroup8",
    "OG0003239primerGroup0",
    "OG0002044primerGroup3",
    "OG0001451primerGroup1",
    "OG0003306primerGroup7",
    "OG0001930primerGroup5",
    "OG0001068primerGroup0",
    "OG0002455primerGroup9",
    "OG0001760primerGroup0",
    "OG0003061primerGroup0",
    "OG0002499primerGroup1",
    "OG0002151primerGroup2",
    "OG0001826primerGroup1",
    "OG0003241primerGroup6",
    "OG0001799primerGroup3",
    "OG0002158primerGroup7",
    "OG0001969primerGroup0",
    "OG0000952primerGroup9",
    "OG0001609primerGroup4",
    "OG0003044primerGroup0",
    "OG0000402primerGroup1",
    "OG0002165primerGroup7",
    "OG0002264primerGroup7",
    "OG0002693primerGroup6",
    "OG0000371primerGroup4",
    "OG0002692primerGroup6",
    "OG0002655primerGroup0",
    "OG0001755primerGroup8",
    "OG0002607primerGroup5",
    "OG0001502primerGroup4",
    "OG0002632primerGroup0",
    "OG0002287primerGroup1",
    "OG0001467primerGroup3",
    "OG0002102primerGroup5",
    "OG0002069primerGroup0",
    "OG0002073primerGroup4",
    "OG0003291primerGroup0",
    "OG0002901primerGroup2",
    "OG0001679primerGroup1",
    "OG0001567primerGroup3",
    "OG0001534primerGroup4",
    "OG0002265primerGroup9",
    "OG0002005primerGroup9",
    "OG0003206primerGroup0",
    "OG0001508primerGroup5",
    "OG0003034primerGroup6",
    "OG0000372primerGroup9",
    "OG0003607primerGroup7",
    "OG0001867primerGroup5",
    "OG0001085primerGroup3",
    "OG0001069primerGroup2",
    "OG0002242primerGroup2",
    "OG0002017primerGroup3",
    "OG0003508primerGroup0",
    "OG0002137primerGroup3",
    "OG0002487primerGroup5",
    "OG0001471primerGroup2",
    "OG0002354primerGroup3",
    "OG0002900primerGroup2",
    "OG0002384primerGroup2",
    "OG0002244primerGroup0",
    "OG0002430primerGroup5",
    "OG0001646primerGroup1",
    "OG0003281primerGroup0",
    "OG0001980primerGroup0",
    "OG0000385primerGroup0",
    "OG0003531primerGroup7",
    "OG0002403primerGroup2",
    "OG0001726primerGroup2",
    "OG0002513primerGroup0",
    "OG0002274primerGroup3",
    "OG0001148primerGroup6",
    "OG0002476primerGroup4",
    "OG0001388primerGroup8",
    "OG0002545primerGroup2",
    "OG0001124primerGroup6",
    "OG0002055primerGroup0",
    "OG0001027primerGroup9",
    "OG0003735primerGroup0",
    "OG0000513primerGroup8",
    "OG0002459primerGroup1",
    "OG0001458primerGroup4",
    "OG0002169primerGroup2",
    "OG0003486primerGroup1",
    "OG0002039primerGroup5",
    "OG0003529primerGroup9",
    "OG0001465primerGroup3",
    "OG0003056primerGroup1",
    "OG0002903primerGroup5",
    "OG0002174primerGroup5",
    "OG0003540primerGroup4",
    "OG0001776primerGroup6",
    "OG0002539primerGroup2",
    "OG0000553primerGroup0",
    "OG0003560primerGroup3",
    "OG0002500primerGroup6",
    "OG0000294primerGroup8",
    "OG0002394primerGroup9"
]

def extract_highest_size_ids(fasta_file, patterns):
    # Dictionary to store the highest size for each pattern
    # we store: seq_ID, abundance_value, seq
    pattern_dict = {pattern: (None, 0, None) for pattern in patterns}

    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Extract the sequence ID and size
        seq_id = record.id
        size = int(seq_id.split(";size=")[-1])
        
        # Check each pattern
        for pattern in patterns:
            if pattern in seq_id:
                # Update if the current size is higher
                if size > pattern_dict[pattern][1]:
                    pattern_dict[pattern] = (seq_id, size, record.seq)
                # If sizes are the same, keep the first one
                # elif size == pattern_dict[pattern][1] and pattern_dict[pattern][0] is None:
                #     pattern_dict[pattern] = (seq_id, size)
    
    # Extract the list of highest size IDs for each pattern
    highest_size_ids = [(pattern_dict[pattern][0],pattern_dict[pattern][2]) for pattern in patterns if pattern_dict[pattern][0] is not None]
    
    return highest_size_ids

def parse_argument():
    parser = argparse.ArgumentParser()
    #HMAS step_mothur pipeline output folder
    parser.add_argument('-ii', '--input_isolate', metavar = '', required = True, help = 'Specify input file folder isolate')
    parser.add_argument('-is', '--input_stool', metavar = '', required = True, help = 'Specify input file folder stool')
    parser.add_argument('-oi', '--output_isolate', metavar = '', required = True, help = 'Specify output file folder isolate')
    parser.add_argument('-os', '--output_stool', metavar = '', required = True, help = 'Specify output file folder stool')

    return parser.parse_args()

# read the sample list file (single column, no header) and convert it to a list
def get_sample_list(csv_file):

    df = pd.read_csv(csv_file, header=None)
    return df[0].to_list()


def make_FP_sequences(in_folder, out_file):

    FP_seq_list = [] # list of false positive primers
    for folder_name in rta.get_folders(in_folder):
        #skip those isolates without any valid unique sequences
        fasta_file = f'{in_folder}/{folder_name}/{folder_name}.final.unique.fasta'
        if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) <= 0:
            continue

        #skip those isolates NOT in the given list (csv file)
        #  our sample_name is probably shorter than 'folder_name'
        if all(sample_name not in folder_name for sample_name in sample_list):
            continue

        highest_size_ids = extract_highest_size_ids(fasta_file, FP_primer_list)

        for seq_id, seq in highest_size_ids:

            #adjust the length of seq_id to make it work with 50 character limit with blast
            # !! this only works with the type of seq_id below
            #>M05688:169:000000000-LGJKB:1:1101:11728:13719=OG0002892primerGroup0=CIMS-MN-064-SH-K-N-MN-M05688-240322_S16;size=7
            # Step 1: Remove 'CIMS-MN-'
            seq_id = seq_id.replace("CIMS-MN-", "")   
            # Step 2: Remove 'M05688-' and everything until the next underscore
            seq_id = re.sub(r'M05688-[^_]*_', '', seq_id)
            # Step 3: Keep only the last 50 characters
            if len(seq_id) > 50:
                seq_id = seq_id[-50:]

            ampliconRec = SeqRecord(seq, id=seq_id)
            FP_seq_list.append(ampliconRec)

    with open(f'{os.getcwd()}/{out_file}.fasta', 'w') as f:
        seqs_list = [f">{rec.id}\n{rec.seq}" for rec in FP_seq_list]
        f.write('\n'.join(seqs_list) + '\n')

if __name__ == "__main__":
    
    args = parse_argument()

    make_FP_sequences(args.input_stool, args.output_stool)
    make_FP_sequences(args.input_isolate, args.output_isolate)

    












