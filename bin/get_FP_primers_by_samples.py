#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import csv
import pandas as pd
import utilities
import settings
import os
import yaml
import run_tripod_analysis as rta
from collections import Counter

'''
running scripts below 
set threshold to 38 for isolate HMAS, and 23 for stool HMAS
python3 get_FP_primers_by_samples.py 
    -i /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/isolate_HMAS/MN_M05688_240618-422066253/MN_M05688_240618_v1.2.1_20240906_163532 
    -o test 
    -r /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/tripod/dev/MN_M05688_240322_v1.0.0_20241007_104321/reference.fasta 
    -p ../MN_M05688_240618-422066253_sample_isolates_IS_mapping.csv 
    -t isolate

python3 get_FP_primers_by_samples.py \
    -i /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/stool_HMAS/MN_M05688_240322-413314107/MN_M05688_240322-413314107_v1.2.1_20240908_131907 \
    -o test \
    -r /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/tripod/dev/MN_M05688_240322_v1.0.0_20241007_104321/reference.fasta \
    -p ../MN_M05688_240618-422066253_sample_isolates_mapping.csv \
    -t stool \
    -s /scicomp/groups/OID/NCEZID/DFWED/EDLB/projects/CIMS/HMAS_pilot/stool_HMAS/MN_M05688_240322-413314107/MN_M05688_240322-413314107_v1.2.1_20240908_131907/MN_M05688_240322-413314107_good_samples.csv
'''

def parse_argument():
    parser = argparse.ArgumentParser()
    #HMAS step_mothur pipeline output folder
    #parser.add_argument('-is', '--input2', metavar = '', required = True, help = 'Specify input file folder isolate')
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify input file folder')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify the tripod report output file')
    #predicted amplicon fasta file (coming from amplicon extraction part, usually named as reference.fasta by default)
    parser.add_argument('-r', '--reference', metavar = '', required = True, help = 'reference_fasta file')
    #look-up table between sample name and all possible WGS isolates in the sample
    parser.add_argument('-p', '--mapping', metavar = '', required = True, help = 'sample isolates mapping file')
    #parser.add_argument('-ps', '--isolate_mapping', metavar = '', required = True, help = 'sample stool mapping file')
    # currently, tag can be either 'stool' or 'isolate'
    parser.add_argument('-t', '--tag', metavar = '', required = True, help = 'either stool or isolate')
    parser.add_argument('-s', '--sample_list', metavar = '', required = True, help = 'selected sample names file')

    return parser.parse_args()

# read the sample list file (single column, no header) and convert it to a list
def get_sample_list(csv_file):

    df = pd.read_csv(csv_file, header=None)
    return df[0].to_list()

if __name__ == "__main__":
    
    args = parse_argument()

    oligo_primers = utilities.Primers(settings.OLIGO_FILE)
    sample_isolate_dict = rta.map_sample_to_isolate(args.mapping)
    predict_amplicon_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))

    FP_dict = {} # key: sample name; value: FP primer list
    TP_dict = {}
    for folder_name in rta.get_folders(args.input):
        #skip those isolates without any valid unique sequences
        fasta_file = f'{args.input}/{folder_name}/{folder_name}.final.unique.fasta'
        if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) <= 0:
            continue

        #skip those isolates NOT in the given list (csv file)
        #  our sample_name is probably shorter than 'folder_name'
        if all(sample_name not in folder_name for sample_name in get_sample_list(args.sample_list)):
            continue

        highest_size_ids = rta.extract_highest_size_ids(fasta_file, oligo_primers.pseqs.keys())

        TP_counter = 0
        FP_counter = 0
        diff_counter = 0

        FP_list = [] # list of false positive primers
        TP_list = []

        for primer, seq in highest_size_ids:
            seq_IDs = [] #hold all amplicon predictions for this primer-folder_name
            for isolate in sample_isolate_dict[folder_name]:
                seq_IDs.extend([key for key in predict_amplicon_dict if f'{primer}-{isolate}' in key])
            if seq_IDs: 
                if any(seq == predict_amplicon_dict[seq_ID].seq \
                    or rta.revcomp(seq) == predict_amplicon_dict[seq_ID].seq for seq_ID in seq_IDs):
                    TP_counter += 1
                    TP_list.append(primer)
                    # break
                else:
                    diff_counter += 1
            else:
                FP_counter += 1
                FP_list.append(primer)

        FP_dict[folder_name] = FP_list 
        TP_dict[folder_name] = TP_list

        

    FP_dict = TP_dict

    # Flatten the lists into a single list and count occurrences
    all_elements = [item for sublist in FP_dict.values() for item in sublist]
    element_counts = Counter(all_elements)

    # Filter elements with occurrences larger than threshold
    threshold = 6
    elements_larger_than = [(element, count) for element, count in element_counts.items() if count > threshold]

    # Print the result
    print(f"Total {len(elements_larger_than)} FP primers")
    print(sorted(elements_larger_than, key=lambda x: x[1], reverse=True))
    # print(f"Total {len(elements_larger_than)} FP primers")






