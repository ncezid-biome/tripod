#!/usr/bin/env python3

from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import argparse
import csv
import pandas as pd
import utilities
import settings
import os
import datetime
import yaml
import time
import glob

def make_report_yaml(output_file, data_df, tag_name):
    '''
    this method generates a custom content yaml file specific for the multiqc report
    this yaml file is for the final combined hmas summary report

    Parameters
    ----------
    output_file: output file (yaml) name
    data_df: data part of the yaml file in the format of dataframe
    tag_name: a tag to differentiate stool and isolate reports

    Returns: None
    ----------
    '''   
    # Create headers dictionary
    headers = {
        f'col1_{tag_name}': {
            'title': '# of FP',
            'description': 'number of False Positives',
            'format': '{:,.0f}',
            "scale": False
        },
        f'col2_{tag_name}': {
            'title': '# of not ident seqs',
            'description': '# of not ident seqs',
            'format': '{:,.0f}',
            "scale": False
        },
        f'col3_{tag_name}': {
            'title': '# of matched most abund. seqs',
            'description': '# of matched most abund. seqs',
            'format': '{:,.0f}',
            "scale": False
        },
        f'col4_{tag_name}': {
            'title': '% of matched most abund. seqs',
            'description': '% of matched most abund. seqs',
            'format': '{:,.3f}',
            "scale": False
        },
        f'col5_{tag_name}': {
            'title': 'Mean read depth',
            'description': 'we include only reads with at least 2 sequence count',
            'format': '{:,.1f}',
            "scale": False
        },
        f'col6_{tag_name}': {
            'title': '% of successful primer-pairs',
            'description': 'primer pairs with at lease 2 amplicons in the sample',
            'format': '{:,.3f}',
            "scale": False
        },
    }

    # Convert the DataFrame to the required format
    data_yaml = data_df.to_dict(orient='index')

    # Create the full YAML dictionary
    yaml_dict = {
        'id': f'tripod_analysis_{tag_name}',
        'section_name': f'tripod report {tag_name}',
        'description': "Combined summary statistics for all the samples in the run, "
            "primer pairs (out of total 2461 in the Salmonella HMAS primer panel).<br>"
            "**note, if not ident seqs and the numerator of next column are 0, it's likely we don't have the matching WGS**",
        'plot_type': 'table',
        'pconfig': {
            'id': f'tripod_analysis_{tag_name}',
            'sort_rows': False
        },
        'headers': headers,
        'data': data_yaml
    }

    # Write to a YAML file
    with open(output_file, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)


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
    highest_size_ids = [(pattern,pattern_dict[pattern][2]) for pattern in patterns if pattern_dict[pattern][0] is not None]
    
    return highest_size_ids

def get_folders(parent_folder):
    
    # Get all folder names in the input folder
    folder_names = [f for f in os.listdir(parent_folder) if os.path.isdir(os.path.join(parent_folder, f))]
    
    return folder_names

def revcomp(myseq):
    rc = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G', 'U' : 'A', 'Y' : 'R', 'R' : 'Y', 'K':'M', 'M':'K','B':'V',\
            'D':'H', 'H':'D', 'V':'B', 'N':'N'}
    seq = [rc[n] if n in rc else n for n in myseq] # allow non-IUPAC code to stay as is
    return("".join(list(reversed(seq))))

def parse_argument():
    parser = argparse.ArgumentParser()
    #HMAS step_mothur pipeline output folder
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify input file folder')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify the tripod report output file')
    #predicted amplicon fasta file (coming from amplicon extraction part, usually named as reference.fasta by default)
    parser.add_argument('-r', '--reference', metavar = '', required = True, help = 'reference_fasta file')
    # currently, tag can be either 'stool' or 'isolate'
    parser.add_argument('-t', '--tag', metavar = '', required = True, help = 'either stool or isolate')
    parser.add_argument('-c', '--correlate', metavar = '', required = True, help = 'Specify the 3-column look-up(stool/isolate/wgs) table')

    return parser.parse_args()

if __name__ == "__main__":
    
    args = parse_argument()

    timestamp = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    base_filename = args.output
    extension = os.path.splitext(base_filename)[1]
    final_filename = f"{os.path.splitext(base_filename)[0]}_{timestamp}{extension}"

    oligo_primers = utilities.Primers(settings.OLIGO_FILE)
    predict_amplicon_dict = SeqIO.to_dict(SeqIO.parse(args.reference, "fasta"))
    
    # Load the 3-column tripod look-up table and strip spaces from all entries
    df = pd.read_csv(args.correlate)#, nrows=15)

    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # Build dictionary: key = 2nd column (stool hmas), value = (1st col -- isolate hmas, 3rd col -- wgs)
    if args.tag == 'stool':
        tripod_mapping = df.set_index(df.columns[1])[df.columns[2]].to_dict()
    else: # it's isolate
        tripod_mapping = df.set_index(df.columns[0])[df.columns[2]].to_dict()

    
    result_dict = {}
    FP_dict = {} # key: sample name; value: FP primer list
    reports_hmas = [] # hold the hmas run report for each sample

    for stool in tripod_mapping:

        # pattern = os.path.join(args.input, "*", f"{stool}*", f"{stool}*.final.unique.fasta")
        # matches = glob.glob(pattern)
        pattern = os.path.join(args.input, "**", f"{stool}*.final.unique.fasta")
        matches = glob.glob(pattern, recursive=True)
        
        if len(matches) <= 0:
            continue
        fasta_file = matches[0]
        #skip those isolates without any valid unique sequences
        if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) <= 0:
            continue

        # Replace extension ".final.unique.fasta" → ".csv"
        csv_file = fasta_file.replace(".final.unique.fasta", ".csv")
        reports_hmas.append(csv_file)

        highest_size_ids = extract_highest_size_ids(fasta_file, oligo_primers.pseqs.keys())

        TP_counter = 0
        FP_counter = 0
        diff_counter = 0

        FP_list = [] # list of false positive primers
        for primer, seq in highest_size_ids:

            seq_IDs = [key for key in predict_amplicon_dict if f'{primer}-{tripod_mapping[stool]}' in key]

            if seq_IDs:
                if any(seq == predict_amplicon_dict[seq_ID].seq \
                    or revcomp(seq) == predict_amplicon_dict[seq_ID].seq for seq_ID in seq_IDs):
                    TP_counter += 1
                else:
                    diff_counter += 1
            else:
                FP_counter += 1
                FP_list.append(primer)

        result_dict[stool] = (TP_counter, len(highest_size_ids), FP_counter, diff_counter)
        FP_dict[stool] = FP_list # not currently using it 

    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(result_dict, orient='index', columns=['item1', 'item2', '# of FP', '# of not ident seqs'])

    # Create the first column with 'item1 / item2' as a string
    df[f'# of matched most abund. seqs'] = ' ' + df['item1'].astype(str) + ' / ' + df['item2'].astype(str)

    # Create the second column with the division result, rounded to 2 decimal places
    df[f'% of matched most abund. seqs'] = (df['item1'] / df['item2']).round(3)

    # Drop the original 'item1' and 'item2' columns
    # df = df[[f'# of matched most abund. seqs', f'% of matched most abund. seqs']]
    df.drop(columns=['item1', 'item2'], inplace=True)

    if reports_hmas:   # only proceed if list is non-empty
        df_hmas = pd.concat([pd.read_csv(f, index_col=0) for f in reports_hmas])
    else:
        df_hmas = pd.DataFrame()  # create empty DataFrame if no files

    # necessary because df_index usually is in the shorter form
    index_mapping = {
        long_id: short_id
        for short_id in df.index
        for long_id in df_hmas.index
        if long_id.startswith(short_id)
    }

    df_hmas = df_hmas.rename(index=index_mapping)

    # 2: Add 2 new columns from DataFrame df_hmas to current df
    # Ensure both DataFrames share the same index type
    common_index = df.index.intersection(df_hmas.index)

    # Ensure the two columns exist, initialize with None
    df['mean read depth'] = None
    df['% successful primers'] = None

    if not df_hmas.empty and df_hmas.shape[1] >= 2 and len(common_index) > 0:
        df.loc[common_index, 'mean read depth'] = df_hmas.loc[common_index].iloc[:, 0]
        df.loc[common_index, '% successful primers'] = df_hmas.loc[common_index].iloc[:, 1]
    else:
        print("Warning: df_hmas is empty or has insufficient columns / no shared indexes.")

    df.to_csv(final_filename)
    
    df.columns = [f'col1_{args.tag}',f'col2_{args.tag}',f'col3_{args.tag}', f'col4_{args.tag}', f'col5_{args.tag}', f'col6_{args.tag}']
    make_report_yaml(f'tripod_{args.tag}_mqc.yaml', df, args.tag)


