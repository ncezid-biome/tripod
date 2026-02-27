#!/usr/bin/env python3

from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import argparse
import pandas as pd
import utilities
import settings
import os
import yaml


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
    # highest_size_ids = [(pattern,pattern_dict[pattern][1]) for pattern in patterns if pattern_dict[pattern][0] is not None]
    highest_size_ids = [pattern_dict[pattern][1] for pattern in patterns]
    
    return highest_size_ids

def get_folders(folder_name):
    
    # Get all folder names in the input folder
    folder_names = [f for f in os.listdir(folder_name) if os.path.isdir(os.path.join(folder_name, f))]
    
    return folder_names

# read the sample list file (single column, no header) and convert it to a list
def get_sample_list(csv_file):

    df = pd.read_csv(csv_file, header=None)
    return df[0].to_list()

def make_report_yaml(output_file, data_df, yaml_id, yaml_description, col_1, col_2):
    '''
    this method generates a custom content yaml file specific for the multiqc report

    Parameters
    ----------
    output_file: output file (yaml) name
    data_df: data part of the yaml file in the format of dataframe

    Returns: None
    ----------
    '''   
    # Create headers dictionary
    headers = {
        col_1: {
            'title': 'Stool sample stats',
            'description': 'min,median,max of sequence read depth',
            # 'format': '{:,.0f}',
            # "scale": False
        },
        col_2: {
            'title': 'Isolate sample stats',
            'description': 'min,median,max of sequence read depth',
            # 'format': '{:,.0f}',
            # "scale": False
        },
    }

    # Convert the DataFrame to the required format
    data_yaml = data_df.to_dict(orient='index')

    # Create the full YAML dictionary
    yaml_dict = {
        'id': yaml_id,
        'section_name': yaml_id,
        'description': yaml_description,
        'plot_type': 'table',
        'pconfig': {
            'id': yaml_id,
            'sort_rows': False,
            'col1_header': 'Primer Name',
            "no_violin": True,
        },
        'headers': headers,
        'data': data_yaml
    }

    # Write to a YAML file
    with open(output_file, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)


def parse_argument():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--isolate', metavar = '', required = True, help = 'Specify isolate input file folder')
    parser.add_argument('-l', '--stool', metavar = '', required = True, help = 'Specify stool input file folder')
    # for now, we have the outputs named as: tripod_primer_performance_1_mqc.yaml, tripod_primer_performance_2_mqc.yaml by default
    # parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify the output file')
    
    #list of good sample names (usually defined as with mean read depth over 30 AND
    #% of successful primer pairs is over 90% for the stool samples)
    #one for the stool samples and another for the matching isolate samples 
    parser.add_argument('-is', '--isolate_samples', metavar = '', required = True, help = 'selected sample names file')
    parser.add_argument('-ls', '--stool_samples', metavar = '', required = True, help = 'selected sample names file')

    return parser.parse_args()

if __name__ == "__main__":
    
    args = parse_argument()

    oligo_primers = utilities.Primers(settings.OLIGO_FILE)

    def get_primer_performance_df(input_folder, sample_list):

        result_dict = {}
        for folder_name in get_folders(input_folder):
            #skip those isolates without any valid unique sequences
            fasta_file = f'{input_folder}/{folder_name}/{folder_name}.final.unique.fasta'
            if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) <= 0:
                continue

            #skip those isolates NOT in the given list (csv file)
            #  our sample_name is probably shorter than 'folder_name'
            if all(sample_name not in folder_name for sample_name in get_sample_list(sample_list)):
                continue

            highest_size_ids = extract_highest_size_ids(fasta_file, sorted(oligo_primers.pseqs.keys()))

            result_dict[folder_name] = highest_size_ids

        # Convert dictionary to DataFrame
        df = pd.DataFrame(result_dict, index=sorted(oligo_primers.pseqs.keys()))
        return df
    

    df_SH = get_primer_performance_df(args.stool, args.stool_samples)
    df_IH = get_primer_performance_df(args.isolate, args.isolate_samples)

    # Step 1: Indexes where SH's average is smaller than IH's average
    SH_avg = df_SH.mean(axis=1)
    IH_avg = df_IH.mean(axis=1)
    smaller_avg_indexes = df_SH.index[SH_avg < IH_avg]
    smaller_avg_reverse_indexes = df_SH.index[SH_avg > IH_avg]

    value_threshold = 10 #read depth threshold
    ratio_threshold = 0.55 #so the majority of samples has the value over the value_threshold

    def get_avg_as_df(filtered_indexes):

        # Compute a tuple (min, median, max) for each row in df_SH and df_IH based on filtered indexes
        SH_filtered_stats = df_SH.loc[filtered_indexes].apply(
            lambda x: f"min: {x.min()}, median: {x.median()}, max: {x.max()}", axis=1
        )

        IH_filtered_stats = df_IH.loc[filtered_indexes].apply(
            lambda x: f"min: {x.min()}, median: {x.median()}, max: {x.max()}", axis=1
        )

        # Combine them into a new DataFrame with side-by-side values
        result_df = pd.DataFrame({
            'SH_avg': SH_filtered_stats,
            'IH_avg': IH_filtered_stats
        }, index=filtered_indexes)

        return result_df


    # Step 2: Indexes where both SH and IH are below the threshold 50% of the time
    SH_below_threshold = (df_SH < value_threshold).mean(axis=1)
    IH_below_threshold = (df_IH < value_threshold).mean(axis=1)
    below_threshold_indexes = df_SH.index[(SH_below_threshold >= ratio_threshold) & (IH_below_threshold >= ratio_threshold)]

    # Step 3: Indexes where SH is below the threshold 50% of the time while IH is above threshold 100% of the time
    SH_below_threshold = (df_SH < value_threshold).mean(axis=1)
    IH_above_threshold = (df_IH > value_threshold).mean(axis=1)
    above_threshold_indexes = df_SH.index[(SH_below_threshold >= ratio_threshold) & (IH_above_threshold == 1)]

    # Step 4: Indexes where both SH and IH are above threshold 100% of the time
    SH_above_threshold = (df_SH > value_threshold).mean(axis=1)
    IH_above_threshold = (df_IH > value_threshold).mean(axis=1)
    both_above_threshold_indexes = df_SH.index[(SH_above_threshold == 1) & (IH_above_threshold == 1)]
    
    df1 = get_avg_as_df(below_threshold_indexes)
    df1.columns = ['pcol1-1','pcol1-2']
    make_report_yaml('tripod_primer_performance_1_mqc.yaml', df1, 
                     'primer_performance_1',
                     'stats for primers that fail in both stool and isolate HMAS samples',
                     'pcol1-1', 'pcol1-2')
    
    df2 = get_avg_as_df(above_threshold_indexes)
    df2.columns = ['pcol2-1','pcol2-2']
    make_report_yaml('tripod_primer_performance_2_mqc.yaml', df2, 
                     'primer_performance_2',
                     'stats for primers that fail in stool but are perfect in isolate HMAS samples',
                     'pcol2-1', 'pcol2-2')
    
    df3 = get_avg_as_df(both_above_threshold_indexes)
    df3.columns = ['pcol3-1','pcol3-2']
    make_report_yaml('tripod_primer_performance_3_mqc.yaml', df3, 
                     'primer_performance_3',
                     'stats for primers that are both perfect in stool and isolate HMAS samples',
                     'pcol3-1', 'pcol3-2')