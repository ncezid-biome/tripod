#!/usr/bin/env python3

from Bio import SeqIO
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

PCENT_PRIMER_SUCCESS = 0.9

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
            'title': '# of ident seqs to isolates',
            'description': '(matched) stool isolate',
            'format': '{:,.0f}',
            "scale": False
        },
        f'col4_{tag_name}': {
            'title': '% of ident seqs to isolates',
            'description': '% of matched most abund. seqs to isolate',
            'format': '{:,.3f}',
            "scale": False
        },
        f'col5_{tag_name}': {
            'title': '# of ident seqs to insilico',
            'description': '(matched) insilico',
            'format': '{:,.0f}',
            "scale": False
        },
        f'col6_{tag_name}': {
            'title': '% of ident seqs to insilico',
            'description': '% of matched most abund. seqs to insilico',
            'format': '{:,.3f}',
            "scale": False
        },
        f'col7_{tag_name}': {
            'title': 'Mean read depth',
            'description': 'we include only reads with at least 2 sequence count',
            'format': '{:,.1f}',
            "scale": False
        },
        f'col8_{tag_name}': {
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
            "**note, if the columns of ident seqs to insilico are zero, it's likely we don't have the matching WGS**",
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


def make_primer_report_yaml(output_file, data_df, yaml_id, yaml_description, col_1, col_2):
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

    if data_df.empty:
        data_yaml = {
            "No primers passed filter": {
                "pcol2-1": "-",
                "pcol2-2": "-"
            }
        }
    else:
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

def extract_all_seqs(fasta_file, patterns):
    """
    Extract all sequences (record.seq) for each pattern.

    Returns:
        A list of tuples in the form:
            [(pattern, [seq1, seq2, ...]), ...]
    """
    # Dictionary: pattern -> list of sequences
    pattern_dict = {pattern: [] for pattern in patterns}

    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id

        # Check each pattern
        for pattern in patterns:
            if pattern in seq_id:
                pattern_dict[pattern].append(record.seq)

    # Convert dictionary to list of tuples (like original output)
    all_seqs = [
        (pattern, pattern_dict[pattern])
        for pattern in patterns
        if pattern_dict[pattern]  # only include non-empty
    ]

    return all_seqs


def extract_highest_size_ids(fasta_file, patterns, return_type="seqs"):
    """
    Extract the highest-abundance sequence (by size) for each pattern.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file.
    patterns : list
        List of string patterns to match in sequence IDs.
    return_type : str, optional
        Determines the output format. Options:
        - "seqs": [(pattern, sequence), ...]
        - "sizes": [size, ...]
    
    Returns
    -------
    list
        Depending on return_type, returns either sequences with their patterns
        or just the size values.
    """
    
    # Dictionary to store the highest size for each pattern
    pattern_dict = {pattern: (None, 0, None) for pattern in patterns}

    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        size = int(seq_id.split(";size=")[-1])

        for pattern in patterns:
            if pattern in seq_id:
                if size > pattern_dict[pattern][1]:
                    pattern_dict[pattern] = (seq_id, size, record.seq)

    # Return format depending on requested output
    if return_type == "seqs":
        return [(pattern, pattern_dict[pattern][2]) 
                for pattern in patterns if pattern_dict[pattern][0] is not None]
    elif return_type == "sizes":
        # Always return one value per pattern, fill with 0 if missing
        return [pattern_dict[pattern][1] if pattern_dict[pattern][0] is not None else 0
                for pattern in patterns]
    else:
        raise ValueError("Invalid return_type. Use 'seqs' or 'sizes'.")


def revcomp(myseq):
    rc = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G', 'U' : 'A', 'Y' : 'R', 'R' : 'Y', 'K':'M', 'M':'K','B':'V',\
            'D':'H', 'H':'D', 'V':'B', 'N':'N'}
    seq = [rc[n] if n in rc else n for n in myseq] # allow non-IUPAC code to stay as is
    return("".join(list(reversed(seq))))


def get_primer_performance_df(input_stool, input_isolate, good_stools, tripod_mapping, oligo_primers):

    result_stool_dict = {}
    result_isolate_dict = {}
    for stool in good_stools:

        pattern = os.path.join(input_stool, "*", f"{stool}*", f"{stool}*.final.unique.fasta")
        matches = glob.glob(pattern)
        if len(matches) <= 0:
            continue
        fasta_file = matches[0]
        #skip those without any valid unique sequences
        if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) <= 0:
            continue
        highest_size_ids = extract_highest_size_ids(fasta_file, sorted(oligo_primers.pseqs.keys()), 'sizes')
        result_stool_dict[stool] = highest_size_ids
        
        ## now get the isolate fasta file
        pattern_isolate = os.path.join(input_isolate, "*", f"{tripod_mapping[stool][0]}*", f"{tripod_mapping[stool][0]}*.final.unique.fasta")
        matches = glob.glob(pattern_isolate)
        if len(matches) <= 0:
            continue
        fasta_file_isolate = matches[0]
        if os.path.exists(fasta_file_isolate) and os.path.getsize(fasta_file_isolate) > 0:
            highest_size_ids = extract_highest_size_ids(fasta_file_isolate, sorted(oligo_primers.pseqs.keys()), 'sizes')
            result_isolate_dict[tripod_mapping[stool][0]] = highest_size_ids

    df_stool = pd.DataFrame(result_stool_dict, index=sorted(oligo_primers.pseqs.keys()))
    df_isolate = pd.DataFrame(result_isolate_dict, index=sorted(oligo_primers.pseqs.keys()))
        
    return df_stool, df_isolate


def parse_argument():
    parser = argparse.ArgumentParser()
    #HMAS step_mothur pipeline output folder
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify stool input file folder')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify the tripod report output file')
    #predicted amplicon fasta file (coming from amplicon extraction part, usually named as reference.fasta by default)
    parser.add_argument('-r', '--reference', metavar = '', required = True, help = 'reference_fasta file')
    # currently, tag can be either 'stool' or 'isolate'
    parser.add_argument('-t', '--tag', metavar = '', required = True, help = 'either stool or isolate') # this will be 'combined'
    parser.add_argument('-s', '--isolates', metavar = '', required = True, help = 'Specify isolates input file folder')
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
    df = pd.read_csv(args.correlate)#, nrows=10)
    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    # Build dictionary: key = 2nd column (stool hmas), value = (1st col -- isolate hmas, 3rd col -- wgs)
    tripod_mapping = df.set_index(df.columns[1])[[df.columns[0], df.columns[2]]].apply(tuple, axis=1).to_dict()
    
    result_dict = {}
    FP_dict = {} # key: sample name; value: FP primer list
    reports_stool = [] # hold the hmas run report for each sample
    reports_isolate = [] # it's not being used for now
    for stool in tripod_mapping:

        pattern = os.path.join(args.input, "**", f"{stool}*.final.unique.fasta")
        matches = glob.glob(pattern, recursive=True)

        if len(matches) <= 0:
            continue
        fasta_file = matches[0]
        #skip those samples without any valid unique sequences
        if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) <= 0:
            continue

        # Replace extension ".final.unique.fasta" → ".csv"
        csv_file = fasta_file.replace(".final.unique.fasta", ".csv")
        reports_stool.append(csv_file)
        
        # highest_size_ids = extract_all_seqs(fasta_file, oligo_primers.pseqs.keys())
        highest_size_ids = extract_highest_size_ids(fasta_file, oligo_primers.pseqs.keys())

        ## now get the isolate fasta file
        pattern_isolate = os.path.join(args.isolates, "*", f"{tripod_mapping[stool][0]}*", f"{tripod_mapping[stool][0]}*.final.unique.fasta")
        matches = glob.glob(pattern_isolate)
        if len(matches) <= 0:
            continue
        fasta_file_isolate = matches[0]
        #skip those samples without any valid unique sequences
        if not os.path.exists(fasta_file_isolate) or os.path.getsize(fasta_file_isolate) <= 0:
            continue
        
        # Replace extension ".final.unique.fasta" → ".csv"
        csv_file_iso = fasta_file_isolate.replace(".final.unique.fasta", ".csv")
        reports_isolate.append(csv_file_iso)
        
        highest_size_ids_isolate = []
        highest_size_ids_isolate = extract_highest_size_ids(fasta_file_isolate, oligo_primers.pseqs.keys())

        TP_counter = 0
        FP_counter = 0
        diff_counter = 0
        match_to_isolate_counter = 0

        def matches_any(seq_list, target_seq):
            return any(s == target_seq or revcomp(s) == target_seq for s in seq_list)

        FP_list = [] # list of false positive primers
        # for primer, seqs in highest_size_ids: #seqs is a list
        for primer, seq in highest_size_ids:

            seq_IDs = [key for key in predict_amplicon_dict if f'{primer}-{tripod_mapping[stool][1]}' in key]

            if seq_IDs:
                if any(seq == predict_amplicon_dict[seq_ID].seq or revcomp(seq) == predict_amplicon_dict[seq_ID].seq for seq_ID in seq_IDs):
                    TP_counter += 1
                else:
                    diff_counter += 1
            else:
                FP_counter += 1
                FP_list.append(primer)
            
            if len(highest_size_ids_isolate) > 0:
                isolate_seq = [x for x in highest_size_ids_isolate if x[0] == primer]
                if len(isolate_seq) > 0:
                    isolate_seq = isolate_seq[0][1]
                
                if isolate_seq == seq or revcomp(isolate_seq) == seq:
                    match_to_isolate_counter += 1

        result_dict[stool] = (
            TP_counter,
            len(highest_size_ids),
            len(highest_size_ids_isolate) if len(highest_size_ids_isolate) > 0 else None,
            FP_counter,
            diff_counter,
            match_to_isolate_counter
        )

        FP_dict[stool] = FP_list # not currently using it 

    # Convert dictionary to DataFrame
    df = pd.DataFrame.from_dict(result_dict, orient='index', columns=['item1', 'item2', 'item3', '# of FP', '# of not ident seqs', '# of ident seqs to isolates'])

    # Create a new column
    df['% of ident seqs to isolates'] = (df['# of ident seqs to isolates'] / df['item2'] * 100).round(1).astype(str) + '%'

    df['# of ident seqs to isolates'] = df.apply(
        lambda row: f" ({row['# of ident seqs to isolates']})  {row['item2']} / {row['item3']:.0f}" 
                    if row['item3'] is not None else None,
        axis=1
    )

    # 2nd new column
    df['# of matched to insilico'] = df['item1'].astype(str)

    # 3rd new column
    df['% of ident seqs to insilico'] = (df['item1'] / df['item2'] * 100).round(1).astype(str) + '%'

    # Drop the original 'item1' and 'item2' columns
    df.drop(columns=['item1', 'item2', 'item3'], inplace=True)

    ### add 2 new columns from HMAS2 report to tripod report
    # 1. concat our hmas report.csv files for stool
    if reports_stool:   # only proceed if list is non-empty
        df_hmas = pd.concat([pd.read_csv(f, index_col=0) for f in reports_stool])
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
    df_hmas = df_hmas[~df_hmas.index.duplicated(keep='first')]

    # 2: Add 2 new columns from DataFrame df_hmas to current df
    
    if not df_hmas.empty and df_hmas.shape[1] >= 2:

        stool_depth = df_hmas.iloc[:, 0]
        stool_succ  = df_hmas.iloc[:, 1]

        df['mean read depth'] = df.index.map(stool_depth)
        df['% successful primers'] = df.index.map(stool_succ)

    else:
        print("Warning: df_hmas is empty or has insufficient columns.")
    
        
    #################################################
    #################################################
    ### add 2 new columns from HMAS2 report.csv files for isolates
    if reports_isolate:
        df_hmas_iso = pd.concat([pd.read_csv(f, index_col=0) for f in reports_isolate])
    else:
        df_hmas_iso = pd.DataFrame()

    # Map long isolate report IDs to short isolate IDs (tripod_mapping[stool][0])
    isolate_ids = [tripod_mapping[s][0] for s in df.index if s in tripod_mapping]

    index_mapping_iso = {
        long_id: short_id
        for short_id in isolate_ids
        for long_id in df_hmas_iso.index
        if long_id.startswith(short_id)
    }
    df_hmas_iso = df_hmas_iso.rename(index=index_mapping_iso)
    # Keep the first occurrence for each isolate id (or use keep='last')
    df_hmas_iso = df_hmas_iso[~df_hmas_iso.index.duplicated(keep='first')]

    # Add isolate columns
    df['mean read depth (isolate)'] = None
    df['% successful primers (isolate)'] = None

    # Create a helper column mapping stool->isolate so we can join cleanly
    df['__isolate_id'] = df.index.map(lambda s: tripod_mapping[s][0] if s in tripod_mapping else None)
  
    if not df_hmas_iso.empty and df_hmas_iso.shape[1] >= 2:

        iso_depth = df_hmas_iso.iloc[:, 0]
        iso_succ  = df_hmas_iso.iloc[:, 1]

        df['mean read depth (isolate)'] = df['__isolate_id'].map(iso_depth)
        df['% successful primers (isolate)'] = df['__isolate_id'].map(iso_succ)

    else:
        print("Warning: df_hmas_iso is empty or has insufficient columns.")
    
    
    ############### end of newly added block #################
    
    df.to_csv(final_filename)
    
    ########################################################
    ########################################################
    # Require both stool and isolate to pass
    stool_pass = (df['% successful primers'] > PCENT_PRIMER_SUCCESS)

    # isolate column might be None for some rows; treat as fail
    iso_pass = (df['% successful primers (isolate)'].fillna(0) > PCENT_PRIMER_SUCCESS)
    df_filtered = df[stool_pass & iso_pass].copy()
    good_stools = list(df_filtered.index)
    
    # All stools from the mapping (from the very beginning)
    all_mapped_stools = sorted(list(tripod_mapping.keys()))

    # Stools that made it into the main results dataframe (were processed)
    processed_stools = set(df.index)

    # Stools that passed final filter
    passed_stools = set(good_stools)

    # Anything mapped but not passed is "filtered out" (including never-processed)
    filtered_out_stools = [s for s in all_mapped_stools if s not in passed_stools]

    filtered_out_file = f"{os.path.splitext(base_filename)[0]}_filtered_out_stools.txt"
    with open(filtered_out_file, "w") as f:
        for s in filtered_out_stools:

            if s not in processed_stools:
                # These were skipped earlier (no fasta, empty fasta, no isolate fasta, etc.)
                f.write(f"{s}\tnot being processed, either the stool or isolates samples is missing in the input folder\n")
            else:
                sp = bool(stool_pass.loc[s]) if s in stool_pass.index else False
                ip = bool(iso_pass.loc[s]) if s in iso_pass.index else False
                f.write(f"{s}\tstool_pass={sp}\tisolate_pass={ip}\n")
    
    ############### end of newly added block #################
    
    avg_row = df_filtered.mean(numeric_only=True).round(2)
    avg_row.name = "Average"
    df_result = pd.concat([df_filtered, avg_row.to_frame().T])

    # If you don't want __isolate_id in outputs:
    df.drop(columns=['__isolate_id','% successful primers (isolate)','mean read depth (isolate)'], inplace=True, errors='ignore')
    df_result.drop(columns=['__isolate_id','% successful primers (isolate)','mean read depth (isolate)'], inplace=True, errors='ignore')
    
    final_filename_new = f"{os.path.splitext(base_filename)[0]}_goodsamples{extension}"
    df_result.to_csv(final_filename_new)    
    
    
    
    # df.columns = ['col1','col2','col3', 'col4']
    df.columns = [f'col1_{args.tag}',f'col2_{args.tag}',f'col3_{args.tag}', f'col4_{args.tag}', f'col5_{args.tag}', f'col6_{args.tag}', f'col7_{args.tag}', f'col8_{args.tag}']
    make_report_yaml(f'tripod_{args.tag}_mqc.yaml', df, args.tag)


    ###################################
    # primer performance part
    df_SH, df_IH = get_primer_performance_df(args.input, args.isolates, good_stools, tripod_mapping, oligo_primers)

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


    # Step 1: Indexes where both SH and IH are below the threshold 50% of the time
    SH_below_threshold = (df_SH < value_threshold).mean(axis=1)
    IH_below_threshold = (df_IH < value_threshold).mean(axis=1)
    below_threshold_indexes = df_SH.index[(SH_below_threshold >= ratio_threshold) & (IH_below_threshold >= ratio_threshold)]

    # Step 2: Indexes where SH is below the threshold 50% of the time while IH is above threshold 99.9% of the time
    SH_below_threshold = (df_SH < value_threshold).mean(axis=1)
    IH_above_threshold = (df_IH >= value_threshold).mean(axis=1)
    above_threshold_indexes = df_SH.index[(SH_below_threshold >= ratio_threshold) & (IH_above_threshold  >= 0.999)]

    # Step 3: Indexes where both SH and IH are above threshold 99.9% of the time
    SH_above_threshold = (df_SH >= value_threshold).mean(axis=1)
    IH_above_threshold = (df_IH >= value_threshold).mean(axis=1)
    both_above_threshold_indexes = df_SH.index[(SH_above_threshold  >= 0.999) & (IH_above_threshold  >= 0.999)]
    
    df1 = get_avg_as_df(below_threshold_indexes)
    df1.columns = ['pcol1-1','pcol1-2']
    make_primer_report_yaml('tripod_primer_performance_1_mqc.yaml', df1, 
                     'primer_performance_1',
                     'stats for primers that fail in both stool and isolate HMAS samples',
                     'pcol1-1', 'pcol1-2')
    
    df2 = get_avg_as_df(above_threshold_indexes)
    df2.columns = ['pcol2-1','pcol2-2']
    make_primer_report_yaml('tripod_primer_performance_2_mqc.yaml', df2, 
                     'primer_performance_2',
                     'stats for primers that fail in stool but are perfect in isolate HMAS samples',
                     'pcol2-1', 'pcol2-2')
    
    df3 = get_avg_as_df(both_above_threshold_indexes)
    df3.columns = ['pcol3-1','pcol3-2']
    make_primer_report_yaml('tripod_primer_performance_3_mqc.yaml', df3, 
                     'primer_performance_3',
                     'stats for primers that are both perfect in stool and isolate HMAS samples',
                     'pcol3-1', 'pcol3-2')
