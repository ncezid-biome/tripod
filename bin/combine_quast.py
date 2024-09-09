#!/usr/bin/env python3

import pandas as pd
import argparse
import yaml

def make_quast_yaml(tsv_files, yaml_output):
    '''
    this method generates a custom content yaml file specific for the multiqc report
    this yaml file is for the quast output

    Parameters
    ----------
    tsv_files: list, a list of tsv files (quast output)
    yaml_output: String, the generated yaml file name

    Returns: None
    ----------
    '''   
    def create_df(tsv_file):

        data = {}
        # Open and read the target file
        with open(tsv_file, 'r') as file:
            for line in file:
                # Remove whitespace characters like `\n` at the end of each line
                line = line.strip()
                
                # Check for the required lines and extract the values
                if line.startswith("Assembly"):
                    assembly_name = line.split('\t')[1]
                    data['Assembly'] = assembly_name
                elif line.startswith("N50"):
                    data['N50'] = int(line.split('\t')[1])*0.001
                elif line.startswith("# contigs (>= 1000 bp)"):
                    data['# contigs (>= 1000 bp)'] = int(line.split('\t')[1])
                elif line.startswith("Largest contig"):
                    data['Largest contig'] = int(line.split('\t')[1])*0.001
                elif line.startswith("Total length (>= 0 bp)"):
                    data['Total length (>= 0 bp)'] = int(line.split('\t')[1])*0.001

                # Create a DataFrame using the extracted values
                column_order = ['Assembly','N50','# contigs (>= 1000 bp)','Largest contig','Total length (>= 0 bp)']
                df = pd.DataFrame([data], columns=column_order)
                df.set_index('Assembly', inplace=True)

        df.columns = ['q_col1', 'q_col2', 'q_col3', 'q_col4']

        return df

    # Create headers dictionary
    headers = {
        'q_col1': {
            'title': 'N50 (Kbp)',
            'description': 'N50 is the contig length such that using longer or equal length contigs produces half (50%) of the bases of the assembly.',
            'format': '{:,.1f}',
        },
        'q_col2': {
            'title': '# of contigs >= 1Kbp',
            'description': 'The number of contigs with size over 1Kbp',
            'format': '{:,.0f}',
            "scale": False
        },
        'q_col3': {
            'title': 'Largest contig (Kbp)',
            'description': 'The size of the largest contig of the assembly',
            'format': '{:,.1f}',
        },
        'q_col4': {
            'title': 'Length (Kbp)',
            'description': 'The total number of bases in the assembly.',
            'format': '{:,.1f}',
            "scale": False
        },
    }

    # Convert the DataFrame to the required format
    data_yaml = pd.concat([create_df(tsv_file) for tsv_file in tsv_files]).to_dict(orient='index')

    # Create the full YAML dictionary
    yaml_dict = {
        'id': 'quast_report',
        'section_name': 'quast_report',
        'description': 'Assembly Statistics using Quast',
        'plot_type': 'table',
        'pconfig': {
            'id': 'quast_report',
            'sort_rows': False,
        },
        'headers': headers,
        'data': data_yaml
    }

    # Write to a YAML file
    with open(yaml_output, 'w') as file:
        yaml.dump(yaml_dict, file, sort_keys=False)

def parse_argument():

    parser = argparse.ArgumentParser(prog = 'combine_quast.py')
    parser.add_argument('-i', '--input', metavar = '', required = True, help = 'Specify input file')
    parser.add_argument('-o', '--output', metavar = '', required = True, help = 'Specify output file')

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_argument()
    
    #read each report csv file as a df
    file_list = [file for file in args.input.split()]
        
    # report_df = pd.concat(df_list)


    # # now generate report yaml file for MultiQC report
    # # change column name
    # report_df.columns = ['col1','col2','col3']
    # #update empty cell to n/a
    # report_df.fillna('n/a', inplace=True)
    make_quast_yaml(file_list, args.output)
