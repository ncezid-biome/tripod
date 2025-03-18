#!/usr/bin/env python
'''
    this script is to replace Sean's validationModule.py
    https://github.com/ncezid-biome/T3Pio/blob/Main/T3Pio_Main/validationModule.py
    
    The script will run primersearch on the given isolate sequence file and a list of primers, then
    parse the output to extract amplicon sequences. The script will generate 4 files: 
    1. output from primersearch
    2. extracted amplicon sequence
    3. a list of not-found primers
    4. a csv file of metasheet table with 3 columns(seqid, primer, isolate)
    
    Note: EMBOSS/6.4.0 is recommended.  As EMBOSS/6.6.0 seems to report un-amplified primers for 
    those primers with the first base being ambiguous code (either on forward or reverse primer) 
    It does work fine with the ambiguous bases at any other places though

    Sean's script, as of 10/14/2022, has the issue that it always assumes the forward primer hits
    the forward sequence (same goes to the reverse primer).  But this is not always the case with 
    PrimerSearch result.
    
    https://emboss.sourceforge.net/apps/cvs/emboss/apps/primersearch.html

    The primer list for primersearch need to have the format of:
    primer_name  forward_primer reverse_primer
'''

'''Extracts PCR amplicons from full-length nucleotide fasta sequences, using Emboss primersearch output. 
Reads in full-length fasta sequences and corresponding Emboss primersearch output for the sequences. 
Parses primersearch output for sequence id and description and uses primer locations to extract amplicon
regions from full-length fasta sequences. If more than one amplimer is present in the primersearch output, 
sequence is extracted using the first amplimer. Prints extracted amplicon sequences to fasta file.
Note: Fasta definition lines should follow the typical format: 
	e.g. >KY925925 A/Santo Antonio da Patrulha/LACENRS-2621/2012 2012/08/10 7 (MP)
Instead of:
	e.g. >A_/_H1N1_|_A/WAKAYAMA/163/2016_|_MP|_|_970779
 
    This script is adapted from Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory
'''

'''usage: python3 extract_amplicon_from_primersearch_output.py 
        -s fastaToParse.fasta -p primers_list_file
    add --metasheet , if you want to create the metasheet table
'''

import re
from Bio import Seq, SeqIO, SeqUtils, SeqRecord
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import argparse
import subprocess
import shutil
import sys
import pandas as pd
import os
import glob
import concurrent.futures
from pathlib import Path
import logging
from logging.handlers import RotatingFileHandler

amplicon_file_extension = '_extractedAmplicons.fasta'
non_match_primer_file_extension = '_not_match_primers.txt'
metasheet_file_extension = '_metasheet.csv'
metasheet_table_columns = ['seq_id', 'primer', 'sample']
mismatch_percent = 6
max_amplicon_len = 375
max_seqid_len = 80
max_worker_num = 10
LOG_FILE = r'primersearch_amplicon_prediction.log'

log_formatter = logging.Formatter('%(asctime)s %(levelname)s %(filename)s(%(lineno)d) - %(message)s')
log_handler = RotatingFileHandler(LOG_FILE, mode='a', maxBytes=5*1024*1024,
                                 backupCount=5, encoding=None, delay=0)
log_handler.setFormatter(log_formatter)
log_handler.setLevel(logging.INFO)

logger = logging.getLogger('root')
logger.setLevel(logging.INFO)
logger.addHandler(log_handler) 

def extractIndex(hit_line_string):
    '''Extract and return hit position from Emboss primersearch forward or reverse strand hit line.'''
    regex = re.compile(r"strand\s+at\s+([0-9]{1,})\s+with")
    clean_str = hit_line_string.replace('[', '').replace(']', '')
    matchArray = regex.findall(clean_str)
    return int(matchArray[0])

def extractPrimerLength(hit_line_string):
    '''Extract and return primer length from Emboss primersearch forward or reverse strand hit line.'''
    regex = re.compile(r"^(\w+)\s+hits")
    clean_str = re.sub(r'\[\w+\]','N',hit_line_string.strip())
    matchArray = regex.findall(clean_str)
    return (len(matchArray[0]))

def extractAmpliconLength(amplimer_length_line):
    '''Extract and return amplicon length from Emboss primersearch amplimer length line.'''
    amplicon_length = ''
    regex = re.compile("[0-9]{1,}")
    matchArray = regex.findall(amplimer_length_line)
    if len(matchArray) > 0:  # check for a match
        amplicon_length = int(matchArray[0])
    return amplicon_length

def parsePrimerSearch(primersearch_results, full_length_dict, file_base):
    # parse primersearch results file
    with open(primersearch_results, 'r') as inputFile:
        extracted_amplicon_list = []  # empty list of extracted amplicon sequence
        not_match_primer_list = []
        
        # these 3 lists are for creating a meta-sheet table
        seqid_list = []
        seqid_primer_list = []
        seqid_isolate_list = []

        for line in inputFile:  # read in each line
            if "Amplimer " in line.strip():  # grab number of amplier
                ampl_count += 1
                line = inputFile.readline().strip() # read the next line
                if "Sequence: " in line:  # grab sequence identifier
                    seq_id = line.replace("Sequence: ", "").strip()
                    description = inputFile.readline().strip()  # grab description from next line

                    forward_hit_line = inputFile.readline().strip()
                    forward_primer_length = extractPrimerLength(forward_hit_line)

                    reverse_hit_line = inputFile.readline().strip()  # next line will be reverse hit
                    reverse_primer_length = extractPrimerLength(reverse_hit_line)

                    if "forward" in forward_hit_line: # check if F primer hits the forward strand
                        forwardHitPosition = extractIndex(forward_hit_line)  # determine F or R hit positions
                    else:
                        forwardHitPosition = extractIndex(reverse_hit_line) # now R primer hits the forward strand
                        logger.info(f'For Amplimer {ampl_count}, FR primer of {primer_name} hit reverse strand ! {seq_id}')

                    amplimer_length_line = inputFile.readline().strip()  # grab next line
                    amplicon_length = extractAmpliconLength(amplimer_length_line)  # determine length of amplicon
                    startIndex = forwardHitPosition - 1  # primersearch results weren't zero-indexed
                    endIndex = startIndex + amplicon_length

                    # exclude primer sequences
                    startIndex = startIndex + forward_primer_length
                    endIndex = endIndex - reverse_primer_length

                    # check if extracted amplicon is over the max length
                    if endIndex - startIndex <= max_amplicon_len:
                        
                        full_length_record = full_length_dict[seq_id]  # created SeqRecord from amplicon and add to list
                        extracted_sequence = full_length_record.seq[startIndex:endIndex]
                        ampliconRec = SeqRecord(extracted_sequence, id=f'{primer_name}-{file_base}-ampl{ampl_count}'[:max_seqid_len])
                        extracted_amplicon_list.append(ampliconRec)
                        
                        seqid_list.append(f'{primer_name}-{file_base}-ampl{ampl_count}'[:max_seqid_len])
                        seqid_primer_list.append(primer_name)
                        seqid_isolate_list.append(file_base)
                        
                        # make sure it's not empty list and we only do this once
                        if not_match_primer_list and ampl_count == 1:
                            not_match_primer_list.pop() #remove this matched primers from the list

                    else:
                        logger.info(f'For {seq_id}, \nAmplimer {ampl_count}, {primer_name} is too long with length {endIndex - startIndex}')

            elif "Primer name " in line:
                # [\w-] is for any word or hyphen (which is not included in word definition)
                # word is any alphanumerics plus underscore
                result = re.search(r'Primer name\s+([\w-]{1,})', line.strip())
                primer_name = result.group(1)
                not_match_primer_list.append(primer_name)
                ampl_count = 0
        
    return (extracted_amplicon_list, not_match_primer_list, 
            [seqid_list, seqid_primer_list, seqid_isolate_list])

def parse_argument():
    # note
    # usage: python3 extract_amplicon_from_primersearch_output.py 
    #       -s fastaToParse.fasta -p primers_list_file
    # 02/10/2023 -s now can also be a directory which holds all the assemblies files
    parser = argparse.ArgumentParser(prog = 'extract_amplicon_from_primersearch_output.py')
    parser.add_argument('-p', '--primersearch', metavar = '', required = True, help = 'Specify primersearch output')
    parser.add_argument('-s', '--sequence', metavar = '', required = True, help = 'Specify isolate sequence file')
    parser.add_argument("--metasheet", action="store_true", help="Enable creating metasheet table")
    return parser.parse_args()


if __name__ == "__main__":
    
    args = parse_argument()
        
    output_dir = f"{Path.cwd()}/primersearch/"
    os.makedirs(output_dir, exist_ok=True)
    file_base = Path(args.sequence).stem # basename of the fasta file without .fasta extension
        
    # read full-length fasta sequences into dict of SeqrRecords with key = record.id
    full_length_dict = SeqIO.to_dict(SeqIO.parse(args.sequence, "fasta"))
    extracted_amplicon_list, not_match_primer_list, _3lists_for_df \
        = parsePrimerSearch(args.primersearch, full_length_dict, file_base)
    
    # not using SeqIO.write because it automatically wraps sequence 60bp per line
    #SeqIO.write(extracted_amplicon_list, ampliconOutputHandle, "fasta")
    with open(f'{output_dir}{file_base}{amplicon_file_extension}', 'w') as f:
        seqs_list = [f">{rec.id}\n{rec.seq}" for rec in extracted_amplicon_list]
        f.write('\n'.join(seqs_list) + '\n')
    
    with open(f'{output_dir}{file_base}{non_match_primer_file_extension}', 'w') as f:
        f.write('\n'.join(not_match_primer_list) + '\n')
    
    # generate the metasheet table 
    # turn off metasheet table by default, metasheet table is only used with confusion_matrix
    if args.metasheet:
        df = pd.DataFrame(_3lists_for_df).T
        df.columns = metasheet_table_columns
        df.to_csv(f"{output_dir}{file_base}{metasheet_file_extension}", index=False)