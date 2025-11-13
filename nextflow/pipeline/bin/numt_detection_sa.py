#!/usr/bin/env python3

import argparse
import pysam
import re
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
    parser.add_argument("-i", "--input", help="File path of input bam.", required=True)
    parser.add_argument("-o", "--output", help="Optional output prefix", required=True)
    parser.add_argument("-r", "--reference", help="Optional reference for cram files", required=False)
    return parser.parse_args()

args = get_args()

def get_cigar_counts(cigar):
    """
    Count all operations in a CIGAR string.
    Args:
        cigar (str): A CIGAR string (e.g., "10M5I3D").
    Returns:
        dict: A dictionary with counts of each operation.
    """
    counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0}
    
    # Regular expression to find all operations and their counts
    matches = re.findall(r'(\d+)([MIDNSH])', cigar)
    
    for count, operation in matches:
        counts[operation] += int(count)
    
    return counts

input_file = args.input
output_file = args.output
ref_file = args.reference

df = pd.DataFrame(columns=['mt_start', 'mt_end', 'nuc_start_chrom', 'nuc_start_pos', 'nuc_end_chrom', 'nuc_end_pos'])

breakend_dict = {}

fr = pysam.AlignmentFile(input_file, 'rc', reference_filename=ref_file)
fw = pysam.AlignmentFile(output_file, 'wb', template=fr)

for read in fr.fetch('chrM'):
    write = False
    if read.has_tag('SA:Z'):
        splits = read.get_tag('SA:Z').split(';')
        for part in splits:
            if part != '':
                chrom = part.split(',')[0]
                if chrom != 'chrM':
                    write = True
        if write:
            fw.write(read)

            cigar_counts = get_cigar_counts(read.cigarstring)
            reference_consuming = cigar_counts['M'] + cigar_counts['D']

            leading_sc = re.match(r'^(\d+)S', read.cigarstring)
            trailing_sc = re.search(r'(\d+)S$', read.cigarstring)

            if trailing_sc:
                mt_breakend_end = read.reference_start + reference_consuming
            else:
                mt_breakend_end = -1

            if leading_sc:
                mt_breakend_start = read.reference_start
            else:
                mt_breakend_start = -1

            print(read.reference_start)
            print(read.cigarstring)
            print(read.get_tag('SA:Z').split(';'))

            new_row_data = [mt_breakend_start, mt_breakend_end, -1, -1, -1, -1 ]
            df.loc[len(df)] = new_row_data 
                
            if (mt_breakend_start, mt_breakend_end) not in breakend_dict:
                breakend_dict[(mt_breakend_start, mt_breakend_end)] = 1
            else:
                breakend_dict[(mt_breakend_start, mt_breakend_end)] += 1

            # if soft clipping at only 1 end (1 breakend)
            # identify where the one softclipped portion maps
            # if softclipping at beginning vs end?
            #splits = read.get_tag('SA:Z').split(';')
            # if soft clipping at 2 ends (1 breakend)


for be in breakend_dict:
    if breakend_dict[be] > 2 and be[0] != 0 and be[1] != 16569:
        print(f'{be}\t{breakend_dict[be]}')


df_collapsed = df.groupby(df.columns.tolist(), dropna=False).size().reset_index(name='count')
df_collapsed = df_collapsed[df_collapsed['count'] > 1]
df_collapsed.to_csv('df.csv', index=False)