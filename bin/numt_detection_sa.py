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
output_prefix = args.output
ref_file = args.reference

df = pd.DataFrame(columns=['nuc_chrom', 'nuc_position1', 'nuc_position2', 'mt_start', 'mt_end'])

breakend_dict = {}

fr = pysam.AlignmentFile(input_file, 'rc', reference_filename=ref_file)
fw = pysam.AlignmentFile(output_prefix + ".numts.SA.bam", 'wb', template=fr)

for read in fr.fetch('chrM'):
    write = False
    if read.has_tag('SA:Z'):
        splits = read.get_tag('SA:Z').split(';')
        splits.remove('')
        for part in splits:
            if part != '':
                chrom = part.split(',')[0]
                if chrom != 'chrM':
                    write = True
        if write:
            fw.write(read)

            cigar_counts = get_cigar_counts(read.cigarstring)
            reference_consuming = cigar_counts['M'] + cigar_counts['D']


            leading_sc = re.match(r'^(\d+)[SH]', read.cigarstring)
            leading_sc_int = -1
            trailing_sc = re.search(r'(\d+)[SH]$', read.cigarstring)
            trailing_sc_int = -1

            if trailing_sc:
                mt_breakend_end = read.reference_start + reference_consuming
                trailing_sc_int = int(trailing_sc[1])
            else:
                trailing_sc_int = 0
                mt_breakend_end = -1
                nuc_end_chrom = -1
                nuc_end_pos = -1

            if leading_sc:
                mt_breakend_start = read.reference_start
                leading_sc_int = int(leading_sc[1])
            else:
                leading_sc_int = 0
                mt_breakend_start = -1
                nuc_start_chrom = -1
                nuc_start_pos = -1

            print(f'trailing sc: {trailing_sc_int}')
            print(f'leading sc: {leading_sc_int}')
            print(read.reference_start)
            print(read.cigarstring)
            print(read.get_tag('SA:Z').split(';'))

            if len(splits) == 1:
                sa_chrom, sa_pos, sa_strand, sa_cigar, *_ = splits[0].split(",")
                sa_cigar_counts = get_cigar_counts(sa_cigar)
                sa_reference_consuming = sa_cigar_counts['M'] + sa_cigar_counts['D']
                
                if leading_sc_int > trailing_sc_int:
                    mt_breakend_end = -1
                    nuc_breakend_1 = int(sa_pos)
                else:
                    mt_breakend_start = -1
                    nuc_breakend_1 = int(sa_pos) + sa_reference_consuming
                    
                new_row_data = [sa_chrom, nuc_breakend_1, -1, mt_breakend_start, mt_breakend_end]
                print(new_row_data)

            elif len(splits) == 2:
                ## first check if both splits to same chrom -- if not skip
                sa_chrom1, *_ = splits[0].split(",")
                sa_chrom2, *_ = splits[0].split(",")

                if sa_chrom1 == sa_chrom2:
                    nuclear_intervals = []

                    for sa in splits:

                        sa_chrom, sa_pos, sa_strand, sa_cigar, *_ = sa.split(",")
                        
                        sa_cigar_counts = get_cigar_counts(sa_cigar)
                        sa_reference_consuming = sa_cigar_counts['M'] + sa_cigar_counts['D']
                        sa_start_pos = int(sa_pos)
                        sa_end_pos = int(sa_pos) + sa_reference_consuming

                        nuclear_intervals.append((sa_chrom, sa_start_pos, sa_end_pos))
                        print((sa_chrom, sa_start_pos, sa_end_pos))

                    nuc_breakend_1  = min([end for _, start, end in nuclear_intervals])
                    nuc_breakend_2 = max([start for _, start, end in nuclear_intervals])

                    new_row_data = [sa_chrom1, nuc_breakend_1, nuc_breakend_2, mt_breakend_start, mt_breakend_end]
                    print(new_row_data)
                else:
                    continue

            elif len(splits) > 2:
                continue

            ### add row to df 
            df.loc[len(df)] = new_row_data 


df_collapsed = df.groupby(df.columns.tolist(), dropna=False).size().reset_index(name='read_count')
df_collapsed = df_collapsed[(df_collapsed['read_count'] >= 4) & (df_collapsed['mt_start'] != 0) & (df_collapsed['mt_end'] != 16569)]
df_collapsed.to_csv(output_prefix + '.numts.SA.tsv', sep='\t', index=False)