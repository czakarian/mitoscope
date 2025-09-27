#!/usr/bin/env python3

"""This script removes reads in an aligned bam file that are representative of: 
(1) fold-backs (split read with alignments on both + and - strands),
(2) NUMTs (reads with significant levels of unaligned softclipping)."""

import argparse
import pysam
import re
import os

def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
    parser.add_argument("-i", "--input", help="File path of input bam.", required=True)
    parser.add_argument("-m", "--max_sc_threshold", type=int, help="Maximum unaligned softclipping allowed in a read.", required=False, default=100)
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
    
    #print(counts)
    return counts

def is_foldback(read_parts):
    """
    Check if a read resembles a foldback with alignments on both forward/reverse strands.
    Args:
        read_parts (list): List of lines from bam (pysam) with same read name.
    Returns:
        bool: True if alignments on both + and - strand, otherwise False.
    """
    strands = {line.is_reverse for line in read_parts}
    return len(strands) > 1

def is_NUMT(read_parts, max_unaligned_threshold):
    """
    Determine whether a read resembles a NUMT sequence based on level of unaligned soft-clipping across alignments.
    Args:
        read_parts (list): List of lines from bam (pysam) with same read name.
        max_unaligned_threshold (int): The maximum amount of unaligned soft-clipping allowed (Default = 1000).
    Returns:
        bool: True if resembles a NUMT sequence, otherwise False.
    """

    total_query_lengths = []
    mapped_query_length = 0

    for line in read_parts:
        cigar_counts = get_cigar_counts(line.cigarstring)
        total_query_lengths.append(cigar_counts['M'] + cigar_counts['I'] + cigar_counts['S'] + cigar_counts['H'])
        mapped_query_length += cigar_counts['M'] + cigar_counts['I']

    if len(set(total_query_lengths)) > 1:
        print('Warning: Query lengths differ between alignments with same read name.')

    unaligned_bases = total_query_lengths[0] - mapped_query_length
    #print(unaligned_bases)
    return unaligned_bases > max_unaligned_threshold

def process_read_group(read_parts, read_counts, keep_file_handle, discard_file_handle):
    """
    Process a group of reads and categorize as foldback, NUMT, or kept read.
    Args:
        read_parts (list): List of lines from bam (pysam) with same read name.
    Returns:
        None
    """
    if is_foldback(read_parts):
        category, target = 'foldback', discard_file_handle
    elif is_NUMT(read_parts, args.max_sc_threshold):
        category, target = 'numt', discard_file_handle
    else:
        category, target = 'kept', keep_file_handle

    read_counts[category] += 1
    for line in read_parts:
        target.write(line)
    #print(line.query_name)

# Define output file names
input_bam = args.input
input_bam_name_sorted = args.input[:-4] + ".sortedbyname.bam"
filtered_bam = args.input[:-4] + ".filtered.bam"
filtered_bam_name_sorted = args.input[:-4] + ".filtered.sortedbyname.bam"
discarded_bam = args.input[:-4] + ".discardReads.bam"
discarded_bam_name_sorted = args.input[:-4] + ".discardReads.sortedbyname.bam"

# Perform sorting by name
pysam.sort("-n", "-o", input_bam_name_sorted, input_bam)

fr = pysam.AlignmentFile(input_bam_name_sorted, 'rb')
fw = pysam.AlignmentFile(filtered_bam_name_sorted, 'wb', template= fr)
fd = pysam.AlignmentFile(discarded_bam_name_sorted, 'wb', template= fr)

read_counts = {"foldback": 0, "numt": 0, "kept": 0}
read_parts = []
current_readname = None

for line in fr:
    if current_readname == line.query_name:
        read_parts.append(line)
    else:
        if read_parts:  # Process the previous group
            process_read_group(read_parts, read_counts, fw, fd)
        read_parts = [line]
        current_readname = line.query_name

# Process the final group
if read_parts:
    process_read_group(read_parts, read_counts, fw, fd)

## close read/write files
fr.close()
fw.close()
fd.close() 

# Perform sorting by name
pysam.sort("-o", filtered_bam, filtered_bam_name_sorted)
pysam.sort("-o", discarded_bam, discarded_bam_name_sorted)

os.remove(input_bam_name_sorted)
os.remove(filtered_bam_name_sorted)
os.remove(discarded_bam_name_sorted)

print("Kept reads: ", read_counts["kept"])
print("Discarded foldback reads: ", read_counts["foldback"])
print("Discarded NUMT reads: ", read_counts["numt"])