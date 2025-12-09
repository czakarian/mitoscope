#!/usr/bin/env python3

"""This script removes reads in an aligned bam file that are representative of: 
(1) fold-backs (split read with alignments on both + and - strands),
(2) NUMTs (reads with significant levels of unaligned softclipping)."""

import argparse
import pysam
import seaborn as sns
import re
import os
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
    parser.add_argument("-i", "--input", help="File path of input bam.", required=True)
    parser.add_argument("-o", "--output", help="Optional output prefix", required=False, default="")
    parser.add_argument("-s", "--max_sc_threshold", type=int, help="Maximum unaligned softclipping allowed in a read.", required=False, default=200)
    parser.add_argument("-m", "--max_meth_threshold", type=float, help="Maximum proportion of read that can be methylated, otherwise filtered as NUMT.", required=False, default=0.5)
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

def methylated_read_portion(line, meth_prob_likelihood, min_probability=0.5, site="C-m"):
    """
    Calculate the proportion of a read's CpG sites that are methylated (currently just 5mC)
    Args:
        line (pysam record): One bam record line 
        meth_prob_likelihood (list) : List of meth probability scores. 
        min_probability (int): Minumum probability to determine whether CpG site is methylated
        site (str): Name of the base modification (currently only suport 'C-m')
    Returns:
        int: Proportion of CpG sites that are methylated in read
    """

    base, code = site.split("-")

    probs = []
    if (base, 1, code) in line.modified_bases:
        probs = [x[1]/255 for x in line.modified_bases[('C', 1, 'm')]]
    elif (base, 0, code) in line.modified_bases:
        probs = [x[1]/255 for x in line.modified_bases[('C', 0, 'm')]]
    
    meth_prob_likelihood.extend(probs)

    if len(probs) > 0:
        percent_read_meth = len([site for site in probs if site > min_probability])/len(probs)
    else:
        percent_read_meth = 0
    
    return percent_read_meth

def methylation_plot(meth_per_read_list, meth_prob_likelihood):
    """
    Plot read methylation distribution to check what proportion of sites on across reads are methylated.
    Args:
        meth_per_read_list (list): List of read methylation levels.
        meth_prob_likelihood (list): List of meth probability scores.
    Returns:
        None
    """

    plt.figure(figsize=(6, 4))
    sns.histplot(meth_per_read_list, bins=50, kde=True, color='blue', edgecolor=None)
    plt.axvline(x=args.max_meth_threshold, color='red', linestyle='--', linewidth=1)
    plt.xlabel('Fraction of methylated CpG sites (5mC)')
    plt.ylabel('# of Reads')
    plt.xlim(0,1)
    #plt.grid(True)
    plt.savefig(prefix + '.methylation_per_read.png', dpi=300)
    plt.close()

    plt.figure(figsize=(6, 4))
    sns.histplot(meth_prob_likelihood, bins=50, kde=True, color='blue', edgecolor=None)
    plt.xlabel('Methylation Likelihood (5mC)')
    plt.ylabel('# of CpG sites')
    plt.tight_layout()
    #plt.grid(True)
    plt.savefig(prefix + '.methylation_likelihood.png', dpi=300)
    plt.close()

def plot_ref_consuming_lengths(ref_cons_lengths):
    """
    Plot 
    """

    plt.figure(figsize=(6, 4))
    sns.histplot(ref_cons_lengths, bins=50, kde=True, color='blue', edgecolor=None)
    plt.axvline(x=16569, color='red', linestyle='--', linewidth=1)
    plt.xlabel('Reference-consuming length')
    plt.ylabel('# of Reads')
    plt.tight_layout()
    plt.savefig(prefix + '.ref_consuming_hist_kde.png', dpi=300)
    plt.close()


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

def is_NUMT(read_parts, max_unaligned_threshold, max_methylation_threshold, meth_per_read_list, meth_prob_likelihood, ref_cons_lengths):
    """
    Determine whether a read resembles a NUMT sequence based on level of unaligned soft-clipping across alignments as well as methylation level of the read.
    Args:
        read_parts (list): List of lines from bam (pysam) with same read name.
        max_unaligned_threshold (int): The maximum amount of unaligned soft-clipping allowed
        max_methylation_threshold (int): The maximum threshold a read can be methylated, otherwise considered a NUMT.
        meth_per_read_list (list): List to append methylation level to (used to plot distribution across bam).
        meth_prob_likelihood (list): List of meth probability scores.
    Returns:
        bool: True if resembles a NUMT sequence, otherwise False.
    """

    ## 1. check unaligned softclipping levels
    total_query_lengths = []
    mapped_query_length = 0
    reference_consuming = 0
    
    for line in read_parts:
        cigar_counts = get_cigar_counts(line.cigarstring)
        total_query_lengths.append(cigar_counts['M'] + cigar_counts['I'] + cigar_counts['S'] + cigar_counts['H'])
        mapped_query_length += cigar_counts['M'] + cigar_counts['I']
        reference_consuming += cigar_counts['M'] + cigar_counts['D']

    ref_cons_lengths.append(reference_consuming)

    if len(set(total_query_lengths)) > 1:
        print('Warning: Query lengths differ between alignments with same read name.')

    unaligned_bases = total_query_lengths[0] - mapped_query_length

    ## 2. check methylation level for read
    meth_level = methylated_read_portion(read_parts[0], meth_prob_likelihood)
    meth_per_read_list.append(meth_level)

    if (unaligned_bases > max_unaligned_threshold) or (meth_level > max_methylation_threshold):
        return True
    else:
        return False

def process_read_group(read_parts, read_counts, meth_per_read_list, keep_file_handle, discard_file_handle, meth_prob_likelihood, ref_cons_lengths):
    """
    Process a group of reads and categorize as foldback, NUMT, or kept read.
    Args:
        read_parts (list): List of lines from bam (pysam) with same read name.
        read_counts (dict): Dictionary of read count categories (foldback, numt, kept)
        meth_per_read_list (list): List of read methylation levels.
        keep_file_handle (pysam AlignmentFile): File to write kept reads to.
        discard_file_handle (pysam AlignmentFile): File to write discarded reads to.
        meth_prob_likelihood (list): List of meth probability scores.
    Returns:
        None
    """
    if is_foldback(read_parts):
        category, target = 'foldback', discard_file_handle
    elif is_NUMT(read_parts, args.max_sc_threshold, args.max_meth_threshold, meth_per_read_list, meth_prob_likelihood, ref_cons_lengths):
        category, target = 'numt', discard_file_handle
    else:
        category, target = 'kept', keep_file_handle

    read_counts[category] += 1

    for line in read_parts:
        target.write(line)

# Define output file names
if args.output != "":
    prefix = args.output 
else:
    prefix = args.input[:-4]

input_bam = args.input
input_bam_name_sorted = prefix + ".sortedbyname.bam"
filtered_bam = prefix + ".mt.bam"
filtered_bam_name_sorted = prefix + ".mt.sortedbyname.bam"
discarded_bam = prefix + ".discard.bam"
discarded_bam_name_sorted = prefix + ".discard.sortedbyname.bam"

# Perform sorting by name
pysam.sort("-n", "-o", input_bam_name_sorted, input_bam)

fr = pysam.AlignmentFile(input_bam_name_sorted, 'rb')
fw = pysam.AlignmentFile(filtered_bam_name_sorted, 'wb', template= fr)
fd = pysam.AlignmentFile(discarded_bam_name_sorted, 'wb', template= fr)

read_counts = {"foldback": 0, "numt": 0, "kept": 0}
read_parts = []
current_readname = None
percent_meth_per_read = []
meth_prob_likelihood = []
ref_consuming_lengths = []

for line in fr:
    if current_readname == line.query_name:
        read_parts.append(line)
    else:
        if read_parts:  # Process the previous group
            process_read_group(read_parts, read_counts, percent_meth_per_read, fw, fd, meth_prob_likelihood, ref_consuming_lengths)
        read_parts = [line]
        current_readname = line.query_name

# Process the final group
if read_parts:
    process_read_group(read_parts, read_counts, percent_meth_per_read, fw, fd, meth_prob_likelihood, ref_consuming_lengths)

## plot distribution of read methylation
methylation_plot(percent_meth_per_read, meth_prob_likelihood)
plot_ref_consuming_lengths(ref_consuming_lengths)

## close read/write files
fr.close()
fw.close()
fd.close() 

# Perform sorting by name
pysam.sort("-o", filtered_bam, filtered_bam_name_sorted)
pysam.sort("-o", discarded_bam, discarded_bam_name_sorted)

pysam.index(filtered_bam)
pysam.index(discarded_bam)

os.remove(input_bam_name_sorted)
os.remove(filtered_bam_name_sorted)
os.remove(discarded_bam_name_sorted)

print("Number of mitochondrial reads: ", read_counts["kept"])
print("Number of reads discarded as foldback: ", read_counts["foldback"])
print("Number of reads discarded as NUMTs: ", read_counts["numt"])