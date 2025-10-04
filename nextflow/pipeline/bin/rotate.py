#!/usr/bin/env python3

import pysam
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

## to dos - check for circularity
## check ~16kb and not 33kb
## more than 1 contig?


def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
    parser.add_argument("-b", "--bam", help="File path of bam for assembly aligned to MT reference.", required=True)
    parser.add_argument("-a", "--assembly_dir", help="File path of flye assembly directory.", required=True)
    return parser.parse_args()
args = get_args()

def process_one_read_contig(read_parts):
    print('hey')

def process_two_read_contig(read_parts):
    reference_start_position = 0  # 0-based in pysam, corresponds to position 1 in SAM
    rotation_point = None
    reverse_strand = False

    ## can i use both reads to confirm the same thing?
    for read in read_parts:
        print(read.is_reverse, read.reference_start, read.reference_end, read.query_alignment_start, read.query_alignment_end)
        # if read.reference_start == reference_start_position:
        #     reverse_strand = read.is_reverse
        #     if reverse_strand:
        #         rotation_point = read.query_alignment_end - read.query_alignment_start
        #     else:
        #         rotation_point = read.query_alignment_start

        #     print(f"Contig: {read.query_name}")
        #     print(f"Rotation point: {rotation_point}")
        #     print(f"Strand: {'reverse' if reverse_strand else 'forward'}")
        #     break

        # if rotation_point is None:
        #     print("Could not find an alignment starting at reference position 1.")
        # else:
        #     record = SeqIO.read(assembly_fasta, "fasta")
        #     seq = str(record.seq)

        #     # If reverse strand, reverse complement and 
        #     if reverse_strand:
        #         rotated_seq = str(Seq(seq[rotation_point:] + seq[:rotation_point]).reverse_complement())
        #     else:
        #         rotated_seq = seq[rotation_point:] + seq[:rotation_point]


        #     rotated_record = SeqRecord(Seq(rotated_seq), 
        #                             id=record.id + "_rotated", 
        #                             description="")
            
        #     out_fasta = assembly_fasta.rsplit(".",1)[0] + "_rotated.fasta"
        #     SeqIO.write(rotated_record, out_fasta, "fasta")

def process_multi_read_contig(read_parts):
    print('hey')


aligned_assembly_bam = args.bam
assembly_fasta = args.assembly_dir + "/assembly.fasta"
assembly_info = args.assembly_dir + "/assembly_info.txt"
reference_start_position = 0  # 0-based in pysam, corresponds to position 1 in SAM
rotation_point = None
reverse_strand = False


## check if any circular contigs, otherwise do something
info_df = pd.read_csv(assembly_info, sep= '\t')
circular_contigs = info_df[info_df['circ.'] == 'Y']
num_circular_contigs = len(circular_contigs)

if(num_circular_contigs == 0):
    exit("No circular contigs found in assembly. Skipping rotation.")
elif(num_circular_contigs > 1):
    print("Multiple circular contigs identified.")
else:
    print("Single circular contig identified.")

##



bam = pysam.AlignmentFile(aligned_assembly_bam, "rb")
read_parts = []

if num_circular_contigs == 1:
    for line in bam:
        read_parts.append(line)

    if len(read_parts) == 1:
        # assembly break is close enough to ref coordinates that it aligned in 1 segment
        process_one_read_contig(read_parts)
    elif len(read_parts) == 2:
        # assembly broke into 2 parts (most common?)
        process_two_read_contig(read_parts)
    else:
        # assembly split into more than 2 parts? 33kb assembly?
        process_multi_read_contig(read_parts)