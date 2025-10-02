#!/usr/bin/env python3

import pysam
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
    parser.add_argument("-b", "--bam", help="File path of bam for assembly aligned to MT reference.", required=True)
    parser.add_argument("-a", "--assembly", help="File path of assembly fasta file.", required=True)
    return parser.parse_args()
args = get_args()

aligned_assembly_bam = args.bam
assembly_fasta = args.assembly
reference_start_position = 0  # 0-based in pysam, corresponds to position 1 in SAM
rotation_point = None
reverse_strand = False

# Open BAM
bam = pysam.AlignmentFile(aligned_assembly_bam, "rb")

for read in bam.fetch():
    if read.reference_start == reference_start_position and not read.is_unmapped:
        reverse_strand = read.is_reverse
        if reverse_strand:
            rotation_point = read.query_alignment_end - read.query_alignment_start
        else:
            rotation_point = read.query_alignment_start

        print(f"Contig: {read.query_name}")
        print(f"Rotation point: {rotation_point}")
        print(f"Strand: {'reverse' if reverse_strand else 'forward'}")
        break

if rotation_point is None:
    print("Could not find an alignment starting at reference position 1.")
else:
    record = SeqIO.read(assembly_fasta, "fasta")
    seq = str(record.seq)

    # If reverse strand, reverse complement and 
    if reverse_strand:
        rotated_seq = str(Seq(seq[rotation_point:] + seq[:rotation_point]).reverse_complement())
    else:
        rotated_seq = seq[rotation_point:] + seq[:rotation_point]


    rotated_record = SeqRecord(Seq(rotated_seq), 
                               id=record.id + "_rotated", 
                               description="")
    
    out_fasta = assembly_fasta.rsplit(".",1)[0] + "_rotated.fasta"
    SeqIO.write(rotated_record, out_fasta, "fasta")

