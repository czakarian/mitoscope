#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input_fasta", help="Input fasta file.", required=True)
    parser.add_argument("-o", "--output_fasta", help="Rotated fasta file.", required=True)
    parser.add_argument("-s", "--shift", type=int, help="Number of bases to shift seqeunce forward.", required=False, default=8000)
    return parser.parse_args()
args = get_args()

shift = args.shift

record = SeqIO.read(args.input_fasta, "fasta")
rotated_seq = record.seq[shift:] + record.seq[:shift]
new_record = SeqRecord(rotated_seq, id=f"{record.id}_rotated_{shift}", description="")
SeqIO.write(new_record, args.output_fasta, "fasta")