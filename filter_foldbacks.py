#!/usr/bin/env python

"""This script removes reads in an aligned bam file that are representative of 'fold-backs',
which appear as a read with multiple alignments appearing on both the + and - strands."""

import argparse
import pysam


def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of fold-backs.")
    parser.add_argument("-i", "--input", help="File path of input bam.", required=True)
    parser.add_argument("-o", "--output", help="File path for filtered output bam.", required=True)
    parser.add_argument("-d", "--discard", help="File path for reads discarded.", required=True)
    return parser.parse_args()

args = get_args()

def get_strand(is_forward):
    if is_forward:
        return '+'
    else:
        return '-'
    
def is_foldback(strands):
    if len(set(strands)) > 1:
        return True
    else:
        return False

## initialize 
read_counts = {"discarded": 0, "kept_sup": 0, "kept_nonsup": 0}

fr = pysam.AlignmentFile(args.input, 'rb')
fw = pysam.AlignmentFile(args.output, 'wb', template= fr)
fd = pysam.AlignmentFile(args.discard, 'wb', template= fr)

## beging parse the input bam
for line in fr:
    if line.has_tag('SA'):
        
        ## get strand for this alignment and rest of strands associated with this read 
        strand = get_strand(line.is_forward)
        all_strands = [strand]
        SA_tag = line.get_tag('SA').split(';')
        for align in SA_tag:
            if align:
                align_strand = align.split(',')[2]
                all_strands.append(align_strand)

        ## check if alignment associated with foldback read
        foldback_status = is_foldback(all_strands)    
        if foldback_status == True:
            read_counts['discarded']+=1
            fd.write(line)
        else:
            read_counts['kept_sup']+=1
            fw.write(line)
    else:
        read_counts["kept_nonsup"]+=1
        fw.write(line)

## close read/write files
fr.close()
fw.close()
fd.close() 

print("Kept non-supplementary alignments: ", read_counts["kept_nonsup"])
print("Kept supplementary alignments: ", read_counts["kept_sup"])
print("Discarded supplementary alignments: ", read_counts["discarded"])