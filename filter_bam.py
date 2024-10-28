#!/usr/bin/env python

"""This script removes reads in an aligned bam file that are representative of: 
(1) fold-backs which appear as a read with multiple alignments appearing on both the + and - strands,
(2) NUMTs which we classify as reads with >1kb softclipping that isn't a result of wrap-around and won't be a split read."""

import argparse
import pysam
import re

def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
    parser.add_argument("-i", "--input", help="File path of input bam.", required=True)
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

def is_NUMT(cigar):

    result_left = re.match('^([0-9]+)[S]', cigar)
    result_right = re.search('([0-9]+)[S]$', cigar)

    sc_left = int(result_left[1]) if result_left else 0
    sc_right = int(result_right[1]) if result_right else 0

    print(sc_left, sc_right)

    ## if >=1kb softclipping discard as NUMT
    if sc_left > 1000 or sc_right > 1000 or sc_left + sc_right > 1000:
        return True
    else:
        return False

read_counts = {"foldback": 0, "numt": 0, "kept_suppl": 0, "kept_primary": 0, }

fr = pysam.AlignmentFile(args.input, 'rb')
fw = pysam.AlignmentFile(args.input[:-4] + ".filtered.bam", 'wb', template= fr)
fd = pysam.AlignmentFile(args.input[:-4] + ".discardReads.bam", 'wb', template= fr)

## begin parse the input bam
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
            read_counts['foldback']+=1
            fd.write(line)
        else:
            read_counts['kept_suppl']+=1
            fw.write(line)
    else:
        ## check if alignment associated with NUMT
        numt_status = is_NUMT(line.cigarstring)

        if numt_status == True:
            read_counts['numt']+=1
            fd.write(line)
        else:
            read_counts['kept_primary']+=1
            fw.write(line)


## close read/write files
fr.close()
fw.close()
fd.close() 

print("Kept primary alignments: ", read_counts["kept_primary"])
print("Kept supplementary alignments: ", read_counts["kept_suppl"])
print("Discarded supplementary alignments (foldback): ", read_counts["foldback"])
print("Discarded NUMT alignments: ", read_counts["numt"])