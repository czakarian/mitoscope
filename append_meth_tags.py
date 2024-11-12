#!/usr/bin/env python

"""This script appends MM/ML tags to the filtered MT reads bam. Upon conversion of the original input bam to fastq for mitoscope processing, 
only sequences are pulled so methylation info is lost in that step. Here we use the reads from the original bam (filtered down to only MT reads) and
append the MM/ML tags to the MT.ref.filtered.bam"""

import argparse
import pysam


def get_args():
    parser = argparse.ArgumentParser(description="Append MM/ML tags back to processed MT read bam.")
    parser.add_argument("-i", "--inputbam", help="File path of bam to append MM/ML tags to.", required=True)
    parser.add_argument("-m", "--methbam", help="File path of bam with MM/ML tags.", required=True)
    return parser.parse_args()

args = get_args()

fr1 = pysam.AlignmentFile(args.inputbam, 'rb')
fr2 = pysam.AlignmentFile(args.methbam, 'rb')
fw = pysam.AlignmentFile(args.inputbam[:-4] + ".withMeth.bam", 'wb', template= fr1)

methtags = {}
## parse the MT subset bam with MM/ML tags to pull read names + MM/ML tags
for line in fr2:
    if line.has_tag("MM") and line.has_tag("ML"):
        readname = line.query_name
        mmtag = line.get_tag("MM")
        mltag = line.get_tag("ML")
        methtags[readname] = {"MM": mmtag, "ML": mltag}

## append MM/ML tags to matching reads in MT.ref.filtered.bam
for line in fr1:
    readname = line.query_name
    if readname in methtags:
        line.set_tag("MM", methtags[readname]["MM"])
        line.set_tag("ML", methtags[readname]["ML"])
    fw.write(line)

## close read/write files
fr1.close()
fr2.close()
fw.close()