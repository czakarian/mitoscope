#!/usr/bin/env python3 

import pysam
import argparse 

def get_args():
    parser = argparse.ArgumentParser(description="Generate a read x variant matrix.")
    parser.add_argument("-b", "--bam", help="Input bam", required=True)
    parser.add_argument("-v", "--vcf", help="Input VCF.", required=True)
    parser.add_argument("-o", "--output", help="Prefix for for output files.", required=True)
    parser.add_argument("-r", "--ref", help="Mitochondrial reference (rCRS) fasta.", required=True)
    return parser.parse_args()
args = get_args()

class SampleVariant:
    """Store variant information for a single variant (position, ref, alt)."""
    def __init__(self, pos, ref, alt):
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.id = f'{pos}:{ref}>{alt}'

        if len(ref) == 1 and len(alt) == 1:
            self.type = 'snv'
        else:
            self.type = 'indel'

def get_read_indels(read, reference_sequence):
    """
    Return indels in a read as list of 'pos:ref>alt'
    """
    indels = []

    ref_pos = read.reference_start
    query_pos = 0

    for op, length in read.cigartuples:
        if op == 0:  # match / mismatch
            ref_pos += length
            query_pos += length
        elif op == 2:  # Deletion from reference
            ref_seq = reference_sequence[ref_pos-1:ref_pos + length]
            alt_seq = reference_sequence[ref_pos-1]
            indels.append(f'{ref_pos}:{ref_seq}>{alt_seq}')
            ref_pos += length
        elif op == 1:  # Insertion in read
            ref_seq = read.query_sequence[query_pos-1]
            alt_seq = read.query_sequence[query_pos-1:query_pos + length]
            indels.append(f'{ref_pos}:{ref_seq}>{alt_seq}')
            query_pos += length
        elif op == 4:
            query_pos += length

    return indels

def get_read_deletions(read, reference_sequence):
    """
    Return large deletions in a read
    """
    dels = []

    ref_pos = read.reference_start
    query_pos = 0

    for op, length in read.cigartuples:
        if op == 0:  # match / mismatch
            ref_pos += length
        elif op == 2:  # Deletion from reference
            if length >= 45:
                ref_seq = reference_sequence[ref_pos-1:ref_pos + length]
                alt_seq = reference_sequence[ref_pos-1]
                dels.append(f'del{ref_pos}_{ref_pos + length - 1}')
            ref_pos += length
        elif op == 1:  # Insertion in read
            if length >= 45:
                ref_seq = read.query_sequence[query_pos-1]
                alt_seq = read.query_sequence[query_pos-1:query_pos + length]
                #dels.append(f'ins{ref_pos}_{length}')
            query_pos += length
        elif op == 4:
            query_pos += length
    return dels

bam = pysam.AlignmentFile(args.bam, "rb")
vcf = pysam.VariantFile(args.vcf, "rb")
snv_indel_out_file = args.output + ".snv_indel.csv"
deletions_out_file = args.output + ".deletions.csv"

mt_ref_fasta = pysam.FastaFile(args.ref)
mt_ref_seq = mt_ref_fasta.fetch('MT')

## store list of variant objects
variants = []
header = ['sample_id']
for variant in vcf.fetch():
    var_obj = SampleVariant(variant.pos,variant.ref,variant.alts[0])
    header.append(var_obj.id)
    variants.append(var_obj)


## store read ids with associated read records (ie for read ids with multiple records // suppl.)
reads = {}
for read in bam.fetch("MT"):
    if read.is_unmapped:
        continue
    reads.setdefault(read.query_name, []).append(read)


## snvs and indels
fw= open(snv_indel_out_file, 'w')
fw.write(','.join(header) + '\n')

for qname, alignments in reads.items():
    calls = []

    for variant in variants:
        ## default value if site not spanned
        value = '-1'

        for read in alignments:
            if variant.type == 'snv':
                ref_positions = read.get_reference_positions(full_length=True)

                ## check if read spans the target variant
                if variant.pos - 1 in ref_positions:
                    idx = ref_positions.index(variant.pos - 1)
                    base = read.query_sequence[idx]

                    if base == variant.ref:
                        value = '0'
                    elif base == variant.alt:
                        value = '1'
                    else:
                        print(variant.id)
                        ## if other alt value
                        value = '2'
                    break

            elif variant.type == 'indel':
                
                ## check if read spans the target variant
                if read.reference_start <= variant.pos <= read.reference_end:

                    read_indels = get_read_indels(read, mt_ref_seq)

                    if variant.id in read_indels:
                        value = '1'
                    else:
                        value = '0'

                    break


        calls.append(value)

    fw.write(qname + "," + ','.join(calls) + '\n')

fw.close()


## larger deletions (>25)
fw_del = open(deletions_out_file, 'w')

for qname, alignments in reads.items():
    dels = set()

    for read in alignments:
            
        read_dels = get_read_deletions(read, mt_ref_seq)
        dels.update(read_dels)

    if len(dels) > 0:
        fw_del.write(qname + "\t" + ','.join(dels) + '\n')

fw_del.close()