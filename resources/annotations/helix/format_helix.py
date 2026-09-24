#!/usr/bin/env python3

import argparse
import pandas as pd
import ast

def get_args():
    parser = argparse.ArgumentParser(description="Reformat Helix table for use with bcftools annotate")
    parser.add_argument("-i", "--input", help="File path to input HelixMTdb table", required=True)
    parser.add_argument("-o", "--output", help="File path for output.", required=True)
    return parser.parse_args()
args = get_args()


helix_df = pd.read_csv(args.input, sep='\t')
helix_df['alleles'] = helix_df['alleles'].apply(ast.literal_eval)
helix_df["allele_count"] = helix_df["alleles"].apply(len)
helix_df = helix_df[helix_df['allele_count'] == 2]
helix_df[['chrom', 'pos']] = helix_df['locus'].str.split(':', expand=True)
helix_df[['ref', 'alt']] = pd.DataFrame(helix_df['alleles'].tolist(), index=helix_df.index)
helix_df['chrom'] = 'MT'

target_columns = ['chrom', 'pos', 'ref', 'alt', 'counts_hom', 'AF_hom', 'counts_het', 'AF_het', 'max_ARF', 'mean_ARF']
helix_df[target_columns].to_csv(args.output, sep='\t', index=False)
