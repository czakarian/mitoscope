#!/usr/bin/env python3

import argparse
import pandas as pd

# def get_args():
#     parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
#     parser.add_argument("-g", "--cpg", help="", required=True)
#     parser.add_argument("-h", "--cph", help="", required=True)
#     parser.add_argument("-o", "--output", help="Optional output prefix", required=True)
#     return parser.parse_args()
# args = get_args()


cpg="/net/nwgc/vol1/nobackup/czaka/mitoscope/smaht/benchmark/ont/output/ST004-1Q_BRAIN-ont-uwsc/methylation/ST004-1Q_BRAIN-ont-uwsc.MT.filtered.minimod.mCG.tsv"
cph="/net/nwgc/vol1/nobackup/czaka/mitoscope/smaht/benchmark/ont/output/ST004-1Q_BRAIN-ont-uwsc/methylation/ST004-1Q_BRAIN-ont-uwsc.MT.filtered.minimod.mCH.tsv"


df_g = pd.read_csv(cpg, sep='\t')
df_h = pd.read_csv(cph, sep='\t')


merged_df = pd.merge(df_g, df_h, on=df_g.columns.to_list(), how='right', indicator=True)

merged_df[merged_df['_merge'] == 'right_only'].to_csv('check.tsv', sep='\t', index=False)
#unique_to_df1_combinations = merged_df[merged_df['_merge'] == 'left_only'][['ColA', 'ColB']]
#print(f"Combinations unique to df1: \n{unique_to_df1_combinations}")