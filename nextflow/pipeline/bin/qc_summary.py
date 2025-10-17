#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--mosdepth_file", help="Mosdepth summary file.", required=True)
    parser.add_argument("-r", "--read_lengths_file", help="File of MT read lengths", required=True)
    parser.add_argument("-m", "--meth_file", help="Methylation output file.", required=True)
    parser.add_argument("-n", "--nuclear_coverage_file", help="Mosdepth summary file for nuclear genome", required=True)
    parser.add_argument("-s", "--sample_id", help="Sample ID to use in output", required=True)
    return parser.parse_args()
args = get_args()

# === Coverage (mosdepth summary) ===
mosdepth_sum_df = pd.read_csv(args.mosdepth_file, sep='\t')
mean_cov = mosdepth_sum_df[mosdepth_sum_df['chrom'] == 'MT']['mean'].iloc[0]

nuc_cov_df = pd.read_csv(args.nuclear_coverage_file, sep='\t')
nuc_cov = nuc_cov_df[nuc_cov_df['chrom'] == 'total']['mean'].iloc[0]

# === Read lengths ===
read_lengths = np.loadtxt(args.read_lengths_file, dtype=int)
read_count = len(read_lengths)
avg_length = read_lengths.mean()

# === Methylation stats ===
meth_df = pd.read_csv(args.meth_file, sep="\t")
avg_meth = meth_df.iloc[1:, 6].mean()      # column 7 (0-based idx=6), skip header
meth_gt_1 = (meth_df.iloc[1:, 6] > 0.01).mean()
meth_gt_5 = (meth_df.iloc[1:, 6] > 0.05).mean()
meth_gt_10 = (meth_df.iloc[1:, 6] > 0.1).mean()
meth_gt_20 = (meth_df.iloc[1:, 6] > 0.2).mean()

# === Write output ===
out_file = f"{args.sample_id}.qc_summary.tsv"
with open(out_file, "w") as out:
    out.write(
        "Sample\t"
        "Mito_Read_Count\t"
        "Mito_Coverage\t"
        "Nuclear_Coverage\t"
        "Mean_Read_Length\t"
        "Mean_Meth_Fraction\t"
        "Fraction_Meth_Sites_GT1\t"
        "Fraction_Meth_Sites_GT5\t"
        "Fraction_Meth_Sites_GT10\t"
        "Fraction_Meth_Sites_GT20\n"
    )
    out.write(
        f"{args.sample_id}\t"
        f"{read_count}\t"
        f"{mean_cov}\t"
        f"{nuc_cov}\t"
        f"{avg_length:.2f}\t"
        f"{avg_meth:.4f}\t"
        f"{meth_gt_1:.4f}\t"
        f"{meth_gt_5:.4f}\t"
        f"{meth_gt_10:.4f}\t"
        f"{meth_gt_20:.4f}\n"
    )