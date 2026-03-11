#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-c", "--mito_coverage_file", help="Mito coverage file.", required=True)
    parser.add_argument("-r", "--read_lengths_file", help="File of MT read lengths", required=True)
    parser.add_argument("-n", "--nuclear_coverage_file", help="Mosdepth summary file for nuclear genome", required=False)
    parser.add_argument("-g", "--haplogroup_file", help="Haplogrep output file.", required=True)
    parser.add_argument("-k", "--kmer_read_count", help="Integer value representing number of k-mer selected reads.", required=True)
    parser.add_argument("-m", "--chrM_read_count", help="Integer value representing total number of chrM derived reads.", required=False)
    parser.add_argument("-s", "--sample_id", help="Sample ID to use in output", required=True)
    return parser.parse_args()
args = get_args()

# === Coverage (mosdepth summary) ===
mito_cov_df = pd.read_csv(args.mito_coverage_file, sep='\t')
mito_cov = mito_cov_df [mito_cov_df ['chrom'] == 'MT']['mean'].iloc[0]

# == Nuc coverage & MtDNA copies ===
if args.nuclear_coverage_file:
    nuc_cov_df = pd.read_csv(args.nuclear_coverage_file, sep='\t')
    nuc_cov = nuc_cov_df[nuc_cov_df['chrom'] == 'total']['mean'].iloc[0]
    mtdna_cn = (mito_cov / nuc_cov) * 2

# === Read lengths ===
read_lengths = np.loadtxt(args.read_lengths_file, dtype=int)
read_count = len(read_lengths)
avg_length = read_lengths.mean()
med_length = np.median(read_lengths)

# == Haplogroup ===
haplogroup = pd.read_csv(args.haplogroup_file, sep ='\t').iloc[0,1]

# === Methylation stats ===
# meth_df = pd.read_csv(args.meth_file, sep="\t")
# avg_meth = meth_df.iloc[1:, 6].mean()      # column 7 (0-based idx=6), skip header

# === Write output ===
out_file = f"{args.sample_id}.qc_summary.tsv"
with open(out_file, "w") as out:
    if args.nuclear_coverage_file:
        out.write(
            "Sample\t"
            "ChrM_Read_Count\t"
            "Kmer_Read_Count\t"
            "Mito_Read_Count\t"
            "Mean_Read_Length\t"
            "Median_Read_Length\t"
            "Mito_Coverage\t"
            "Nuclear_Coverage\t"
            "mtDNA_CN\t"
            "Haplogroup\n"
        )
        out.write(
            f"{args.sample_id}\t"
            f"{args.chrM_read_count}\t"
            f"{args.kmer_read_count}\t"
            f"{read_count}\t"
            f"{avg_length:.2f}\t"
            f"{med_length:.2f}\t"
            f"{mito_cov}\t"
            f"{nuc_cov}\t"
            f"{mtdna_cn}\t"
            f"{haplogroup}\n"
        )
    else:
        out.write(
            "Sample\t"
            "Kmer_Read_Count\t"
            "Mito_Read_Count\t"
            "Mean_Read_Length\t"
            "Median_Read_Length\t"
            "Mito_Coverage\t"
            "Haplogroup\n"
        )
        out.write(
            f"{args.sample_id}\t"
            f"{args.kmer_read_count}\t"
            f"{read_count}\t"
            f"{avg_length:.2f}\t"
            f"{med_length:.2f}\t"
            f"{mito_cov}\t"
            f"{haplogroup}\n"
        )
        