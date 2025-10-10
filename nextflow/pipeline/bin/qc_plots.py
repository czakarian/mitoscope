#!/usr/bin/env python3
"""QC plots"""

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-p", "--plot", help="Plot type <coverage|read_length|methylation>", required=True)
    parser.add_argument("-i", "--input", help="Input file", required=True)
    parser.add_argument("-o", "--outprefix", help="Output prefix", required=True)
    return parser.parse_args()
args = get_args()

def plot_coverage(input, outprefix):
    # Load mosdepth per-base bedGraph (merged intervals)
    df = pd.read_csv(
        input, 
        sep='\t', 
        header=None, 
        names=['chr', 'start', 'end', 'depth'],
        compression='gzip'
    )

    # For plotting, use mid-point of each interval as x-axis
    df['midpoint'] = (df['start'] + df['end']) / 2

    plt.figure(figsize=(12, 4))
    plt.step(df['midpoint'], df['depth'], where='mid', color='blue')
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.title('Mitochondrial coverage from mosdepth')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outprefix + '.mitochondrial_coverage.png', dpi=300)
    plt.close()


def plot_read_lengths(input, outprefix):
    # Load read lengths
    read_lengths = pd.read_csv(input, header=None, names=['length'])

    plt.figure(figsize=(10,6))
    sns.histplot(read_lengths['length'], bins=100, kde=True, color='blue', edgecolor=None)
    plt.title('Read Length Distribution')
    plt.xlabel('Read Length (bp)')
    plt.ylabel('Count')
    #plt.tight_layout()
    plt.savefig(outprefix + '.read_length_distribution.png', dpi=300)
    plt.close()

def plot_methylation(input, outprefix):

    # Load minimod file
    df = pd.read_csv(input, sep='\t').sort_values('start')

    df_aligned = df.copy()
    # shift the '-' strand positions by -1 to merge across strands 
    df_aligned.loc[df_aligned["strand"] == "-", "start"] -= 1
    df_combined = df_aligned.groupby("start", as_index=False).agg({"n_called": "sum", "n_mod": "sum"})
    # recompute methylation frequency
    df_combined["freq"] = df_combined["n_mod"] / df_combined["n_called"]

    # smoothed frequency
    df_combined = df_combined.sort_values("start")
    df_combined["freq_smooth"] = df_combined["freq"].rolling(window=10, center=True, min_periods=1).mean()
    df_combined['index'] = range(0, len(df_combined))

    ## output strand merged tsv?

    ## plot %meth by position and by CpG index
    fig, axes = plt.subplots(3, 1, figsize=(15, 8), sharey=False)
    sns.scatterplot(data=df_combined,x="start",y="n_called",s=8,ax=axes[0])
    axes[0].set_ylabel("Coverage")
    axes[0].set_xlabel("MT Genome Position")
    #axes[0].set_ylim(0, 1)
    sns.lineplot(ax=axes[1], data=df_combined, x="index", y="freq", linewidth=1.5)
    axes[1].set_xlabel("CpG Sites")
    axes[1].set_ylabel("% Methylation")
    axes[1].set_ylim(0, 1)
    sns.lineplot(ax=axes[2], data=df_combined, x="index", y="freq_smooth", linewidth=1.5)
    axes[2].set_xlabel("CpG Sites")
    axes[2].set_ylabel("% Methylation (smoothed)")
    axes[2].set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(outprefix + ".meth_by_site.png", dpi=300)
    plt.show()

    ## plot histogram of %meth frequencies
    plt.figure(figsize=(6, 4))
    sns.histplot(df['freq'], bins=50, kde=True, color='blue', edgecolor=None)
    plt.xlabel('Fraction mtDNA methylation')
    plt.ylabel('# of CpG Sites')
    plt.xlim(0,1)
    plt.tight_layout()
    plt.savefig(outprefix + '.meth_freq_histogram.png', dpi=300)
    plt.close()




if args.plot == 'coverage':
    plot_coverage(args.input, args.outprefix)
elif args.plot == 'read_length':
    plot_read_lengths(args.input, args.outprefix)
elif args.plot == 'methylation':
    plot_methylation(args.input, args.outprefix)
else:
    exit('--plot argument invalid. Must choose appropriate plot type.')

