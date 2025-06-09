#!/usr/bin/env python
"""QC plots"""

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-p", "--plot", help="Plot type <coverage|read_length>", required=True)
    parser.add_argument("-i", "--input", help="File path of read lengths", required=True)
    parser.add_argument("-o", "--outprefix", help="", required=True)
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
    plt.savefig(outprefix + '.mitochondrial_coverage.png', dpi=150)
    plt.show()


def plot_read_lengths(input, outprefix):
    # Load read lengths
    read_lengths = pd.read_csv(input, header=None, names=['length'])

    plt.figure(figsize=(10,6))
    sns.histplot(read_lengths['length'], bins=100, kde=True, color='skyblue')

    plt.title('Read Length Distribution')
    plt.xlabel('Read Length (bp)')
    plt.ylabel('Count')
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(outprefix + '.read_length_distribution.png', dpi=150)
    plt.show()

if args.plot == 'coverage':
    plot_coverage(args.input, args.outprefix)
elif args.plot == 'read_length':
    plot_read_lengths(args.input, args.outprefix)
else:
    exit('--plot argument invalid. Must choose appropriate plot type.')

