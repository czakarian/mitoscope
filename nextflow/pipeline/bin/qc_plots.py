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

    plt.savefig(outprefix + '.read_length_distribution.png', dpi=300)
    plt.show()

def plot_methylation(input, outprefix):
    # Load minimod bedmethyl file 
    df = pd.read_csv(input, sep='\t').sort_values('start')
    df['index'] = range(0, len(df))

    plt.figure(figsize=(15, 4))
    plt.scatter(df['start'], df['freq'], color='blue', s=10)
    plt.xlabel('Position')
    plt.ylabel('% mtDNA methylation')
    plt.ylim(0,1)
    #plt.title('Average methylation levels across MT genome')
    #plt.grid(True)
    plt.tight_layout()
    plt.savefig(outprefix + '.MT_methylation_by_pos.png', dpi=300)
    plt.show()

    plt.figure(figsize=(15, 4))
    plt.bar(df['index'], df['freq'], color='blue', width=1)
    plt.xlabel('CpG Sites')
    plt.ylabel('% mtDNA methylation')
    plt.ylim(0,1)
    #plt.title('Average methylation levels across MT genome')
    #plt.grid(True)
    plt.tight_layout()
    #plt.xticks([])
    plt.savefig(outprefix + '.MT_methylation_by_site.png', dpi=300)
    plt.show()

    plt.figure(figsize=(6, 4))
    plt.hist(df['freq'], bins=50, color='blue')
    plt.xlabel('Fraction mtDNA methylation')
    plt.ylabel('# of CpG Sites')
    plt.xlim(0,1)
    #plt.grid(True)
    plt.savefig(outprefix + '.MT_methylation_freq_histogram.png', dpi=300)
    plt.show()


if args.plot == 'coverage':
    plot_coverage(args.input, args.outprefix)
elif args.plot == 'read_length':
    plot_read_lengths(args.input, args.outprefix)
elif args.plot == 'methylation':
    plot_methylation(args.input, args.outprefix)
else:
    exit('--plot argument invalid. Must choose appropriate plot type.')

