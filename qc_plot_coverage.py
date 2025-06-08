#!/usr/bin/env python
""""""

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="File path of coverage", required=True)
    parser.add_argument("-o", "--outprefix", help="", required=True)
    return parser.parse_args()
args = get_args()

# Load mosdepth per-base bedGraph (merged intervals)
df = pd.read_csv(
    args.input, 
    sep='\t', 
    header=None, 
    names=['chr', 'start', 'end', 'depth'],
    compression='gzip'
)

# For plotting, use mid-point of each interval as x-axis
df['midpoint'] = (df['start'] + df['end']) / 2

plt.figure(figsize=(12, 4))
plt.step(df['midpoint'], df['depth'], where='mid', color='blue')
plt.xlabel('Position on mitochondrial genome')
plt.ylabel('Coverage depth')
plt.title('Mitochondrial coverage from mosdepth bedGraph intervals')
plt.grid(True)
plt.tight_layout()
plt.savefig(args.outprefix + '.mitochondrial_coverage.png', dpi=150)
plt.show()

