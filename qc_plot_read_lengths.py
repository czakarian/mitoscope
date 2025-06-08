#!/usr/bin/env python
""""""

import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-i", "--input", help="File path of read lengths", required=True)
    parser.add_argument("-o", "--outprefix", help="", required=True)
    return parser.parse_args()
args = get_args()

# Load read lengths
read_lengths = pd.read_csv(args.input, header=None, names=['length'])

plt.figure(figsize=(10,6))
sns.histplot(read_lengths['length'], bins=100, kde=True, color='skyblue')

plt.title('Read Length Distribution')
plt.xlabel('Read Length (bp)')
plt.ylabel('Count')
plt.grid(True)
plt.tight_layout()

plt.savefig(args.outprefix + '.read_length_distribution.png', dpi=150)
plt.show()
