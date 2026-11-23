#!/usr/bin/env python3

from pycirclize import Circos
from pycirclize.utils import ColorCycler, load_eukaryote_example_dataset
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(description="Filters out reads in bam file representative of NUMTs and fold-backs (split reads with +/- alignments).")
    parser.add_argument("-i", "--input", help="File path of input bam.", required=True)
    parser.add_argument("-o", "--output", help="Optional output prefix", required=True)
    parser.add_argument("-b", "--bed", help="Bed file for circos plot", required=True)
    return parser.parse_args()
args = get_args()

mt_length = 16569 
scale_factor=200000
input_file=args.input
output_file=args.output

# Initialize Circos with space between chromosomes
circos = Circos.initialize_from_bed(args.bed, space=2)

# Assign Colors
chr_names = [s.name for s in circos.sectors]

#ColorCycler.set_cmap("gist_rainbow")
colors = ColorCycler.get_color_list(len(chr_names))

colors = list(plt.cm.tab20(np.linspace(0, 1, 20)))  # 20 colors
colors += list(plt.cm.tab20b(np.linspace(0, 1, 4)))  # 4 more colors
chr_name2color = {name: color for name, color in zip(chr_names, colors)}
#chr_name2color = {name: colors[i % 24] for i, name in enumerate(chr_names)}

#---- Label chromosomes and draw outer track ----
tick_interval = 2000
positions = [i * scale_factor for i in range(0, mt_length, tick_interval)]
labels = [f"{i//1000}k" for i in range(0, mt_length, tick_interval)]

for sector in circos.sectors:
    label = sector.name.replace("chr", "")
    sector.text(label, size=10)
    track = sector.add_track((95, 100))

    # Add coordinate labels only for chrMT
    if sector.name in ["chrM", "chrMT"]:
        track.xticks(
            positions,
            labels=labels,
            outer=True,
            label_orientation="vertical"
        )
        track.axis(fc="#6baed6", lw=1)
    else:
        color = chr_name2color.get(sector.name, "gray")
        track.axis(fc=color, lw=1)


# ---- Load mappings ----
links_df = pd.read_csv(input_file, sep="\t")
links_df = links_df[links_df['nuc_chrom'] != 'chrM']
links_df = links_df[links_df['nuc_chrom'].isin(chr_names)]

new_rows = []
for _, row in links_df.iterrows():
    start = row["mt_start"]
    end = row["mt_end"]

    if start > end:
        row1 = row.copy()
        row2 = row.copy()

        row1["mt_end"] = mt_length
        row2["mt_start"] = 0

        new_rows.extend([row1, row2])
    else:
        new_rows.append(row)

new_df = pd.DataFrame(new_rows)

# ---- Plot links ----
for _, row in new_df.iterrows():
    region_mt = ("chrMT", int(row["mt_start"]*scale_factor), int(row["mt_end"]*scale_factor))
    region_nuc = (row["nuc_chrom"], int(row["nuc_position"]), int(row["nuc_position"]) + 1)
    # Color based on nuclear chromosome
    color = chr_name2color.get(row["nuc_chrom"], "blue")
    circos.link(region_mt, region_nuc, color=color, alpha=0.4, lw=0.5)

circos.savefig(output_file, dpi=300)
