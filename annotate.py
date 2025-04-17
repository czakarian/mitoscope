#!/usr/bin/env python

"""This script appends annotations from MITOMAP to a VCF file based on position, REF, ALT fields."""

import argparse
import re
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def get_args():
    parser = argparse.ArgumentParser(description="Appends annotations from MITOMAP to a VCF file based on position, REF, ALT fields.")
    parser.add_argument("-i", "--input", help="File path of input VCF (assumes gzipped).", required=True)
    parser.add_argument("-a", "--annotations", help="File path to MITOMAP annotation csv.", required=True)
    parser.add_argument("-c", "--caller", help="Specify VCF caller used to generated VCF - mutserve, mutect2, or baldur.", required=True)
    return parser.parse_args()
args = get_args()

def get_variant_status(source):
    if source == "disease_df":
        return "Disease Variant"
    elif source in {"variant_df", "disease_df,variant_df"}:
        return "General Variant"
    else:
        return "Unknown Variant"

def resolve_af(row):
    if ',' not in row['AF']:
        return row['AF']
    
    af_values = row['AF'].split(',')
    gt = row['GT']
    
    if gt == '0/1':
        return af_values[1]
    elif gt == '1/0':
        return af_values[0]
    else:
        raise ValueError(f"Unexpected genotype {gt} for multiallelic variant.")

def create_heteroplasmy_plot(df):

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.scatterplot(data=df, x='POS', y='AF', hue='DiseaseVariantStatus',
        palette={"General Variant": "#00BFC4", "Disease Variant": "#F8766D", "Unknown Variant": "grey"},
        s=200, alpha=0.7, ax=ax)

    label_data = df[df['DiseaseVariantStatus'].isin(['Disease Variant', 'Unknown Variant'])]
    for _, row in label_data.iterrows():
        ax.text(row['POS'], row['AF'] + 0.03,  # Adjust vertical position (nudge)
            f"{row['REF']}{row['POS']}{row['ALT']}\n{row['AF']*100:.2f}%", 
            fontsize=14, color='black', ha='center')

    ax.set_xlim(-1000, 17500)
    ax.set_xticks([0, 4000, 8000, 12000, 16000])
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

    # Add axis labels
    ax.set_xlabel("Position in mtDNA genome", fontsize=24)
    ax.set_ylabel("Variant Level", fontsize=24)

    # Customize grid, legend, and text sizes
    ax.grid(False)
    ax.legend(title='', fontsize=18, bbox_to_anchor=(1, 0.5), loc='upper left')
    ax.tick_params(axis='both', which='major', labelsize=20)

    return fig


input_file = args.input
anno_file = args.annotations
caller = args.caller

## read in and format MITOMAP anno file
anno_df = pd.read_csv(anno_file)
anno_df['POS'] = anno_df['Position']
anno_df = anno_df.drop(columns=['Position'])

## read in and format input VCF
input_df = pd.read_csv(input_file,
    comment='#', sep='\t', compression='gzip',
    names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]  
)

# remove blacklisted 3107 row
input_df = input_df[input_df.POS != 3107]

if caller == 'mutserve':
    format_fields = ['GT', 'AF', 'BQ', 'DP']
    input_df[format_fields] = input_df['SAMPLE'].str.split(':', expand=True)
elif caller == 'mutect2':
    format_fields = ['GT', 'AD', 'AF', 'DP']
    input_df[format_fields] = input_df['SAMPLE'].str.split(':', expand=True).iloc[:,:4]
elif caller == 'baldur':
    format_fields = ['GT', 'ADF', 'ADR', 'AF', 'FQSE', 'AQ', 'AFLT', 'QAVG', 'FSB', 'QBS']
    input_df[format_fields] = input_df['SAMPLE'].str.split(':', expand=True)
else:
    # Handle case when caller is not defined or another case
    format_column_fields = []

# Handle multiallelics where AF field was not split // currently only handles up to 2 multiallleic variants, revisit for cases of more
input_df['AF'] = input_df.apply(resolve_af, axis=1)

print(input_df)

input_df = input_df.drop(columns=['FORMAT', 'SAMPLE'])
input_df[['POS', 'AF']] = input_df[['POS', 'AF']].apply(pd.to_numeric, errors='raise')


## merge input df with anno df
merged_df = pd.merge(input_df, anno_df, how="left", on=["POS", "REF", "ALT"])
merged_df['DiseaseVariantStatus'] = merged_df['Source'].apply(get_variant_status)

## output annotated file and heteroplasmy plot
input_prefix = os.path.join(os.path.dirname(input_file), re.sub(r"\.vcf.gz$", "", os.path.basename(input_file)))  
merged_df.to_csv(f"{input_prefix}.annotated.txt", sep='\t', index=False)
fig = create_heteroplasmy_plot(merged_df)
fig.savefig(f"{input_prefix}.heteroplasmy.png", dpi=300, bbox_inches='tight')

