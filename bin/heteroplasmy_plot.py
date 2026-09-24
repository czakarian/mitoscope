#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pysam import VariantFile

def get_args():
    parser = argparse.ArgumentParser(description="Plots heteroplasmy levels of called variants across mtDNA gneome.")
    parser.add_argument("-i", "--input", help="File path of input VCF.", required=True)
    parser.add_argument("-o", "--output_prefix", help="Output prefix.", required=True)
    return parser.parse_args()
args = get_args()

def get_variant_status(row):
    if row['Disease_Status'] in ["Cfrm-[P]", "Cfrm-[LP]", "Cfrm-[VUS*]"]:
        return "MITOMAP Cfrm [VUS/LP/P]"
    elif not row['MITOMAP_AF']:
        return "Not in MITOMAP"
    else:
        return "MITOMAP Variant"

def create_heteroplasmy_plot(df):

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.scatterplot(data=df, x='POS', y='AF', hue='DiseaseVariantStatus',
        palette={"MITOMAP Variant": "#00BFC4", "MITOMAP Cfrm [VUS/LP/P]": "#F8766D", "Not in MITOMAP": "grey"},
        s=200, alpha=0.7, ax=ax)

    label_data = df[df['DiseaseVariantStatus'].isin(['MITOMAP Cfrm [VUS/LP/P]', 'Not in MITOMAP'])]
    for _, row in label_data.iterrows():
        ax.text(row['POS'], row['AF'] + 0.03,  # Adjust vertical position
            f"{row['REF']}{row['POS']}{row['ALT']}\n{row['AF']*100:.2f}%", 
            fontsize=14, color='black', ha='center')

    ax.set_xlim(-1000, 17500)
    ax.set_xticks([0, 4000, 8000, 12000, 16000])
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])

    # Add axis labels
    ax.set_xlabel("Position in mtDNA genome", fontsize=24)
    ax.set_ylabel("Heteroplasmy Level", fontsize=24)

    # Customize grid, legend, and text sizes
    ax.grid(False)
    ax.legend(title='', fontsize=18, bbox_to_anchor=(1, 0.5), loc='upper left')
    ax.tick_params(axis='both', which='major', labelsize=20)

    return fig


def read_vcf(input_file):
    rows = []
    vcf = VariantFile(input_file, 'r')

    for record in vcf.fetch():
        sample_data = record.samples[0]
        
        row = {
            "CHROM": record.chrom,
            "POS": record.pos,
            "REF": record.ref,
            "ALT": record.alts[0],
            "Disease": record.info["MITOMAP_Disease"][0] if "MITOMAP_Disease" in record.info else "",
            "Disease_Status": record.info["MITOMAP_Disease_Status"][0] if "MITOMAP_Disease_Status" in record.info else "",
            "MITOMAP_AC": record.info["MITOMAP_GenBank_AC"][0] if "MITOMAP_GenBank_AC" in record.info else "",
            "MITOMAP_AF": record.info["MITOMAP_GenBank_AF"][0] if "MITOMAP_GenBank_AF" in record.info else "",
            "Helix_Counts_Hom": record.info["Helix_Counts_Hom"][0] if "Helix_Counts_Hom" in record.info else "",
            "Helix_Counts_Het": record.info["Helix_Counts_Het"][0] if "Helix_Counts_Het" in record.info else "",
            "Helix_AF_Hom": record.info["Helix_AF_Hom"][0] if "Helix_AF_Hom" in record.info else "",
            "Helix_AF_Het": record.info["Helix_AF_Het"][0] if "Helix_AF_Het" in record.info else "",
            "CSQ": record.info["CSQ"][0],
            "AF": sample_data["HPL"][0]
        }
        rows.append(row)
    df = pd.DataFrame(rows)
    return df


input_file = args.input
output_prefix = args.output_prefix

input_df = read_vcf(input_file)
input_df['DiseaseVariantStatus'] = input_df.apply(get_variant_status, axis=1)

fig = create_heteroplasmy_plot(input_df)
fig.savefig(f"{output_prefix}.heteroplasmy.png", dpi=300, bbox_inches='tight')


# input_df[['gnomAD_AC_het', 'gnomAD_AC_hom', 'gnomAD_AN']] = input_df['CSQ'].str.split('|', expand=True)[[25,26,29]]
# cols = ['gnomAD_AC_het', 'gnomAD_AC_hom', 'gnomAD_AN', 'Helix_Counts_Hom', 'Helix_Counts_Het']
# input_df[cols] = input_df[cols].apply(pd.to_numeric, errors='coerce')
# input_df['gnomAD_AF'] = (input_df['gnomAD_AC_het'] + input_df['gnomAD_AC_hom']) / input_df['gnomAD_AN']
# input_df['Helix_AN'] =  195983
# input_df['Helix_AF'] = (input_df['Helix_Counts_Hom'] + input_df['Helix_Counts_Het']) / input_df['Helix_AN']

