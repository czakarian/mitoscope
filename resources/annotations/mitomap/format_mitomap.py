#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
from pysam import VariantFile

def get_args():
    parser = argparse.ArgumentParser(description="Appends annotations from MITOMAP to a VCF file based on position, REF, ALT fields.")
    parser.add_argument("-d", "--annotation_disease", help="File path to MITOMAP disease VCF.", required=True)
    parser.add_argument("-p", "--annotations_polymorphisms", help="File path to MITOMAP polymorphisms VCF.", required=True)
    parser.add_argument("-o", "--output", help="File path for output.", required=True)
    return parser.parse_args()
args = get_args()

disease_file = args.annotation_disease
poly_file = args.annotations_polymorphisms
output_file = args.output

def read_mitomap_disease_vcf(input_file):
    rows = []
    vcf = VariantFile(input_file, 'r')

    for record in vcf.fetch():
        row = {
            "CHROM": record.chrom,
            "POS": record.pos,
            "REF": record.ref,
            "ALT": record.alts[0],
            "GenBank_AC": record.info["AC"][0],
            "GenBank_AF": record.info["AF"][0],
            "Disease": record.info["Disease"][0],
            "Disease_Status": record.info["DiseaseStatus"][0],
        }

        rows.append(row)

    df = pd.DataFrame(rows)
    return df

def read_mitomap_poly_vcf(input_file):
    rows = []
    vcf = VariantFile(input_file, 'r')

    for record in vcf.fetch():
        row = {
            "CHROM": record.chrom,
            "POS": record.pos,
            "REF": record.ref,
            "ALT": record.alts[0],
            "GenBank_AC": record.info["AC"][0],
            "GenBank_AF": record.info["AF"][0],
        }

        rows.append(row)

    df = pd.DataFrame(rows)
    return df


disease_df = read_mitomap_disease_vcf(disease_file)
poly_df = read_mitomap_poly_vcf(poly_file)

anno_df = pd.merge(disease_df, poly_df, how="outer", on=["CHROM", "POS", "REF", "ALT"])
anno_df['GenBank_AC'] = np.where(anno_df['GenBank_AC_x'].isna(), anno_df['GenBank_AC_y'], anno_df['GenBank_AC_x'])
anno_df['GenBank_AF'] = np.where(anno_df['GenBank_AF_x'].isna(), anno_df['GenBank_AF_y'], anno_df['GenBank_AF_x'])
anno_df = anno_df.drop(columns=['GenBank_AC_x', 'GenBank_AC_y', 'GenBank_AF_x','GenBank_AF_y'])

anno_df.to_csv(output_file, sep='\t', index=False)