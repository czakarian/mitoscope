#!/bin/bash

REF="/net/nwgc/vol1/home/czaka/tools/mitoscope/resources/MT.fasta"
INDIR="/net/nwgc/vol1/home/czaka/tools/mitoscope/resources/annotations/mitomap/"
DISEASE_VCF="${INDIR}/disease.vcf"
POLY_VCF="${INDIR}/polymorphisms.vcf"


DISEASE_NORM_VCF="${INDIR}/$(basename ${DISEASE_VCF} .vcf).norm.vcf"
bcftools norm --force -cx -f ${REF} --multiallelics -any ${DISEASE_VCF} \
| bcftools norm -cx -f ${REF} --atomize --atom-overlaps . -Oz -o ${DISEASE_NORM_VCF}

POLY_NORM_VCF="${INDIR}/$(basename ${POLY_VCF} .vcf).norm.vcf"
bcftools norm -cx -f ${REF} --multiallelics -any ${POLY_VCF} \
| bcftools norm -cx -f ${REF} --atomize --atom-overlaps . -Oz -o ${POLY_NORM_VCF}

/net/nwgc/vol1/home/czaka/tools/mitoscope/resources/annotations/mitomap/format_mitomap.py \
    -d ${DISEASE_NORM_VCF} \
    -p ${POLY_NORM_VCF} \
    -o ${INDIR}/disease_and_polymorphisms.mitomap.tsv

tail +2 ${INDIR}/disease_and_polymorphisms.mitomap.tsv > ${INDIR}/disease_and_polymorphisms.mitomap.noheader.tsv

bgzip -c ${INDIR}/disease_and_polymorphisms.mitomap.noheader.tsv > ${INDIR}/disease_and_polymorphisms.mitomap.noheader.tsv.gz
tabix -f -s 1 -b 2 -e 2 ${INDIR}/disease_and_polymorphisms.mitomap.noheader.tsv.gz