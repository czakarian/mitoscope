#!/bin/bash

module load modules modules-init modules-gs
module load bcftools/1.21 htslib/1.21

helix_file="/net/nwgc/vol1/home/czaka/tools/mitoscope/resources/annotations/helix/HelixMTdb_20200327.tsv"

OUTDIR=$(dirname ${helix_file})
FILENAME=$(basename ${helix_file} .tsv)

/net/nwgc/vol1/home/czaka/tools/mitoscope/resources/annotations/helix/format_helix.py \
    -i ${helix_file} \
    -o ${OUTDIR}/${FILENAME}.formatted.tsv

tail +2 ${OUTDIR}/${FILENAME}.formatted.tsv > ${OUTDIR}/${FILENAME}.formatted.noheader.tsv

bgzip -c ${OUTDIR}/${FILENAME}.formatted.noheader.tsv > ${OUTDIR}/${FILENAME}.formatted.noheader.tsv.gz
tabix -f -s 1 -b 2 -e 2 ${OUTDIR}/${FILENAME}.formatted.noheader.tsv.gz