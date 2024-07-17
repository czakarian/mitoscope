#!/bin/bash
set -eu -o pipefail

##
## mitoscope
## Copyright (C) 2019-2024 Dr Chia-Lin Wei Laboratory, NWGC
##

echo $(date) - Downloading sif files to singularity subfolder..
[ ! -d "singularity" ] && mkdir "singularity"
curl -o singularity/kmc_3.2.1.sif https://depot.galaxyproject.org/singularity/kmc%3A3.2.1--hf1761c0_2 -C -
curl -o singularity/pigz_2.3.4.sif https://depot.galaxyproject.org/singularity/pigz%3A2.3.4 -C -
curl -o singularity/flye_2.9.1.sif https://depot.galaxyproject.org/singularity/flye%3A2.9.1--py39h6935b12_0 -C -
curl -o singularity/minimap2_2.24.sif https://depot.galaxyproject.org/singularity/minimap2%3A2.24--h7132678_1 -C -
curl -o singularity/samtools_v1.15.1.sif https://depot.galaxyproject.org/singularity/samtools%3A1.15.1--h1170115_0 -C -
curl -o singularity/seqtk_1.3.sif https://depot.galaxyproject.org/singularity/seqtk%3A1.3--h7132678_4 -C -
curl -o singularity/bedtools_2.31.0.sif https://depot.galaxyproject.org/singularity/bedtools%3A2.31.0--hf5e1c6e_2 -C -
curl -o singularity/ucsc-bedgraphtobigwig_445.sif https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig%3A445--h954228d_0 -C -
echo $(date) - End of download.

echo $(date) - Set up sif execution permission..
chmod a+rx singularity/
chmod a+rx singularity/*.sif
echo $(date) - End of permission set up.
