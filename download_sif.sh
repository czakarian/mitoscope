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
curl -o singularity/flye_2.9.6.sif https://depot.galaxyproject.org/singularity/flye%3A2.9.6--py311h2de2dd3_0 -C -
curl -o singularity/minimap2_2.24.sif https://depot.galaxyproject.org/singularity/minimap2%3A2.24--h7132678_1 -C -
curl -o singularity/samtools_1.21.sif https://depot.galaxyproject.org/singularity/samtools%3A1.21--h96c455f_1 -C -
curl -o singularity/bcftools_1.21.sif https://depot.galaxyproject.org/singularity/bcftools%3A1.21--h3a4d415_1 -C - 
curl -o singularity/bedtools_2.31.0.sif https://depot.galaxyproject.org/singularity/bedtools%3A2.31.0--hf5e1c6e_2 -C -
curl -o singularity/ucsc-bedgraphtobigwig_445.sif https://depot.galaxyproject.org/singularity/ucsc-bedgraphtobigwig%3A445--h954228d_0 -C -
curl -o singularity/sniffles_2.6.2.sif https://depot.galaxyproject.org/singularity/sniffles%3A2.6.2--pyhdfd78af_0 -C - 
curl -o singularity/gatk_4.5.0.0.sif https://depot.galaxyproject.org/singularity/gatk4%3A4.5.0.0--py36hdfd78af_0 -C - 
curl -o singularity/mutserve_2.0.3.sif https://depot.galaxyproject.org/singularity/mutserve%3A2.0.3--hdfd78af_0 -C -
curl -o singularity/haplogrep_3.2.2.sif https://depot.galaxyproject.org/singularity/haplogrep3%3A3.2.2--hdfd78af_0 -C -
curl -o singularity/haplocheck_1.3.3.sif https://depot.galaxyproject.org/singularity/haplocheck%3A1.3.3--h2a3209d_2 -C -
curl -o singularity/mosdepth_0.3.8.sif https://depot.galaxyproject.org/singularity/mosdepth%3A0.3.8--hd299d5a_0 -C -
curl -o https://github.com/warp9seq/minimod/releases/download/v0.4.0/minimod-v0.4.0-x86_64-linux-binaries.tar.gz

echo $(date) - End of download.

echo $(date) - Pull singularity image for baldur from quay.io
singularity pull singularity/baldur_1.2.2.sif oras://quay.io/czakarian/mitoscope:baldur_1.2.2
echo $(date) - Done.

echo $(date) - Set up sif execution permission..
chmod a+rx singularity/
chmod a+rx singularity/*.sif
echo $(date) - End of permission set up.

echo $(date) - Installing minimod..
MINIMOD_VERSION=0.4.0
mkdir -p tools/minimod
wget -P tools/minimod https://github.com/warp9seq/minimod/releases/download/v${MINIMOD_VERSION}/minimod-v${MINIMOD_VERSION}-x86_64-linux-binaries.tar.gz
tar -xvf tools/minimod/minimod-v${MINIMOD_VERSION}-x86_64-linux-binaries.tar.gz -C tools/minimod/ && \
 rm tools/minimod/minimod-v${MINIMOD_VERSION}-x86_64-linux-binaries.tar.gz
echo $(date) - End of download.



