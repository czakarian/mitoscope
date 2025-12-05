# mitoScope

mitoScope is long-read specific mitochondrial analysis tool that performs NUMT filtration, de novo assembly, variant calling (SNV, indel, SV), and methylation calling. It is built as a fully containerized Nextflow workflow using Singularity images and can be run on either ONT or PacBio data. Raw sequencing data or aligned reads are both accepted as input.

**mitoScope is currently in active development. Open an [issue](https://github.com/czakarian/mitoscope) if you find a problem or have a suggestion.**


# Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)

# Installation
## Pre-requisites
- Nextflow
- Singularity
## Setup
```bash
VERSION=0.1.0
wget https://github.com/czakarian/mitoscope/archive/refs/tags/v${VERSION}.tar.gz
tar xvf mitoscope-${VERSION}.tar.gz
cd mitoscope-${VERSION}/
```

# Usage
```bash
# basic usage
nextflow run main.nf -params-file example_params.yaml

```

*example_params.yaml*
```
samplesheet: "samples.csv" # path to 2-column headerless csv file with: sample_id,sample_file
input_type: 'cram' # set to either 'bam', 'cram', 'fastq', or 'fastq.gz'
is_aligned: true # set to true if aligned bam or cram, otherwise false for raw input
platform: "pb" # pb or ont
outdir: "/path/to/output/dir/result" # path to an existing output directory 
reference: "/path/to/ref/for/cram/GRCH38.fa" # path to reference genome (only required for bam/cram inputs, otherwise set to null)

```

*samples.csv*
```
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
sample3,/path/to/sample3.bam
sample4,/path/to/sample4.bam
```
