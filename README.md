# Mitoscope

Mitoscope is a long-read specific mitochondrial analysis workflow to assembly mtDNA, perform heteroplasmic variant calling (of SNVs, indels, and SVs), and characterize non-reference NUMTs. It is built as a fully containerized Nextflow workflow using Singularity images and can be run on either ONT or PacBio data. Raw sequencing data or aligned reads are both accepted as input.

**mitoScope is currently in active development. Open an [issue](https://github.com/czakarian/mitoscope/issues/new) if you find a problem or have a suggestion.**


# Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Outputs](#outputs)

# Installation
## Pre-requisites
- Nextflow
- Singularity
## Setup
```bash
VERSION=0.2.2
wget https://github.com/czakarian/mitoscope/archive/refs/tags/v${VERSION}.tar.gz
tar xvf v${VERSION}.tar.gz && rm v${VERSION}.tar.gz
cd mitoscope-${VERSION}/
```

# Usage
```bash
nextflow run /path/to/mitoscope-${VERSION}/main.nf -params-file example_params.yaml
```

### Parameters Specifications 
The required `-params-file` argument accepts a yaml file with the following parameter specifications. See below or `example_params.yaml` for example.

| name                  | description                                                                    |
|-----------------------|------------------------------------------------------------------------------- |
| samplesheet           | Path to a headerless csv file with 1) sample ids & 2) input file paths  |
| input_type            | Input type (options are: 'bam', 'cram', 'fastq', 'fastq.gz').                   |
| is_aligned            | Set to true if input is aligned, otherwise set to false.                       |
| platform              | Specify 'pb' for PacBio or 'ont' for ONT.                                       |
| outdir                | Path to directory to output results.                                            |
| reference             | Path to reference genome.                                                       |

example_params.yaml
```yaml
samplesheet: "samples.csv" 
input_type: "cram" 
is_aligned: true 
platform: "pb" 
outdir: "/path/to/output/dir/result" 
reference: "/path/to/ref/for/cram/GRCH38.fa" 
```

#### Optional Parameters

To call non-reference NUMT candidates, add `numt_profiling: true` to the yaml file. This functionality is turned off by default to otherwise reduce computational load. 

### Sample Manifest
One of the required parameters in the yaml file is the path to a sample manifest (`samplesheet`). See below or `example_sample_manifest.csv` for example.

samples.csv
```
sample1,/path/to/sample1.bam
sample2,/path/to/sample2.bam
sample3,/path/to/sample3.bam
sample4,/path/to/sample4.bam
```

# Outputs

Mitoscope outputs a directory for each sample with the following subdirectories.

| name                 | description                                                     |
|----------------------|---------------------------------------------------------------- |
| alignments/          | Bam files for retained mtDNA reads and discarded reads.         |
| assembly/            | Assembly output from meta-flye.                                 |
| methylation/         | Methylation calls from modkit/pb-CpG-tools.                     |
| numts/               | Non-reference NUMT candidates.                                  |
| qc/                  | QC metrics.                                                     |
| variants/            | Variant calls (snv/indel, SV)                                   |
| logs/                | Log files.                                                      |
