process FILTER_NUMTS {

    publishDir path: "${params.outdir}/${sample_id}/alignments/to_ref", pattern: "*.{bam,bai}", mode: 'symlink'
    publishDir path: "${params.outdir}/${sample_id}/methylation", pattern: "*.png", mode: 'symlink'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'symlink'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id),path(bam_file), path(bam_file_index)

    output:
    tuple val(sample_id), path("${bam_file.baseName}.filtered.bam"), path("${bam_file.baseName}.filtered.bam.bai"), emit: filtered_bam
    tuple val(sample_id), path("${bam_file.baseName}.discardReads.bam"), path("${bam_file.baseName}.discardReads.bam.bai"), emit: numt_bam
    path("${bam_file.baseName}.methylation_per_read.png")
    path("${bam_file.baseName}.methylation_likelihood.png")
    path("filter_bam.log")

    script:
    """
    set -euo pipefail
    
    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    filter_bam.py -i ${bam_file} \
    --max_sc_threshold ${params.max_sc_threshold} \
    --max_meth_threshold ${params.max_meth_threshold_per_read} \
    > filter_bam.log
    """
}

process FILTERED_BAM_TO_FASTQ {

    publishDir "${params.outdir}/${sample_id}", mode: 'symlink'
    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.fastq.gz")

    script:
    """
    set -euo pipefail

    samtools fastq -T MM,ML -@ ${task.cpus} -0 ${input_bam.getBaseName()}.fastq.gz ${input_bam}
    """
}