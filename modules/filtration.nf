process FILTER_NUMTS {

    publishDir path: "${params.outdir}/${sample_id}/alignments/", pattern: "*.{bam,bai}", mode: 'copy'
    publishDir path: "${params.outdir}/${sample_id}/qc/methylation/", pattern: "*.png", mode: 'copy'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id),path(bam_file), path(bam_file_index)

    output:
    tuple val(sample_id), path("${sample_id}.mt.bam"), path("${sample_id}.mt.bam.bai"), emit: filtered_bam
    tuple val(sample_id), path("${sample_id}.discard.bam"), path("${sample_id}.discard.bam.bai")
    path("${sample_id}.methylation_per_read.png")
    path("${sample_id}.methylation_likelihood.png")
    //path("${sample_id}.unaligned_sc_hist.png")
    //path("${sample_id}.ref_consuming_hist_kde.png")
    path("filter_bam.log")

    script:
    """
    set -euo pipefail
    
    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    filter_bam.py -i ${bam_file} \
    --max_sc_threshold ${params.max_sc_threshold} \
    --max_meth_threshold ${params.max_meth_threshold_per_read} > filter_bam.log
    """
}

process FILTERED_BAM_TO_FASTQ {

    // publishDir "${params.outdir}/${sample_id}", mode: 'copy'
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