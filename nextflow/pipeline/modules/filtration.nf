process FILTER_NUMTS {

    publishDir path: "${params.outdir}/alignments", pattern: "*.{bam,bai}", mode: 'symlink'
    publishDir path: "${params.outdir}/methylation", pattern: "*.png", mode: 'symlink'
    publishDir path: "${params.outdir}/logs", pattern: "*.log", mode: 'symlink'
    container params.python
    tag params.sample_id

    input:
    tuple path(bam_file), path(bam_file_index)

    output:
    tuple path("${bam_file.baseName}.filtered.bam"),path("${bam_file.baseName}.filtered.bam.bai"), emit: mt_filtered_bam
    tuple path("${bam_file.baseName}.discardReads.bam"),path("${bam_file.baseName}.discardReads.bam.bai"), emit: mt_numt_bam
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
    > filter_bam.log 2>&1
    """
}