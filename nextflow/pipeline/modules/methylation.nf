process METH_FREQ {

    publishDir "${params.outdir}/${sample_id}/methylation", mode: 'copy'
    container params.minimod
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)
    tuple path(mt_ref), path(mt_ref_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.tsv"), emit: minimod_tsv
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.bedGraph"), emit: minimod_bedgraph

    script:
    """
    set -euo pipefail

    minimod freq -t ${task.cpus} \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.tsv \
    ${mt_ref} ${input_bam}
    
    sort -k2 -n ${input_bam.getBaseName()}.minimod.tsv -o ${input_bam.getBaseName()}.minimod.tsv
    tail +2 ${input_bam.getBaseName()}.minimod.tsv | cut -f 1-3,7  > "${input_bam.getBaseName()}.minimod.bedGraph"

    """
}

process METH_PLOT {

    publishDir "${params.outdir}/${sample_id}/methylation", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(meth_freq_tsv)

    output:
    path("${meth_freq_tsv.getBaseName()}.meth_by_site.png")
    path("${meth_freq_tsv.getBaseName()}.meth_freq_histogram.png")

    script:
    """
    set -euo pipefail

    qc_plots.py --plot methylation --input ${meth_freq_tsv} --outprefix ${meth_freq_tsv.getBaseName()}
    """


}