process METH_FREQ {

    publishDir "${params.outdir}/methylation", mode: 'symlink'
    container params.minimod

    tag "${params.sample_id}"

    input:
    tuple path(input_bam), path(input_bam_index)
    path mt_ref

    output:
    path("${input_bam.getBaseName()}.minimod.tsv"), emit: minimod_tsv
    path("${input_bam.getBaseName()}.minimod.bedGraph"), emit: minimod_bedgraph

    script:
    """
    set -euo pipefail

    minimod freq -t ${task.cpus} -m 0.5 -o ${input_bam.getBaseName()}.minimod.tsv ${mt_ref} ${input_bam}
    sort -k2 -n ${input_bam.getBaseName()}.minimod.tsv -o ${input_bam.getBaseName()}.minimod.tsv
    tail +2 ${input_bam.getBaseName()}.minimod.tsv | cut -f 1-3,7  > "${input_bam.getBaseName()}.minimod.bedGraph"

    """
}

process METH_PLOT {

    publishDir "${params.outdir}/methylation", mode: 'symlink'

    tag "${params.sample_id}"

    input:
    path meth_freq_tsv

    output:
    path("${meth_freq_tsv.getBaseName()}.MT_methylation_by_pos.png")
    path("${meth_freq_tsv.getBaseName()}.MT_methylation_by_site.png")
    path("${meth_freq_tsv.getBaseName()}.MT_methylation_freq_histogram.png")

    script:
    """
    set -euo pipefail

    qc_plots.py --plot methylation --input ${meth_freq_tsv} --outprefix ${meth_freq_tsv.getBaseName()}

    """


}