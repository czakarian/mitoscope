process METH_FREQ {

    publishDir "${params.outdir}/${sample_id}/methylation", mode: 'copy'
    container params.minimod
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)
    tuple path(mt_ref), path(mt_ref_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.mCG.tsv"), emit: minimod_tsv_m
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.hCG.tsv"), emit: minimod_tsv
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.aA.tsv"), emit: minimod_tsv_aA
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.aT.tsv"), emit: minimod_tsv_aT
    // tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.bedGraph"), emit: minimod_bedgraph
    path('minimod_summary.txt')

    script:
    """
    set -euo pipefail

    minimod freq -t ${task.cpus} -c "m[CG]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.mCG.tsv \
    ${mt_ref} ${input_bam}

    minimod freq -t ${task.cpus} -c "h[CG]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.hCG.tsv \
    ${mt_ref} ${input_bam}

    minimod freq -t ${task.cpus} -c "a[A]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.aA.tsv \
    ${mt_ref} ${input_bam}

    minimod freq -t ${task.cpus} -c "a[T]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.aT.tsv \
    ${mt_ref} ${input_bam}
    
    minimod summary ${input_bam} > minimod_summary.txt

    sort -k2 -n ${input_bam.getBaseName()}.minimod.mCG.tsv -o ${input_bam.getBaseName()}.minimod.mCG.tsv
    sort -k2 -n ${input_bam.getBaseName()}.minimod.hCG.tsv -o ${input_bam.getBaseName()}.minimod.hCG.tsv
    sort -k2 -n ${input_bam.getBaseName()}.minimod.aA.tsv -o ${input_bam.getBaseName()}.minimod.aA.tsv
    sort -k2 -n ${input_bam.getBaseName()}.minimod.aT.tsv -o ${input_bam.getBaseName()}.minimod.aT.tsv

    #tail +2 ${input_bam.getBaseName()}.minimod.tsv | cut -f 1-3,7  > "${input_bam.getBaseName()}.minimod.bedGraph"

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
    path("meth.csv")

    script:
    """
    set -euo pipefail

    qc_plots.py --plot methylation --input ${meth_freq_tsv} --outprefix ${meth_freq_tsv.getBaseName()}
    """

}