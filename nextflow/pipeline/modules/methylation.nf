process METH_FREQ {

    publishDir "${params.outdir}/${sample_id}/methylation", mode: 'copy'
    container params.minimod
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)
    tuple path(mt_ref), path(mt_ref_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.mCG.tsv"), emit: minimod_tsv_mcg
    tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.hCG.tsv"), emit: minimod_tsv_hcg
    //tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.mCH.tsv"), emit: minimod_tsv_mch
    //tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.hCH.tsv"), emit: minimod_tsv_hch
    //tuple val(sample_id), path("${input_bam.getBaseName()}.minimod.aA.tsv"), emit: minimod_tsv_aA
    //path('minimod_summary.txt')

    script:
    """
    set -euo pipefail

    minimod summary ${input_bam} > minimod_summary.txt

    minimod freq -t ${task.cpus} -c "m[CG]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.mCG.tsv \
    ${mt_ref} ${input_bam}

    minimod freq -t ${task.cpus} -c "m[*]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.mCH.tsv \
    ${mt_ref} ${input_bam}

    minimod freq -t ${task.cpus} -c "h[CG]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.hCG.tsv \
    ${mt_ref} ${input_bam}

    minimod freq -t ${task.cpus} -c "h[*]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.hCH.tsv \
    ${mt_ref} ${input_bam}

    minimod freq -t ${task.cpus} -c "a[A]" \
    -m ${params.meth_likelihood_threshold} \
    -o ${input_bam.getBaseName()}.minimod.aA.tsv \
    ${mt_ref} ${input_bam}
    
    sort -k2 -n ${input_bam.getBaseName()}.minimod.mCG.tsv -o ${input_bam.getBaseName()}.minimod.mCG.tsv
    sort -k2 -n ${input_bam.getBaseName()}.minimod.mCH.tsv -o ${input_bam.getBaseName()}.minimod.mCH.tsv
    sort -k2 -n ${input_bam.getBaseName()}.minimod.hCG.tsv -o ${input_bam.getBaseName()}.minimod.hCG.tsv
    sort -k2 -n ${input_bam.getBaseName()}.minimod.hCH.tsv -o ${input_bam.getBaseName()}.minimod.hCH.tsv
    sort -k2 -n ${input_bam.getBaseName()}.minimod.aA.tsv -o ${input_bam.getBaseName()}.minimod.aA.tsv

    """
}

process METH_PLOT {

    publishDir "${params.outdir}/${sample_id}/methylation/plots", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(minimod_tsv_mcg),  path(minimod_tsv_hcg)

    output:
    path("${minimod_tsv_mcg.getBaseName(2)}.meth_by_site.png")
    path("${minimod_tsv_mcg.getBaseName(2)}.meth_freq_histogram.png")

    script:
    """
    set -euo pipefail

    cat ${minimod_tsv_mcg} > combined_minimod_output.tsv
    tail +2 ${minimod_tsv_hcg} >> combined_minimod_output.tsv

    qc_plots.py --plot methylation --input combined_minimod_output.tsv --outprefix ${minimod_tsv_mcg.getBaseName(2)}
    """

}