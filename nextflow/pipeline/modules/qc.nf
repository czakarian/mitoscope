process MT_COVERAGE {

    publishDir "${params.outdir}/${sample_id}/qc/coverage", mode: 'copy'
    container params.mosdepth
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.per-base.bed.gz"), emit: per_base_bed
    tuple val(sample_id), path("${input_bam.getBaseName()}.mosdepth.summary.txt"), emit: mosdepth_summary

    script:
    """
    set -euo pipefail

    mosdepth ${input_bam.getBaseName()} ${input_bam}
    """

}

process MT_READ_LENGTH {

    publishDir "${params.outdir}/${sample_id}/qc/read_length", mode: 'copy'
    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.read_lengths.txt")

    script:
    """
    set -euo pipefail

    samtools view ${input_bam} | awk '{print length(\$10)}' > "${input_bam.getBaseName()}.read_lengths.txt"
    """

}

process COVERAGE_PLOT {

    publishDir "${params.outdir}/${sample_id}/qc/coverage", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(per_base_bed)

    output:
    path("${per_base_bed.getBaseName(3)}.mitochondrial_coverage.png")

    script:
    """
    set -euo pipefail

    qc_plots.py --plot coverage --input ${per_base_bed} --outprefix ${per_base_bed.getBaseName(3)}

    """
}


process READ_LENGTH_PLOT {

    publishDir "${params.outdir}/${sample_id}/qc/read_length", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(read_length_file)

    output:
    path("${read_length_file.getBaseName(2)}.read_length_distribution.png")

    script:
    """
    set -euo pipefail

    qc_plots.py --plot read_length --input ${read_length_file} --outprefix ${read_length_file.getBaseName(2)}
    """
}

process QC_SUMMARY {
    
    // publishDir "${params.outdir}/${sample_id}/qc/", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mosdepth_summary_file)
    tuple val(sample_id), path(read_lengths_file)
    tuple val(sample_id), path(minimod_file)

    output:
    path("${sample_id}.qc_summary.tsv")

    """
    qc_summary.py \
    -c ${mosdepth_summary_file} \
    -r ${read_lengths_file} \
    -m ${minimod_file} \
    -s ${sample_id}
    """
}


process COMBINE_QC_SUMMARY {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    path summary_files

    output:
    path "qc_summary.tsv"

    script:
    """
    # grab header from the first file
    head -n 1 \$(ls ${summary_files} | head -n 1) > qc_summary.tsv

    # append all data rows (skip header)
    for f in ${summary_files}; do
        tail -n +2 \$f >> qc_summary.tsv
    done
    """
}
