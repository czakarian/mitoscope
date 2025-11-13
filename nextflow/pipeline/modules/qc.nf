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

process NUCLEAR_COVERAGE {

    publishDir "${params.outdir}/${sample_id}/qc/coverage/nuclear", mode: 'copy'
    container params.mosdepth
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_cram), path(input_cram_index)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.mosdepth.summary.txt"), emit: mosdepth_summary

    script:
    """
    set -euo pipefail

    mosdepth ${sample_id} ${input_cram} --fasta ${ref} --no-per-base --fast-mode
    """
}

process MT_READ_LENGTH {

    publishDir "${params.outdir}/${sample_id}/qc/read_length", mode: 'copy'
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_fastq)

    output:
    tuple val(sample_id), path("${input_fastq.getBaseName(2)}.read_lengths.txt")

    script:
    """
    set -euo pipefail

    zcat ${input_fastq} | sed -n '2~4p'  | awk '{print length(\$0)}' > "${input_fastq.getBaseName(2)}.read_lengths.txt"
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
    
    publishDir "${params.outdir}/${sample_id}/qc/", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mosdepth_summary_file), path(read_lengths_file), path(minimod_file), path(nuclear_coverage_file), path(assembly_info_file)

    output:
    path("${sample_id}.qc_summary.tsv")

    """
    qc_summary.py \
    -c ${mosdepth_summary_file} \
    -r ${read_lengths_file} \
    -m ${minimod_file} \
    -n ${nuclear_coverage_file} \
    -a ${assembly_info_file} \
    -d ${params.num_downsampled_reads} \
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
