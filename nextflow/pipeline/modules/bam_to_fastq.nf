process BAM_TO_FASTQ {

    container params.samtools

    tag "${params.sample_id}"

    input:
    path bam_file
    val sample_id

    output:
    path "${sample_id}.fastq"

    script:
    """
    set -euo pipefail
    samtools fastq -T MM,ML -@ ${task.cpus} ${bam_file} | tr '\t' ' ' > ${sample_id}.fastq
    """
}

process CRAM_TO_FASTQ {

    container params.samtools

    tag "${params.sample_id}"

    input:
    path cram_file
    path ref
    val sample_id

    output:
    path "${sample_id}.fastq"

    script:
    """
    set -euo pipefail
    samtools view -b -T ${ref} ${cram_file} | samtools fastq -T MM,ML -@ ${task.cpus} | tr '\t' ' ' > ${sample_id}.fastq
    """
}


process COMPRESS_FASTQ {

    publishDir "${params.outdir}", mode: 'symlink'
    container params.pigz

    tag "${params.sample_id}"

    input:
    path input_fastq

    output:
    path "${input_fastq}.gz"

    script:
    """
    set -euo pipefail
    pigz -p ${task.cpus} -c ${input_fastq} > "${input_fastq}.gz"
    """
}

