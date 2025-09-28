process ALIGNED_BAM_TO_FASTQ {
    // for aligned bams, pull reads from chrM before running kmer selection 
    container params.samtools
    tag params.sample_id

    input:
    tuple path(bam_file), path(bam_file_index)
    val sample_id

    output:
    path "${sample_id}.fastq"

    script:
    """
    set -euo pipefail
    samtools view -h -@ ${task.cpus} ${bam_file} chrM | samtools fastq -T MM,ML -@ ${task.cpus} | tr '\t' ' ' > ${sample_id}.fastq
    """
}

process UNALIGNED_BAM_TO_FASTQ {
    // for unaligned bams, convert full bam to fastq then run kmer selection
    container params.samtools
    tag params.sample_id

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

process ALIGNED_CRAM_TO_FASTQ {

    container params.samtools
    tag params.sample_id

    input:
    tuple path(cram_file), path(cram_file_index)
    val sample_id
    path ref

    output:
    path("${sample_id}.fastq")

    script:
    """
    set -euo pipefail
    samtools view -h -@ ${task.cpus} -T ${ref} ${cram_file} chrM | samtools fastq -T MM,ML -@ ${task.cpus} | tr '\t' ' ' > ${sample_id}.fastq
    """
}

process UNALIGNED_CRAM_TO_FASTQ {

    container params.samtools
    tag params.sample_id

    input:
    path cram_file
    val sample_id
    path ref

    output:
    path("${sample_id}.fastq")

    script:
    """
    set -euo pipefail
    samtools fastq --reference ${ref} -T MM,ML -@ ${task.cpus} ${cram_file} | tr '\t' ' ' > ${sample_id}.fastq
    """
}

process COMPRESS_FASTQ {

    publishDir "${params.outdir}", mode: 'symlink'
    container params.pigz
    tag params.sample_id

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

