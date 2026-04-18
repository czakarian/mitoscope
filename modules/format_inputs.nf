process ALIGNED_BAM_TO_FASTQ {

    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam_file), path(bam_file_index)

    output:
    tuple val(sample_id), path("${sample_id}.fastq"), emit:fastq
    tuple val(sample_id), path("${sample_id}.chrM_read_count.txt"), emit:chrM_read_count

    script:
    """
    set -euo pipefail

    samtools view -h -@ ${task.cpus} ${bam_file} chrM | samtools fastq -T MM,ML -@ ${task.cpus} | tr '\t' ' ' > ${sample_id}.fastq
    samtools view -@ ${task.cpus} -c -F 2308 ${bam_file} chrM > ${sample_id}.chrM_read_count.txt
    """
}

process UNALIGNED_BAM_TO_FASTQ {
    // for unaligned bams, convert full bam to fastq then run kmer selection
    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.fastq"), emit:fastq

    script:
    """
    set -euo pipefail

    samtools fastq -T MM,ML -@ ${task.cpus} ${bam_file} | tr '\t' ' ' > ${sample_id}.fastq
    """
}

process ALIGNED_CRAM_TO_FASTQ {
    
    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(cram_file), path(cram_file_index)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.fastq"), emit:fastq
    tuple val(sample_id), path("${sample_id}.chrM_read_count.txt"), emit:chrM_read_count


    script:
    """
    set -euo pipefail

    samtools view -h -@ ${task.cpus} -T ${ref} ${cram_file} chrM | samtools fastq -T MM,ML -@ ${task.cpus} | tr '\t' ' ' > ${sample_id}.fastq
    samtools view -@ ${task.cpus} -T ${ref} -c -F 2308 ${cram_file} chrM > ${sample_id}.chrM_read_count.txt
    """
}

process UNALIGNED_CRAM_TO_FASTQ {

    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(cram_file)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.fastq"), emit:fastq

    script:
    """
    set -euo pipefail
    samtools fastq --reference ${ref} -T MM,ML -@ ${task.cpus} ${cram_file} | tr '\t' ' ' > ${sample_id}.fastq
    """
}

process COMPRESS_FASTQ {

    container params.pigz
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_fastq)

    output:
    tuple val(sample_id), path("${input_fastq}.gz")

    script:
    """
    set -euo pipefail
    pigz -p ${task.cpus} -c ${input_fastq} > "${input_fastq}.gz"
    """
}

