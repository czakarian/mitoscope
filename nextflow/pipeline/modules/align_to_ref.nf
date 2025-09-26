process ALIGN_TO_REF {

    // publishDir "${params.outdir}", mode: 'symlink'
    container params.minimap

    tag "${params.sample_id}"

    input:
    path fastq_file
    val platform
    path minimap_index

    output:
    path "${fastq_file.getBaseName(2)}.ref.sam", emit: mt_sam

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -y -t ${task.cpus} ${minimap_index} ${fastq_file} > ${fastq_file.getBaseName(2)}.ref.sam
    """
}


process SAMTOOLS_SAM_TO_BAM {

    publishDir "${params.outdir}/alignments", mode: 'symlink'
    container params.samtools

    tag "${params.sample_id}"

    input:
    path input_sam

    output:
    path("${input_sam.baseName}.bam")

    script:
    """
    set -euo pipefail
    
    samtools view -@${task.cpus} -o ${input_sam.baseName}.bam ${input_sam} 
    """
}

process SAMTOOLS_SORT {

    publishDir "${params.outdir}/alignments", mode: 'symlink'
    container params.samtools

    tag "${params.sample_id}"

    input:
    path input_bam

    output:
    path("${input_bam.baseName}.sorted.bam")

    script:
    """
    set -euo pipefail
    
    samtools sort -@${task.cpus} -o "${input_bam.baseName}.sorted.bam" "${input_bam}"
    """
}


process SAMTOOLS_INDEX {

    publishDir "${params.outdir}/alignments", mode: 'symlink'
    container params.samtools

    tag "${params.sample_id}"

    input:
    path input_bam

    output:
    tuple path("${input_bam}"), path("${input_bam.baseName}.bam.bai")

    script:
    """
    set -euo pipefail
    
    samtools index -@${task.cpus} ${input_bam}
    """
}







