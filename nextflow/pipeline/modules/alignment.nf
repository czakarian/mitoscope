process ALIGN_TO_REF {

    container params.minimap
    tag params.sample_id

    input:
    path fastq_file
    val platform
    path minimap_index

    output:
    path "${fastq_file.getBaseName(2)}.sam", emit: mt_sam

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -y -t ${task.cpus} ${minimap_index} ${fastq_file} > ${fastq_file.getBaseName(2)}.sam
    """
}

process ALIGN_TO_ASSEMBLY {

    container params.minimap
    tag params.sample_id

    input:
    path fastq_file
    val platform
    path assembly_dir

    output:
    path "${fastq_file.getBaseName(2)}.assembly.sam", emit: mt_sam

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -y -t ${task.cpus} ${assembly_dir}/assembly.fasta ${fastq_file} > ${fastq_file.getBaseName(2)}.assembly.sam
    """
}


process SAMTOOLS_SAM_TO_BAM {

    publishDir "${params.outdir}/alignments", mode: 'symlink'
    container params.samtools
    tag params.sample_id

    input:
    path input_sam

    output:
    path("${input_sam.baseName}.bam")

    script:
    """
    set -euo pipefail
    
    samtools sort -@${task.cpus} -o "${input_sam.baseName}.bam" "${input_sam}"
    """
}

process SAMTOOLS_INDEX {

    publishDir "${params.outdir}/alignments", mode: 'symlink'
    container params.samtools
    tag params.sample_id

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







