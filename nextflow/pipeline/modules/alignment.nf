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
    path assembly_fasta

    output:
    path "${fastq_file.getBaseName(2)}.${assembly_fasta.getBaseName()}.sam", emit: mt_sam

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -y -t ${task.cpus} ${assembly_fasta} ${fastq_file} > ${fastq_file.getBaseName(2)}.${assembly_fasta.getBaseName()}.sam
    """
}

process ALIGN_ASSEMBLY_TO_REF {

    container params.minimap
    tag params.sample_id

    input:
    path assembly_fasta
    path minimap_index
    val platform
    val sample_id

    output:
    path "${sample_id}.${assembly_fasta.getBaseName()}.ref.sam", emit: mt_sam

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -t ${task.cpus} ${minimap_index} ${assembly_fasta} > ${sample_id}.${assembly_fasta.getBaseName()}.ref.sam
    """
}

process SAM_TO_BAM {

    publishDir "${params.outdir}/alignments", mode: 'symlink'
    container params.samtools
    tag params.sample_id

    input:
    path input_sam

    output:
    tuple path("${input_sam.getBaseName()}.bam"), path("${input_sam.getBaseName()}.bam.bai")

    script:
    """
    set -euo pipefail
    
    samtools sort -@${task.cpus} -o "${input_sam.baseName}.bam" "${input_sam}"
    samtools index -@${task.cpus} "${input_sam.baseName}.bam"
    """
}

// process SAMTOOLS_INDEX {

//     publishDir "${params.outdir}/alignments", mode: 'symlink'
//     container params.samtools
//     tag params.sample_id

//     input:
//     path input_bam

//     output:
//     tuple path("${input_bam}"), path("${input_bam.baseName}.bam.bai")

//     script:
//     """
//     set -euo pipefail
    
//     samtools index -@${task.cpus} ${input_bam}
//     """
// }







