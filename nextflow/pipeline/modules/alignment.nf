process ALIGN_TO_REF {

    publishDir path: "${params.outdir}/${sample_id}/alignments/to_ref", pattern: "*.{bam,bai}", mode: 'symlink'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'symlink'
    container params.minimap
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_file)
    val platform
    path minimap_index

    output:
    tuple val(sample_id), path("${fastq_file.getBaseName(2)}.bam"), path("${fastq_file.getBaseName(2)}.bam.bai"), emit: bam
    path("align_to_ref.log")

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -y -t ${task.cpus} ${minimap_index} ${fastq_file} 2> align_to_ref.log | samtools sort -@${task.cpus} -o "${fastq_file.getBaseName(2)}.bam"
    samtools index -@${task.cpus} "${fastq_file.getBaseName(2)}.bam"
    """
}

process ALIGN_TO_ASSEMBLY {

    publishDir path: "${params.outdir}/${sample_id}/alignments/to_assembly", pattern: "*.{bam,bai}", mode: 'symlink'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'symlink'
    container params.minimap
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_file)
    tuple val(sample_id), path(assembly_fasta), path(assembly_fasta_fai)
    val platform

    output:
    tuple val(sample_id), 
        path("${fastq_file.getBaseName(2)}.${assembly_fasta.getBaseName()}.bam"), 
        path("${fastq_file.getBaseName(2)}.${assembly_fasta.getBaseName()}.bam.bai")
    path("align_to_${assembly_fasta.getBaseName()}.log")

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -y -t ${task.cpus} ${assembly_fasta} ${fastq_file} 2> align_to_${assembly_fasta.getBaseName()}.log | samtools sort -@${task.cpus} -o "${fastq_file.getBaseName(2)}.${assembly_fasta.getBaseName()}.bam"
    samtools index -@${task.cpus} "${fastq_file.getBaseName(2)}.${assembly_fasta.getBaseName()}.bam"
    """
}

process ALIGN_ASSEMBLY_TO_REF {

    publishDir path: "${params.outdir}/${sample_id}/alignments/to_ref", pattern: "*.{bam,bai}", mode: 'symlink'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'symlink'
    container params.minimap
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly_fasta), path(assembly_fasta_fai)
    path minimap_index
    val platform

    output:
    tuple val(sample_id),
        path("${sample_id}.${assembly_fasta.getBaseName()}.ref.bam"),
        path("${sample_id}.${assembly_fasta.getBaseName()}.ref.bam.bai"), emit: bam
    path("align_${assembly_fasta.getBaseName()}_to_ref.log")

    script:
    def preset = platform == 'ont' ? 'map-ont' :
                 platform == 'pb'  ? 'map-hifi' :
                 null

    """
    set -euo pipefail

    minimap2 -ax ${preset} -Y -t ${task.cpus} ${minimap_index} ${assembly_fasta} 2> align_${assembly_fasta.getBaseName()}_to_ref.log | samtools sort -@${task.cpus} -o "${sample_id}.${assembly_fasta.getBaseName()}.ref.bam"
    samtools index -@${task.cpus} "${sample_id}.${assembly_fasta.getBaseName()}.ref.bam"
    """
}

// process SAM_TO_BAM {

//     publishDir "${params.outdir}/alignments", mode: 'symlink'
//     container params.samtools
//     tag params.sample_id

//     input:
//     path input_sam

//     output:
//     tuple path("${input_sam.getBaseName()}.bam"), path("${input_sam.getBaseName()}.bam.bai")

//     script:
//     """
//     set -euo pipefail
    
//     samtools sort -@${task.cpus} -o "${input_sam.baseName}.bam" "${input_sam}"
//     samtools index -@${task.cpus} "${input_sam.baseName}.bam"
//     """
// }

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







