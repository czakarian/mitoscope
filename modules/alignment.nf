process ALIGN_TO_REF {

    // publishDir path: "${params.outdir}/${sample_id}/alignments/", pattern: "*.{bam,bai}", mode: 'copy'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'copy'
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

    // publishDir path: "${params.outdir}/${sample_id}/alignments/", pattern: "*.{bam,bai}", mode: 'copy'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'copy'
    container params.minimap
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_file), path(assembly_fasta), path(assembly_fasta_fai)
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

    publishDir path: "${params.outdir}/${sample_id}/alignments/", pattern: "*.{bam,bai}", mode: 'copy'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern: "*.log", mode: 'copy'
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