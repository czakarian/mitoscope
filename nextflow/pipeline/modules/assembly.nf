process MT_ASSEMBLY {

    publishDir path: "${params.outdir}/${sample_id}/", pattern:"MT_assembly",  mode: 'copy'
    container params.flye
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(filtered_fastq)
    val platform

    output:
    tuple val(sample_id), path("MT_assembly"), emit: assembly_dir

    script:
    def flye_preset = platform == 'ont' ? '--nano-hq' : 
                      platform == 'pb'  ? '--pacbio-hifi':
                      null
    """
    flye --threads ${task.cpus} --meta ${flye_preset} ${filtered_fastq} --out-dir MT_assembly -m ${params.flye_min_overlap}
    """
}


process INDEX_ASSEMBLY {

    publishDir "${params.outdir}/${sample_id}/MT_assembly", mode: 'copy'
    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly_fasta)

    output:
    tuple val(sample_id), path(assembly_fasta), path("*.fai")

    script:
    """
    samtools faidx ${assembly_fasta}
    """
}

// process CHECK_CIRCULARITY {



// }

process ROTATE_ASSEMBLY {

    publishDir "${params.outdir}/${sample_id}/MT_assembly", pattern: "*.fasta", mode: 'symlink'
    publishDir path: "${params.outdir}/logs", pattern: "*.log", mode: 'symlink'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly_fasta), path(assembly_fasta_fai)
    tuple val(sample_id), path(assembly_align_ref_bam), path(assembly_align_ref_bam_index)
    tuple val(sample_id), path(assembly_dir)

    output:
    tuple val(sample_id), path("assembly_rotated.fasta"), optional: true
    // path("rotate_assembly.log"), emit: log

    script:
    """
    rotate.py --bam ${assembly_align_ref_bam} --assembly ${assembly_fasta}
    """
}