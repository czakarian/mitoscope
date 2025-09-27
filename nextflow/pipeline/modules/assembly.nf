process BAM_TO_FASTQ_FOR_ASSEMBLY {

    publishDir "${params.outdir}", mode: 'symlink'
    container params.samtools
    tag params.sample_id

    input:
    tuple path(input_bam), path(input_bam_index)

    output:
    path("${input_bam.getBaseName()}.fastq.gz"), emit: filtered_fastq

    script:
    """
    set -euo pipefail

    samtools fastq -T MM,ML -@ ${task.cpus} -0 ${input_bam.getBaseName()}.fastq.gz ${input_bam}
    """
}


process MT_ASSEMBLY {

    publishDir "${params.outdir}/assembly", mode: 'symlink'
    container params.flye
    tag params.sample_id

    input:
    path filtered_fastq
    val platform

    output:
    path("MT_assembly"), emit: assembly_dir

    script:
    def flye_preset = platform == 'ont' ? '--nano-hq' : 
                      platform == 'pb'  ? '--pacbio-hifi':
                      null
    """
    flye --threads ${task.cpus} --meta ${flye_preset} ${filtered_fastq} --out-dir MT_assembly -m ${params.flye_min_overlap}
    """
}
