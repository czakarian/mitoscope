process MT_ASSEMBLY {

    publishDir "${params.outdir}/${sample_id}/", mode: 'symlink'
    container params.flye
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(filtered_fastq)
    val platform

    output:
    tuple val(sample_id), path("MT_assembly")

    script:
    def flye_preset = platform == 'ont' ? '--nano-hq' : 
                      platform == 'pb'  ? '--pacbio-hifi':
                      null
    """
    flye --threads ${task.cpus} --meta ${flye_preset} ${filtered_fastq} --out-dir MT_assembly -m ${params.flye_min_overlap}
    """
}

process ROTATE_ASSEMBLY {

    publishDir "${params.outdir}/${sample_id}/MT_assembly", pattern: "*.fasta", mode: 'symlink'
    publishDir path: "${params.outdir}/logs", pattern: "*.log", mode: 'symlink'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly_fasta), path(assembly_fasta_fai)
    tuple val(sample_id), path(assembly_align_ref_bam), path(assembly_align_ref_bam_index)

    output:
    tuple val(sample_id), path("assembly_rotated.fasta"), emit: fastq, optional: true
    path("rotate_assembly.log")

    script:
    """
    rotate.py --bam ${assembly_align_ref_bam} --assembly ${assembly_fasta}
    """
}

process INDEX_ASSEMBLY {

    publishDir "${params.outdir}/${sample_id}/MT_assembly", mode: 'symlink'
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


// process CHECK_CIRCULAR_GENOME{

//     publishDir "${params.outdir}/MT_assembly/", mode: 'symlink'
//     tag params.sample_id

//     input:
//     path mt_assembly_dir
//     val sample_id
//     tuple path(mt_assembly_to_ref_bam), path(mt_assembly_to_ref_bam_index)

//     output:
//     path("sieved_graph"), emit: sieved_graph_dir

//     script:
//     """
//     # check assembly graph for circular genome..
//     mkdir -p sieved_graph

//     ecLegov2.pl sievegraph \
//     --gv ${mt_assembly_dir}/assembly_graph.gv \
//     --diploidcov 24 \
//     --oprefix sieved_graph/${sample_id}.MT.assembly \
//     --bam ${mt_assembly_to_ref_bam}
//     """
// }

// process SUBPOPULATIONS {

//     publishDir "${params.outdir}/MT_assembly/", mode: 'symlink'

//     tag params.sample_id

//     input:
//     path mt_assembly_dir
//     val sample_id
//     mt_assembly_align_to_ref_bam

//     output:
//     path("MT_assembly/*"), emit: mt_assembly_dir


//     script:
//     """
//     echo '==' $(date) '==' subpopulation script setup..
//     cgSupPop.pl setup \
//     --assembly ${sample_id}.MT.assembly.cn24.disjointcyclic.gv.overview.xls \
//     --vcf ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf \
//     --dicncov 4 --subgraph 1 --rsg 1 --sample ${sample_id} \
//     --scriptPath .
//     #
//     mv process.sh process.sh.tmp
//     cat process.sh.tmp | sed s?t2tv2.fasta?MT.fasta? > process.sh
//     chmod u+x *.sh
//     echo '==' $(date) '==' writing possible subpopulations..
//     ./setup.sh
//     echo '==' $(date) '==' variation calling on possible subpopulations..
//     ./process.sh
//     echo '==' $(date) '==' refining alignments possible subpopulations..
//     ./refine.sh
//     """
// }