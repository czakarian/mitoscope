
include { MT_ASSEMBLY; ROTATE_ASSEMBLY; INDEX_ASSEMBLY } from '../modules/assembly.nf'

workflow ASSEMBLY {

    take:
        filtered_fastq

    main:

        // Assemble mito
        MT_ASSEMBLY(filtered_fastq, params.platform)

        // if need to index assembly, do within MT_assembly process (add tabix)
        // INDEX_ASSEMBLY(MT_ASSEMBLY.out.assembly_dir
        //     .map { sample_id, dir -> tuple(sample_id, file("${dir}/assembly.fasta"))})
        //     .set { assembly_fasta }

        // Align reads to assembly and align assembly to ref
        // ALIGN_TO_ASSEMBLY(FILTERED_BAM_TO_FASTQ.out.join(assembly_fasta), params.platform)
        // ALIGN_ASSEMBLY_TO_REF(assembly_fasta, minimap_index_ch, params.platform)


}