include { NUMT_DETECTION_SNIFFLES; NUMT_DETECTION_INSERTIONS_TO_FASTA; NUMT_DETECTION_INSERTIONS_BLAST; 
          NUMT_DETECTION_SUPPLEMENTARY; NUMT_DETECTION_PLOT; NUMT_DETECTION_MTDNA_INSERTIONS_TO_FASTA; 
          NUMT_DETECTION_MTDNA_INSERTIONS_BLAST_CHECK; MAKE_BLAST_DBS; MAKE_ROTATED_MT_REF } from '../modules/numt_detection.nf'

workflow NUMT_PIPELINE {

    take:
        samples

    main:

        MAKE_ROTATED_MT_REF(params.mt_ref)
        rotated_mt_ref_ch = MAKE_ROTATED_MT_REF.out

        MAKE_BLAST_DBS(params.reference, params.mt_ref, rotated_mt_ref_ch)
        blast_db_ch = MAKE_BLAST_DBS.out.blast_db

        NUMT_DETECTION_SNIFFLES(samples, params.reference)

        NUMT_DETECTION_INSERTIONS_TO_FASTA(NUMT_DETECTION_SNIFFLES.out.vcf)

        NUMT_DETECTION_INSERTIONS_BLAST(
            NUMT_DETECTION_INSERTIONS_TO_FASTA.out.fasta,
            blast_db_ch
        )

        NUMT_DETECTION_MTDNA_INSERTIONS_TO_FASTA(
            NUMT_DETECTION_INSERTIONS_TO_FASTA.out.fasta
            .join(NUMT_DETECTION_INSERTIONS_BLAST.out.mt_blast)
        )

        NUMT_DETECTION_MTDNA_INSERTIONS_BLAST_CHECK(
            NUMT_DETECTION_MTDNA_INSERTIONS_TO_FASTA.out,
            blast_db_ch
        )

        NUMT_DETECTION_PLOT(
            NUMT_DETECTION_MTDNA_INSERTIONS_BLAST_CHECK.out.blast_output,
            params.circos_bed
        )

        NUMT_DETECTION_SUPPLEMENTARY(samples, params.reference)
}