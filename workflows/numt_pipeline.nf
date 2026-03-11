include { NUMT_DETECTION_SNIFFLES; NUMT_DETECTION_INSERTIONS_TO_FASTA; NUMT_DETECTION_INSERTIONS_BLAST; 
          NUMT_DETECTION_SUPPLEMENTARY; NUMT_DETECTION_PLOT; NUMT_DETECTION_MTDNA_INSERTIONS_TO_FASTA; 
          NUMT_DETECTION_MTDNA_INSERTIONS_BLAST_CHECK } from '../modules/numt_detection.nf'

workflow NUMT_PIPELINE {

    take:
        samples

    main:

        Channel.fromPath(params.blast_db)
        .first()
        .set { blast_db_ch }

        NUMT_DETECTION_SNIFFLES(samples, params.reference)

        NUMT_DETECTION_INSERTIONS_TO_FASTA(NUMT_DETECTION_SNIFFLES.out.vcf)

        NUMT_DETECTION_INSERTIONS_BLAST(
            NUMT_DETECTION_INSERTIONS_TO_FASTA.out.fasta,
            blast_db_ch
        )

        NUMT_DETECTION_MTDNA_INSERTIONS_TO_FASTA(
            NUMT_DETECTION_INSERTIONS_TO_FASTA.out.fasta
            .join(NUMT_DETECTION_INSERTIONS_BLAST.out)
        )

        NUMT_DETECTION_MTDNA_INSERTIONS_BLAST_CHECK(
            NUMT_DETECTION_MTDNA_INSERTIONS_TO_FASTA.out,
            blast_db_ch
        )

        NUMT_DETECTION_PLOT(
            NUMT_DETECTION_MTDNA_INSERTIONS_BLAST_CHECK.out,
            params.circos_bed
        )

        NUMT_DETECTION_SUPPLEMENTARY(samples, params.reference)
}