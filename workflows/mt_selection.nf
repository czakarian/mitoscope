include { KMER_SELECTION } from '../modules/kmer_selection.nf'
include { ALIGN_TO_REF } from '../modules/alignment.nf'
include { FILTER_NUMTS; FILTERED_BAM_TO_FASTQ; } from '../modules/filtration.nf'

workflow MT_SELECTION {

    take:
        fastq_gz
        minimap_index_ch

    main:

        // Set up channels for reusable parameters
        Channel.fromPath(params.kmc_pre, checkIfExists: true)
            .first()
            .set { kmc_pre_ch }
        
        Channel.fromPath(params.kmc_suf, checkIfExists: true)
            .first()
            .set { kmc_suf_ch }

        // Select MT reads via kmer selection
        KMER_SELECTION(fastq_gz, kmc_pre_ch, kmc_suf_ch)
        
        // Generate filtered bam without NUMTs and discarded NUMT bam
        ALIGN_TO_REF(KMER_SELECTION.out.fastq, params.platform, minimap_index_ch)
        FILTER_NUMTS(ALIGN_TO_REF.out.bam)
        FILTERED_BAM_TO_FASTQ(FILTER_NUMTS.out.filtered_bam)

    emit:
        kmer_read_count_ch = KMER_SELECTION.out.kmer_read_count
        filtered_bam_ch = FILTER_NUMTS.out.filtered_bam
        filtered_fastq_ch = FILTERED_BAM_TO_FASTQ.out


}