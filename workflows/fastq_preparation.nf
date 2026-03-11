
include { ALIGNED_BAM_TO_FASTQ; ALIGNED_CRAM_TO_FASTQ; UNALIGNED_BAM_TO_FASTQ; UNALIGNED_CRAM_TO_FASTQ; COMPRESS_FASTQ } from '../modules/format_inputs.nf'

workflow FASTQ_PREP {

    take:
        samples_ch

    main:
        // Generate gzipped fastq file from input 
        if (params.is_aligned) {
            if (params.input_type == 'bam') {
                ALIGNED_BAM_TO_FASTQ(samples_ch)
                chrM_read_count = ALIGNED_BAM_TO_FASTQ.out.chrM_read_count
                fastq_gz_out = COMPRESS_FASTQ(ALIGNED_BAM_TO_FASTQ.out.fastq)
            } else if (params.input_type == 'cram') {
                ALIGNED_CRAM_TO_FASTQ(samples_ch, params.reference)
                chrM_read_count = ALIGNED_CRAM_TO_FASTQ.out.chrM_read_count
                fastq_gz_out = COMPRESS_FASTQ(ALIGNED_CRAM_TO_FASTQ.out.fastq)
            } 

        } else {
            if (params.input_type == 'bam') {
                UNALIGNED_BAM_TO_FASTQ(samples_ch)
                fastq_gz_out = COMPRESS_FASTQ(UNALIGNED_BAM_TO_FASTQ.out.fastq)
            } else if (params.input_type == 'cram') {
                UNALIGNED_CRAM_TO_FASTQ(samples_ch, params.reference)
                fastq_gz_out = COMPRESS_FASTQ(UNALIGNED_CRAM_TO_FASTQ.out.fastq)
            } else if (params.input_type == 'fastq') {
                fastq_gz_out = COMPRESS_FASTQ(samples_ch)
            } else if (params.input_type == 'fastq.gz') {
                fastq_gz_out = samples_ch
            } 
        }

    emit:
        fastq_gz_out_ch = fastq_gz_out
        chrM_read_count_ch = chrM_read_count

}