include { MT_COVERAGE; NUCLEAR_COVERAGE; MT_READ_LENGTH; COVERAGE_PLOT; READ_LENGTH_PLOT; 
          QC_SUMMARY; QC_SUMMARY_FOR_UNALIGNED_INPUT; COMBINE_QC_SUMMARY} from '../modules/qc.nf'

workflow QUALITY_CONTROL {

    take:
        samples
        filtered_bam
        filtered_fastq
        haplogrep_output
        kmer_read_counts
        chrM_read_counts

    main:

        MT_COVERAGE(filtered_bam)
        mito_coverage = MT_COVERAGE.out.mosdepth_summary
        COVERAGE_PLOT(MT_COVERAGE.out.per_base_bed)

        MT_READ_LENGTH(filtered_fastq)
        mito_read_lengths = MT_READ_LENGTH.out
        READ_LENGTH_PLOT(mito_read_lengths)

        if (params.check_nuclear_coverage && params.is_aligned) {
            NUCLEAR_COVERAGE(samples, params.reference, params.nuclear_intervals)
            nuc_coverage = NUCLEAR_COVERAGE.out.mosdepth_summary
        }

        // Generate summary output table for samples in run
        if (params.is_aligned && params.check_nuclear_coverage) {
            QC_SUMMARY(mito_coverage
                .join(mito_read_lengths)
                .join(nuc_coverage)
                .join(haplogrep_output)
                .join(kmer_read_counts)
                .join(chrM_read_counts)
            )
            COMBINE_QC_SUMMARY(QC_SUMMARY.out.collect())
        } else {
            QC_SUMMARY_FOR_UNALIGNED_INPUT(mito_coverage.
                join(mito_read_lengths).
                join(haplogrep_output).
                join(kmer_read_counts)
            )
            COMBINE_QC_SUMMARY(QC_SUMMARY_FOR_UNALIGNED_INPUT.out.collect())
        }
}