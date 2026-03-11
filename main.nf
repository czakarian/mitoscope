#!/usr/bin/env nextflow

/*
========================================================================================
    czakarian/mitoscope
========================================================================================
    Github : https://github.com/czakarian/mitoscope
    Author: Christina Zakarian
    ---------------------------
*/

nextflow.enable.dsl = 2

include { NUMT_PIPELINE } from './workflows/numt_pipeline.nf'
include { FASTQ_PREP } from './workflows/fastq_preparation.nf'
include { QUALITY_CONTROL } from './workflows/quality_control.nf'
include { VARIANT_CALLING } from './workflows/variant_calling.nf'
include { ASSEMBLY } from './workflows/assembly.nf'
include { MT_SELECTION } from './workflows/mt_selection.nf'

workflow {

    // Check for required parameters
    if (!params.samplesheet || !params.input_type || !params.outdir || !params.platform || !params.reference) {
        error "Missing one of required parameters: samplesheet, input_type, is_aligned, outdir, platform, reference"
    }
    
    Channel.fromPath(params.mt_ref)
        .map { ref -> tuple(ref, file("${ref}.fai")) }
        .first()
        .set { mt_ref_ch }
    
    if (params.platform == "pb") {
        Channel.fromPath(params.mt_mmi_hifi, checkIfExists: true)
            .first()
            .set { minimap_index_ch }
    } else if (params.platform == "ont") {
        Channel.fromPath(params.mt_mmi_ont, checkIfExists: true)
            .first()
            .set { minimap_index_ch }
    } else {
        error "Invalid value for --platform. Must be 'ont' or 'pb'."
    }

    // read in sample sheet
    Channel
        .fromPath(params.samplesheet) 
        .splitCsv(header: ['sample_id', 'sample_file'])
        .map { row -> 
            def index_file  = ""
            def sample_file = file(row.sample_file)
            // Only assign index file if aligned BAM or CRAM
            if (params.is_aligned) {
                if (sample_file.toString().endsWith(".bam")) {
                    index_file = file("${sample_file}.bai")
                } else if (sample_file.toString().endsWith(".cram")) {
                    index_file = file("${sample_file}.crai")
                }
            }
            tuple(row.sample_id, sample_file, index_file)
    }.set { samples_ch }

    // Standarize inputs to gzipped fastq 
    FASTQ_PREP(samples_ch)

    MT_SELECTION(
        FASTQ_PREP.out.fastq_gz_out_ch,
        minimap_index_ch
    )
 
    ASSEMBLY(MT_SELECTION.out.filtered_fastq_ch)

    VARIANT_CALLING(
        MT_SELECTION.out.filtered_bam_ch,
        mt_ref_ch,
    )

    QUALITY_CONTROL(
        samples_ch, 
        MT_SELECTION.out.filtered_bam_ch, 
        MT_SELECTION.out.filtered_fastq_ch, 
        VARIANT_CALLING.out.haplogrep_output_ch, 
        MT_SELECTION.out.kmer_read_count_ch, 
        FASTQ_PREP.out.chrM_read_count_ch
    )

    if (params.numt_profiling && params.is_aligned) {
        NUMT_PIPELINE(samples_ch)
    }

}