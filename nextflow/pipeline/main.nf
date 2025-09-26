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

include { BAM_TO_FASTQ; CRAM_TO_FASTQ; COMPRESS_FASTQ } from './modules/bam_to_fastq.nf'
include { KMER_SELECTION } from './modules/kmer_selection.nf'
include { ALIGN_TO_REF; SAMTOOLS_SAM_TO_BAM; SAMTOOLS_SORT } from './modules/align_to_ref.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX; SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED; SAMTOOLS_INDEX as SAMTOOLS_INDEX_NUMT} from './modules/align_to_ref.nf'
include { FILTER_NUMTS } from './modules/filter_numts.nf'
include { METH_FREQ; METH_PLOT} from './modules/methylation.nf'
include { MT_COVERAGE; MT_READ_LENGTH; COVERAGE_PLOT; READ_LENGTH_PLOT} from './modules/qc.nf'


workflow {

    // Check for required parameters
    if (!params.inputfile || !params.sample_id || !params.outdir || !params.platform) {
        error "Missing required params. Usage: --inputfile --sample_id --outdir --platform"
    }
    if (!file(params.inputfile).exists()) error "Input file does not exist: ${params.inputfile}"

    def bam_ch = null
    def cram_ch = null
    def fastq_ch = null

    if (params.inputfile.endsWith('.bam')) {
        bam_ch = Channel.fromPath(params.inputfile)
    } else if (params.inputfile.endsWith('.cram')) {
        cram_ch = Channel.fromPath(params.inputfile)
    } else if (params.inputfile.endsWith('.fastq') || params.inputfile.endsWith('.fastq.gz')) {
        fastq_ch = Channel.fromPath(params.inputfile)
    } else {
        error "Unsupported input format: ${params.inputfile}, must be .bam, .cram, .fastq, or .fastq.gz"
    }

    // check if cram input file that ref file is provided 
    if (cram_ch) {
        if (!params.reference) {
            error "Missing genome reference fasta (--reference) to accompany .cram input"
        }
        else if (!file(params.reference).exists()) {
            error "Genome reference fasta file does not exist: ${params.reference}"
        }
        else {
            ref_ch = Channel.fromPath(params.reference)
        }
    }
    else {
        ref_ch = null
    }


    platform_ch = Channel.value(params.platform)
    sample_id_ch = Channel.value(params.sample_id)

    if (params.platform == "pb") {
        minimap_index_ch = Channel.fromPath(params.mt_mmi_hifi)
    } else if (params.platform == "ont") {
        minimap_index_ch = Channel.fromPath(params.mt_mmi_ont)
    } else {
        error "Invalid value for --platform. Must be 'ont' or 'pb'."
    }

    // Other static parameters
    mt_ref_ch = Channel.fromPath(params.mt_ref)
    kmc_pre_ch = Channel.fromPath(params.kmc_pre)
    kmc_suf_ch = Channel.fromPath(params.kmc_suf)

    // Workflow 
    def fastq_out
    def fastq_gz_out

    if (bam_ch) {
        fastq_out = BAM_TO_FASTQ(bam_ch, sample_id_ch)
        fastq_gz_out = COMPRESS_FASTQ(fastq_out)
    } else if (cram_ch) {
        fastq_out = CRAM_TO_FASTQ(cram_ch, sample_id_ch, ref_ch)
        fastq_gz_out = COMPRESS_FASTQ(fastq_out)
    } else if (params.inputfile.endsWith('.fastq')) {
        fastq_gz_out = COMPRESS_FASTQ(fastq_ch)
    } else {
        fastq_gz_out = fastq_ch
    }

    mt_fastq = KMER_SELECTION(fastq_gz_out, kmc_pre_ch, kmc_suf_ch)
    mt_align_bam = ALIGN_TO_REF(mt_fastq, platform_ch, minimap_index_ch) | SAMTOOLS_SAM_TO_BAM | SAMTOOLS_SORT | SAMTOOLS_INDEX
    
    FILTER_NUMTS(mt_align_bam)
    mt_filtered_bam = SAMTOOLS_INDEX_FILTERED(FILTER_NUMTS.out.mt_filtered_bam)
    mt_numt_bam = SAMTOOLS_INDEX_NUMT(FILTER_NUMTS.out.mt_numt_bam)

    // methylation
    METH_FREQ(mt_filtered_bam, mt_ref_ch)
    METH_PLOT(METH_FREQ.out.minimod_tsv)

    // QC - coverage 
    MT_COVERAGE(mt_filtered_bam)
    COVERAGE_PLOT(MT_COVERAGE.out.per_base_bed)

    // QC - ead length dist
    MT_READ_LENGTH(mt_filtered_bam)
    READ_LENGTH_PLOT(MT_READ_LENGTH.out.read_length_file)

}


