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

include { BAM_TO_FASTQ; CRAM_TO_FASTQ; COMPRESS_FASTQ } from './modules/format_inputs.nf'
include { KMER_SELECTION } from './modules/kmer_selection.nf'
include { ALIGN_TO_REF; ALIGN_TO_ASSEMBLY; 
          SAMTOOLS_SAM_TO_BAM as SAMTOOLS_SAM_TO_BAM; 
          SAMTOOLS_SAM_TO_BAM as SAMTOOLS_SAM_TO_BAM_ASSEMBLY;
          SAMTOOLS_INDEX as SAMTOOLS_INDEX; 
          SAMTOOLS_INDEX as SAMTOOLS_INDEX_FILTERED; 
          SAMTOOLS_INDEX as SAMTOOLS_INDEX_NUMT; 
          SAMTOOLS_INDEX as SAMTOOLS_INDEX_ASSEMBLY} from './modules/alignment.nf'
include { FILTER_NUMTS } from './modules/filtration.nf'
include { MT_ASSEMBLY; BAM_TO_FASTQ_FOR_ASSEMBLY} from './modules/assembly.nf'
include { METH_FREQ; METH_PLOT} from './modules/methylation.nf'
include { MT_COVERAGE; MT_READ_LENGTH; COVERAGE_PLOT; READ_LENGTH_PLOT} from './modules/qc.nf'
include { VARIANT_CALLS_BALDUR; NORMALIZE_BALDUR_VCF; ANNOTATE_BALDUR_INDELS; ANNOTATE_BALDUR_SNVS} from './modules/variant_calling.nf'
include { VARIANT_CALLS_MUTSERVE; NORMALIZE_MUTSERVE_VCF; ANNOTATE_MUTSERVE_VCF} from './modules/variant_calling.nf'
include { VARIANT_CALLS_SNIFFLES; FILTER_SNIFFLES_VCF_MINSUPPORT} from './modules/variant_calling.nf'
include { HAPLOGREP; HAPLOCHECK} from './modules/haplo.nf'

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
    mitomap_anno_file_ch = Channel.fromPath(params.mitomap_anno_file)

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
    mt_align_bam = ALIGN_TO_REF(mt_fastq, platform_ch, minimap_index_ch) | SAMTOOLS_SAM_TO_BAM | SAMTOOLS_INDEX
    
    FILTER_NUMTS(mt_align_bam)
    mt_filtered_bam = SAMTOOLS_INDEX_FILTERED(FILTER_NUMTS.out.mt_filtered_bam)
    mt_numt_bam = SAMTOOLS_INDEX_NUMT(FILTER_NUMTS.out.mt_numt_bam)

    // Assembly
    mt_filtered_fastq = BAM_TO_FASTQ_FOR_ASSEMBLY(mt_filtered_bam)
    assembly_dir = MT_ASSEMBLY(mt_filtered_fastq, platform_ch)
    ALIGN_TO_ASSEMBLY(mt_filtered_fastq, platform_ch, assembly_dir) | SAMTOOLS_SAM_TO_BAM_ASSEMBLY | SAMTOOLS_INDEX_ASSEMBLY

    // Methylation
    METH_FREQ(mt_filtered_bam, mt_ref_ch)
    METH_PLOT(METH_FREQ.out.minimod_tsv)

    // QC - coverage 
    MT_COVERAGE(mt_filtered_bam)
    COVERAGE_PLOT(MT_COVERAGE.out.per_base_bed)

    // QC - ead length dist
    MT_READ_LENGTH(mt_filtered_bam)
    READ_LENGTH_PLOT(MT_READ_LENGTH.out.read_length_file)

    // SNV/indel/del variant calling - baldur
    VARIANT_CALLS_BALDUR(mt_filtered_bam, mt_ref_ch, sample_id_ch)
    NORMALIZE_BALDUR_VCF(VARIANT_CALLS_BALDUR.out.baldur_vcf, mt_ref_ch)
    ANNOTATE_BALDUR_INDELS(NORMALIZE_BALDUR_VCF.out.baldur_norm_indels_vcf, mitomap_anno_file_ch)
    ANNOTATE_BALDUR_SNVS(NORMALIZE_BALDUR_VCF.out.baldur_norm_snvs_vcf, mitomap_anno_file_ch)

    // SNV variant calling - mutserve
    VARIANT_CALLS_MUTSERVE(mt_filtered_bam, mt_ref_ch, sample_id_ch)
    NORMALIZE_MUTSERVE_VCF(VARIANT_CALLS_MUTSERVE.out.mutserve_vcf, mt_ref_ch)
    ANNOTATE_MUTSERVE_VCF(NORMALIZE_MUTSERVE_VCF.out.mutserve_norm_vcf, mitomap_anno_file_ch)

    // Haplogrep/haplocheck
    HAPLOGREP(VARIANT_CALLS_MUTSERVE.out.mutserve_vcf)
    HAPLOCHECK(VARIANT_CALLS_MUTSERVE.out.mutserve_vcf)

    // SV Calling
    VARIANT_CALLS_SNIFFLES(mt_filtered_bam)
    FILTER_SNIFFLES_VCF_MINSUPPORT(VARIANT_CALLS_SNIFFLES.out.sniffles_vcf)


}


