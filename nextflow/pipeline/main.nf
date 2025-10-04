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

include { ALIGNED_BAM_TO_FASTQ; ALIGNED_CRAM_TO_FASTQ; UNALIGNED_BAM_TO_FASTQ; UNALIGNED_CRAM_TO_FASTQ; COMPRESS_FASTQ } from './modules/format_inputs.nf'
include { KMER_SELECTION } from './modules/kmer_selection.nf'
include { ALIGN_TO_REF; 
          ALIGN_TO_ASSEMBLY as ALIGN_TO_ASSEMBLY; 
          ALIGN_TO_ASSEMBLY as ALIGN_TO_ROTATED_ASSEMBLY; 
          ALIGN_ASSEMBLY_TO_REF as ALIGN_ASSEMBLY_TO_REF;
          ALIGN_ASSEMBLY_TO_REF as ALIGN_ROTATED_ASSEMBLY_TO_REF} from './modules/alignment.nf'
include { FILTER_NUMTS; FILTERED_BAM_TO_FASTQ; } from './modules/filtration.nf'
include { MT_ASSEMBLY; ROTATE_ASSEMBLY; 
          INDEX_ASSEMBLY as INDEX_ASSEMBLY;
          INDEX_ASSEMBLY as INDEX_ROTATED_ASSEMBLY} from './modules/assembly.nf'
include { METH_FREQ; METH_PLOT} from './modules/methylation.nf'
include { MT_COVERAGE; MT_READ_LENGTH; COVERAGE_PLOT; READ_LENGTH_PLOT; QC_SUMMARY; COMBINE_QC_SUMMARY} from './modules/qc.nf'
include { VARIANT_CALLS_BALDUR as VARIANT_CALLS_BALDUR;
          VARIANT_CALLS_BALDUR as VARIANT_CALLS_BALDUR_ASSEMBLY; 
          NORMALIZE_BALDUR_VCF as NORMALIZE_BALDUR_VCF;
          NORMALIZE_BALDUR_VCF as NORMALIZE_BALDUR_VCF_ASSEMBLY;
          ANNOTATE_BALDUR_INDELS as ANNOTATE_BALDUR_INDELS;
          ANNOTATE_BALDUR_INDELS as ANNOTATE_BALDUR_INDELS_ASSEMBLY;
          ANNOTATE_BALDUR_SNVS as ANNOTATE_BALDUR_SNVS;
          ANNOTATE_BALDUR_SNVS as ANNOTATE_BALDUR_SNVS_ASSEMBLY} from './modules/variant_calling.nf'
include { VARIANT_CALLS_MUTSERVE as VARIANT_CALLS_MUTSERVE;
          VARIANT_CALLS_MUTSERVE as VARIANT_CALLS_MUTSERVE_ASSEMBLY;
          NORMALIZE_MUTSERVE_VCF as NORMALIZE_MUTSERVE_VCF;
          NORMALIZE_MUTSERVE_VCF as NORMALIZE_MUTSERVE_VCF_ASSEMBLY;
          ANNOTATE_MUTSERVE_VCF as ANNOTATE_MUTSERVE_VCF; 
          ANNOTATE_MUTSERVE_VCF as ANNOTATE_MUTSERVE_VCF_ASSEMBLY} from './modules/variant_calling.nf'
include { VARIANT_CALLS_SNIFFLES as VARIANT_CALLS_SNIFFLES; 
          VARIANT_CALLS_SNIFFLES as VARIANT_CALLS_SNIFFLES_ASSEMBLY;
          VARIANT_CALLS_SNIFFLES as VARIANT_CALLS_SNIFFLES_ASSEMBLY_TO_REF;
          FILTER_SNIFFLES_VCF_MINSUPPORT as FILTER_SNIFFLES_VCF_MINSUPPORT;
          FILTER_SNIFFLES_VCF_MINSUPPORT as FILTER_SNIFFLES_VCF_MINSUPPORT_ASSEMBLY;
          FILTER_SNIFFLES_VCF_MINSUPPORT as FILTER_SNIFFLES_VCF_MINSUPPORT_ASSEMBLY_TO_REF} from './modules/variant_calling.nf'
include { HAPLOGREP; HAPLOCHECK} from './modules/haplo.nf'

workflow {

    // Check for required parameters
    if (!params.samplesheet || !params.input_type || !params.is_aligned || !params.outdir || !params.platform) {
        error "Missing one of required parameters: samplesheet, input_type, is_aligned, outdir, platform"
    }

    // Set up channels for reusable parameters
    Channel.fromPath(params.kmc_pre, checkIfExists: true)
        .first()
        .set { kmc_pre_ch }
    
    Channel.fromPath(params.kmc_suf, checkIfExists: true)
        .first()
        .set { kmc_suf_ch }
    
    Channel.fromPath(params.mt_ref)
        .map { ref -> tuple(ref, file("${ref}.fai")) }
        .first()
        .set { mt_ref_ch }

    Channel.fromPath(params.mitomap_anno_file)
        .first()
        .set { mitomap_anno_file_ch }
    
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
        .fromPath(params.samplesheet) // Path to your samplesheet
        .splitCsv(header: true) // Split into a channel of maps, using header
        .map { row -> 
            def index_file  = null
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

    // Generate gzipped fastq file from input 
    if (params.is_aligned) {
        if (params.input_type == 'bam') {
            fastq_out = ALIGNED_BAM_TO_FASTQ(samples_ch)
            fastq_gz_out = COMPRESS_FASTQ(fastq_out)
        } else if (params.input_type == 'cram') {
            fastq_out = ALIGNED_CRAM_TO_FASTQ(samples_ch, params.reference)
            fastq_gz_out = COMPRESS_FASTQ(fastq_out)
        } 
    } else {
        if (params.input_type == 'bam') {
            fastq_out = UNALIGNED_BAM_TO_FASTQ(samples_ch)
            fastq_gz_out = COMPRESS_FASTQ(fastq_out)
        } else if (params.input_type == 'cram') {
            fastq_out = UNALIGNED_CRAM_TO_FASTQ(samples_ch, params.reference)
            fastq_gz_out = COMPRESS_FASTQ(fastq_out)
        } else if (params.input_type == 'fastq') {
            fastq_gz_out = COMPRESS_FASTQ(samples_ch)
        } else if (params.input_type == 'fastq.gz') {
            fastq_gz_out = samples_ch
        } 
    }

    // Select MT reads via kmer selection
    KMER_SELECTION(fastq_gz_out, kmc_pre_ch, kmc_suf_ch)
    
    // Generate filtered bam without NUMTs and discarded NUMT bam
    ALIGN_TO_REF(KMER_SELECTION.out, params.platform, minimap_index_ch)
    FILTER_NUMTS(ALIGN_TO_REF.out.bam)
    FILTERED_BAM_TO_FASTQ(FILTER_NUMTS.out.filtered_bam)

    // Assemble mito
    MT_ASSEMBLY(FILTERED_BAM_TO_FASTQ.out, params.platform)
    INDEX_ASSEMBLY(MT_ASSEMBLY.out.assembly_dir
        .map { sample_id, dir -> tuple(sample_id, file("${dir}/assembly.fasta"))})
        .set { assembly_fasta }

    // Align reads to assembly and align assembly to ref
    ALIGN_TO_ASSEMBLY(FILTERED_BAM_TO_FASTQ.out, assembly_fasta, params.platform)
    ALIGN_ASSEMBLY_TO_REF(assembly_fasta, minimap_index_ch, params.platform)

    // Methylation
    METH_FREQ(FILTER_NUMTS.out.filtered_bam, mt_ref_ch)
    METH_PLOT(METH_FREQ.out.minimod_tsv)

    // QC - coverage and read length
    MT_COVERAGE(FILTER_NUMTS.out.filtered_bam)
    COVERAGE_PLOT(MT_COVERAGE.out.per_base_bed)
    MT_READ_LENGTH(FILTER_NUMTS.out.filtered_bam)
    READ_LENGTH_PLOT(MT_READ_LENGTH.out)

    // QC summary 
    QC_SUMMARY(MT_COVERAGE.out.mosdepth_summary, MT_READ_LENGTH.out, METH_FREQ.out.minimod_tsv)
    COMBINE_QC_SUMMARY(QC_SUMMARY.out.collect())

    // SNV/indel/del variant calling (baldur) on reference
    VARIANT_CALLS_BALDUR(FILTER_NUMTS.out.filtered_bam, mt_ref_ch)
    NORMALIZE_BALDUR_VCF(VARIANT_CALLS_BALDUR.out.vcf)
    ANNOTATE_BALDUR_INDELS(NORMALIZE_BALDUR_VCF.out.norm_indels_vcf, mitomap_anno_file_ch)
    ANNOTATE_BALDUR_SNVS(NORMALIZE_BALDUR_VCF.out.norm_snvs_vcf, mitomap_anno_file_ch)

    // SNV variant calling (mutserve) on reference
    VARIANT_CALLS_MUTSERVE(FILTER_NUMTS.out.filtered_bam, mt_ref_ch, Channel.value("MT"))
    NORMALIZE_MUTSERVE_VCF(VARIANT_CALLS_MUTSERVE.out.vcf)
    ANNOTATE_MUTSERVE_VCF(NORMALIZE_MUTSERVE_VCF.out.norm_vcf, mitomap_anno_file_ch)

    // Haplogrep/haplocheck
    HAPLOGREP(VARIANT_CALLS_MUTSERVE.out.vcf)
    HAPLOCHECK(VARIANT_CALLS_MUTSERVE.out.vcf)

    // SV Calling on reference
    VARIANT_CALLS_SNIFFLES(FILTER_NUMTS.out.filtered_bam)
    FILTER_SNIFFLES_VCF_MINSUPPORT(VARIANT_CALLS_SNIFFLES.out.vcf)
    
    // SV Calling of assembly to reference
    VARIANT_CALLS_SNIFFLES_ASSEMBLY_TO_REF(ALIGN_ASSEMBLY_TO_REF.out.bam)

    // Rotate assembly to match reference coordinates and realign MT reads to it
    //ROTATE_ASSEMBLY(assembly_fasta, ALIGN_ASSEMBLY_TO_REF.out.bam, MT_ASSEMBLY.out) 
    //INDEX_ROTATED_ASSEMBLY(ROTATE_ASSEMBLY.out).set{ rotated_assembly_fasta }
    
    // // Align reads to rotated assembly and align rotated assembly to ref
    // ALIGN_TO_ROTATED_ASSEMBLY(FILTERED_BAM_TO_FASTQ.out, rotated_assembly_fasta, params.platform).set{ mt_align_rotated_assembly_bam }
    // ALIGN_ROTATED_ASSEMBLY_TO_REF(rotated_assembly_fasta, minimap_index_ch, params.platform)

    // // variant calling (baldur) on rotated assembly
    // VARIANT_CALLS_BALDUR_ASSEMBLY(mt_align_rotated_assembly_bam, rotated_assembly_fasta, sample_id_ch)
    // NORMALIZE_BALDUR_VCF_ASSEMBLY(VARIANT_CALLS_BALDUR_ASSEMBLY.out.baldur_vcf)
    // ANNOTATE_BALDUR_INDELS_ASSEMBLY(NORMALIZE_BALDUR_VCF_ASSEMBLY.out.baldur_norm_indels_vcf, mitomap_anno_file_ch)
    // ANNOTATE_BALDUR_SNVS_ASSEMBLY(NORMALIZE_BALDUR_VCF_ASSEMBLY.out.baldur_norm_snvs_vcf, mitomap_anno_file_ch)

    // // SNV variant calling (mutserve) on rotated assembly
    // VARIANT_CALLS_MUTSERVE_ASSEMBLY(mt_align_rotated_assembly_bam, rotated_assembly_fasta, sample_id_ch, Channel.value("contig_1_rotated"))
    // NORMALIZE_MUTSERVE_VCF_ASSEMBLY(VARIANT_CALLS_MUTSERVE_ASSEMBLY.out.mutserve_vcf)
    // ANNOTATE_MUTSERVE_VCF_ASSEMBLY(NORMALIZE_MUTSERVE_VCF_ASSEMBLY.out.mutserve_norm_vcf, mitomap_anno_file_ch)

    // // SV Calling on rotated assembly
    // VARIANT_CALLS_SNIFFLES_ASSEMBLY(mt_align_rotated_assembly_bam)
    // FILTER_SNIFFLES_VCF_MINSUPPORT_ASSEMBLY(VARIANT_CALLS_SNIFFLES_ASSEMBLY.out.sniffles_vcf)


}


