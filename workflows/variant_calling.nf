include { VARIANT_CALLS_BALDUR; NORMALIZE_BALDUR_VCF; ANNOTATE_BALDUR; VEP_BALDUR_VCF; ADD_MITOMAP_TO_BALDUR_VCF; MERGE_BALDUR_VCFS} from '../modules/variant_calling.nf'
include { VARIANT_CALLS_MUTSERVE } from '../modules/variant_calling.nf'
include { VARIANT_CALLS_SNIFFLES; FILTER_SNIFFLES_VCF; COMBINE_SV_CALLS} from '../modules/variant_calling.nf'
include { HAPLOGREP; HAPLOCHECK} from '../modules/haplo.nf'
include { METH_FREQ_ONT; METH_FREQ_PB} from '../modules/methylation.nf'

workflow VARIANT_CALLING {


    take:
        filtered_bam

    main:

        Channel.fromPath(params.mt_ref)
        .map { ref -> tuple(ref, file("${ref}.fai")) }
        .first()
        .set { mt_ref_ch }

        Channel.fromPath(params.mitomap_anno_file)
        .first()
        .set { mitomap_anno_file_ch }

        Channel.fromPath(params.vep_cache_dir)
        .first()
        .set { vep_cache_ch }

        // SNV/indel/del variant calling (baldur)
        VARIANT_CALLS_BALDUR(filtered_bam, mt_ref_ch)
        NORMALIZE_BALDUR_VCF(VARIANT_CALLS_BALDUR.out.vcf)
        ANNOTATE_BALDUR(NORMALIZE_BALDUR_VCF.out.norm_vcf, mitomap_anno_file_ch)
        VEP_BALDUR_VCF(NORMALIZE_BALDUR_VCF.out.norm_vcf, mt_ref_ch, vep_cache_ch)
        ADD_MITOMAP_TO_BALDUR_VCF(VEP_BALDUR_VCF.out.vep_vcf.join(ANNOTATE_BALDUR.out.mitomap_txt))
        baldur_vcfs = ADD_MITOMAP_TO_BALDUR_VCF.out.anno_vcf.map { it[1]}.collect()
        baldur_indexes = ADD_MITOMAP_TO_BALDUR_VCF.out.anno_vcf.map { it[2]}.collect()
        MERGE_BALDUR_VCFS(baldur_vcfs, baldur_indexes)

        // SNV variant calling with mutserve for haplogrep/haplocheck
        VARIANT_CALLS_MUTSERVE(filtered_bam, mt_ref_ch)
        HAPLOGREP(VARIANT_CALLS_MUTSERVE.out.vcf)
        HAPLOCHECK(VARIANT_CALLS_MUTSERVE.out.vcf)

        // SV Calling on reference
        VARIANT_CALLS_SNIFFLES(filtered_bam, mt_ref_ch)
        FILTER_SNIFFLES_VCF(VARIANT_CALLS_SNIFFLES.out.vcf)
        
        // Methylation
        if (params.platform == "pb") {
            METH_FREQ_PB(filtered_bam, mt_ref_ch)
        } else if (params.platform == "ont") {
            METH_FREQ_ONT(filtered_bam, mt_ref_ch)
        }

    emit:
        haplogrep_output_ch = HAPLOGREP.out.haplogrep_txt
    

}