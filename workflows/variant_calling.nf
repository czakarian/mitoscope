include { VARIANT_CALLS_BALDUR; NORMALIZE_BALDUR_VCF; VEP_BALDUR_VCF; ADD_MITOMAP_TO_BALDUR_VCF; ADD_HELIX_TO_BALDUR_VCF; PLOT_HETEROPLASMY; MERGE_BALDUR_VCFS} from '../modules/variant_calling.nf'
include { READ_VARIANT_MATRIX } from '../modules/variant_calling.nf'
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

        Channel.fromPath(params.gnomad_vcf)
        .map { vcf -> tuple(vcf, file("${vcf}.tbi")) }
        .first()
        .set { gnomad_anno_ch }

        Channel.fromPath(params.mitomap_anno_file)
        .map { tsv -> tuple(tsv, file("${tsv}.tbi"))}
        .first()
        .set { mitomap_anno_file_ch }

        Channel.fromPath(params.helix_anno_file)
        .map { tsv -> tuple(tsv, file("${tsv}.tbi"))}
        .first()
        .set { helix_anno_file_ch }

        Channel.fromPath(params.vep_cache_dir)
        .first()
        .set { vep_cache_ch }

        Channel.fromPath(params.vep_plugin_dir)
        .first()
        .set { vep_plugin_ch }

        // SNV/indel/del variant calling (baldur)
        VARIANT_CALLS_BALDUR(filtered_bam, mt_ref_ch)
        NORMALIZE_BALDUR_VCF(VARIANT_CALLS_BALDUR.out.vcf)
        VEP_BALDUR_VCF(NORMALIZE_BALDUR_VCF.out.norm_vcf, mt_ref_ch, vep_cache_ch, vep_plugin_ch, gnomad_anno_ch)
        ADD_MITOMAP_TO_BALDUR_VCF(VEP_BALDUR_VCF.out.vep_vcf, mitomap_anno_file_ch)
        ADD_HELIX_TO_BALDUR_VCF(ADD_MITOMAP_TO_BALDUR_VCF.out.anno_vcf, helix_anno_file_ch)
        PLOT_HETEROPLASMY(ADD_HELIX_TO_BALDUR_VCF.out.anno_vcf)

        // Generated merged baldur VCF across all samples
        baldur_vcfs = ADD_HELIX_TO_BALDUR_VCF.out.anno_vcf.map { it[1]}.collect()
        baldur_indexes = ADD_HELIX_TO_BALDUR_VCF.out.anno_vcf.map { it[2]}.collect()
        MERGE_BALDUR_VCFS(baldur_vcfs, baldur_indexes)

        if (params.generate_variant_matrix) {
            READ_VARIANT_MATRIX(filtered_bam.join(NORMALIZE_BALDUR_VCF.out.norm_vcf), mt_ref_ch)
        }

        // SNV variant calling with mutserve for haplogrep/haplocheck
        VARIANT_CALLS_MUTSERVE(filtered_bam, mt_ref_ch)
        HAPLOGREP(VARIANT_CALLS_MUTSERVE.out.vcf)
        HAPLOCHECK(VARIANT_CALLS_MUTSERVE.out.vcf)

        // SV Calling on reference
        VARIANT_CALLS_SNIFFLES(filtered_bam, mt_ref_ch)
        FILTER_SNIFFLES_VCF(VARIANT_CALLS_SNIFFLES.out.vcf)
        
        // Methylation
        if (params.platform == "pb" && params.skip_methylation == false) {
            METH_FREQ_PB(filtered_bam, mt_ref_ch)
        } else if (params.platform == "ont" && params.skip_methylation == false) {
            METH_FREQ_ONT(filtered_bam, mt_ref_ch)
        }

    emit:
        haplogrep_output_ch = HAPLOGREP.out.haplogrep_txt
    

}