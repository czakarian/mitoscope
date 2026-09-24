process VARIANT_CALLS_BALDUR {

    publishDir "${params.outdir}/${sample_id}/variants/", pattern: "*.baldur_del.txt", mode: 'copy'
    container params.baldur
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)
    tuple path(mt_ref), path(mt_ref_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.baldur.vcf.gz"), emit: vcf
    tuple val(sample_id), path("${input_bam.getBaseName()}.baldur_del.txt"), emit: dels, optional: true

    script:
    """
    set -euo pipefail

    baldur -l debug \
    --output-deletions \
    --small-deletion-limit ${params.deletion_size_limit} \
    --large-deletion-limit ${params.deletion_size_limit} \
    --indel-thresholds ${params.indel_threshold} ${params.indel_threshold} \
    --snv-thresholds ${params.snv_hard_threshold} ${params.snv_soft_threshold}  \
    --homopolymer-limit ${params.homopolymer_limit} \
    -q${params.mapq_threshold} -Q${params.base_qual_threshold} -I${params.max_indel_base_qual} \
    -T ${mt_ref} \
    -n ${sample_id} \
    -o "${input_bam.getBaseName()}.baldur" \
    ${input_bam} 

    """
}

process NORMALIZE_BALDUR_VCF {
    //publishDir "${params.outdir}/${sample_id}/variants/", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_vcf)

    output:
    tuple val(sample_id), path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz.tbi"), emit: norm_vcf

    script:
    """
    set -euo pipefail

    bcftools norm --no-version --multiallelics -both ${baldur_vcf} \
    | bcftools norm --no-version --atomize --atom-overlaps . \
    | bcftools annotate -x INFO/DP,INFO/INDEL \
    | bcftools view --no-version -i 'FMT/HPL>=${params.snv_soft_threshold}' -f PASS -Oz -o ${baldur_vcf.getBaseName(2)}.norm.vcf.gz 
    
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.vcf.gz  

    """
}

process VEP_BALDUR_VCF {
    //publishDir "${params.outdir}/${sample_id}/variants/", mode: 'copy'
    container params.vep
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_norm_vcf), path(baldur_norm_vcf_index)
    tuple path(mt_ref), path(mt_ref_index)
    path vep_cache_dir
    path vep_plugin_dir
    tuple path(gnomad_vcf), path(gnomad_vcf_index)

    output:
    tuple val(sample_id), path("${baldur_norm_vcf.getBaseName(2)}.vep.vcf.gz"), path("${baldur_norm_vcf.getBaseName(2)}.vep.vcf.gz.tbi"), emit: vep_vcf

    script:
    """
    set -euo pipefail

    vep --input_file ${baldur_norm_vcf} \
    --safe \
    --cache \
    --dir_cache ${vep_cache_dir} \
    --offline \
    --fork ${task.cpus} \
    --fasta ${mt_ref} \
    --vcf \
    --compress_output bgzip \
    --output_file ${baldur_norm_vcf.getBaseName(2)}.vep.vcf.gz \
    --hgvs \
    --protein \
    --symbol \
    --biotype \
    --no_stats \
    --allow_non_variant \
    --distance 0 \
    --dir_plugins ${vep_plugin_dir} \
    --plugin gnomADMt,file=${gnomad_vcf},fields=AN:AC_hom:AC_het:AF_hom:AF_het:max_hl

    tabix -p vcf ${baldur_norm_vcf.getBaseName(2)}.vep.vcf.gz

    """

}

process ADD_MITOMAP_TO_BALDUR_VCF {
    //publishDir "${params.outdir}/${sample_id}/variants/", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_index)
    path mitomap_annos

    output:
    tuple val(sample_id), path("${sample_id}.mt.baldur.mitomap.vcf.gz"), path("${sample_id}.mt.baldur.mitomap.vcf.gz.tbi"), emit: anno_vcf

    script:
    """

    echo '##INFO=<ID=MITOMAP_Disease,Number=A,Type=String,Description="MITOMAP associated diseases.">
    ##INFO=<ID=MITOMAP_Disease_Status,Number=A,Type=String,Description="MITOMAP disease status.">
    ##INFO=<ID=MITOMAP_GenBank_AC,Number=A,Type=Float,Description="MITOMAP GenBank sample count.">
    ##INFO=<ID=MITOMAP_GenBank_AF,Number=A,Type=Float,Description="MITOMAP GenBank population allele frequency.">' > header.txt

    bcftools annotate --no-version ${input_vcf} \
        -c CHROM,POS,REF,ALT,INFO/MITOMAP_Disease,INFO/MITOMAP_Disease_Status,INFO/MITOMAP_GenBank_AC,INFO/MITOMAP_GenBank_AF \
        --annotations ${mitomap_annos} \
        --header-lines header.txt \
        --threads ${task.cpus} \
        --output-type z \
        --write-index=tbi \
        --output ${sample_id}.mt.baldur.mitomap.vcf.gz 

    """

}


process ADD_HELIX_TO_BALDUR_VCF {
    publishDir "${params.outdir}/${sample_id}/variants/", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_index)
    path helix_annos

    output:
    tuple val(sample_id), path("${sample_id}.mt.baldur.annotated.vcf.gz"), path("${sample_id}.mt.baldur.annotated.vcf.gz.tbi"), emit: anno_vcf

    script:
    """

    echo '##INFO=<ID=Helix_Counts_Hom,Number=A,Type=String,Description="HelixMTdb homoplasmic population sample count.">
    ##INFO=<ID=Helix_AF_Hom,Number=A,Type=String,Description="HelixMTdb homoplasmic population allele frequency.">
    ##INFO=<ID=Helix_Counts_Het,Number=A,Type=String,Description="HelixMTdb heteroplasmic population sample count.">
    ##INFO=<ID=Helix_AF_Het,Number=A,Type=String,Description="HelixMTdb heteroplasmic population allele frequency.">
    ##INFO=<ID=Helix_max_ARF,Number=A,Type=String,Description="HelixMTdb maximum heteroplasmy fraction across heteroplasmic samples.">
    ##INFO=<ID=Helix_mean_ARF,Number=A,Type=String,Description="HelixMTdb average heteroplasmy fraction across heteroplasmic samples.">' > header.txt

    bcftools annotate --no-version ${input_vcf} \
        -c CHROM,POS,REF,ALT,INFO/Helix_Counts_Hom,INFO/Helix_AF_Hom,INFO/Helix_Counts_Het,INFO/Helix_AF_Het,INFO/Helix_max_ARF,INFO/Helix_mean_ARF \
        --annotations ${helix_annos} \
        --header-lines header.txt \
        --threads ${task.cpus} \
        --output-type z \
        --write-index=tbi \
        --output ${sample_id}.mt.baldur.annotated.vcf.gz 

    """

}


process PLOT_HETEROPLASMY {
    
    publishDir "${params.outdir}/${sample_id}/variants/", pattern: "*.heteroplasmy.png", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(anno_vcf), path(anno_vcf_index)

    output:
    path("${anno_vcf.getBaseName(3)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail

    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    heteroplasmy_plot.py \
    --input ${anno_vcf} \
    --output_prefix ${anno_vcf.getBaseName(3)}
    """
}

process MERGE_BALDUR_VCFS {
    publishDir "${params.outdir}/", mode: 'copy'
    container params.bcftools

    input:
    path baldur_vcfs
    path baldur_indexes

    output:
    tuple path("merged.mt.baldur.annotated.vcf.gz"), path("merged.mt.baldur.annotated.vcf.gz.tbi")

    script:
    """
    set -euo pipefail

    bcftools merge --force-single -m none ${baldur_vcfs} -Oz -o merged.mt.baldur.annotated.vcf.gz
    bcftools index --tbi merged.mt.baldur.annotated.vcf.gz

    """
}

process VARIANT_CALLS_MUTSERVE {

    // publishDir "${params.outdir}/${sample_id}/variants/mutserve/", mode: 'copy'
    container params.mutserve
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)
    tuple path(mt_ref), path(mt_ref_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.mutserve.vcf.gz"), emit: vcf

    script:
    """
    set -euo pipefail

    mutserve call ${input_bam} \
    --output ${input_bam.getBaseName()}.mutserve.vcf.gz \
    --reference ${mt_ref} \
    --contig-name MT \
    --threads ${task.cpus} --no-ansi \
    --alignQ 0 --mapQ 20 --baseQ 20 --level 0.005

    """
}

process NORMALIZE_MUTSERVE_VCF {
    
    publishDir "${params.outdir}/${sample_id}/variants/mutserve/", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mutserve_vcf)

    output:
    tuple val(sample_id), 
        path("${mutserve_vcf.getBaseName(2)}.norm.vcf.gz"), 
        path("${mutserve_vcf.getBaseName(2)}.norm.vcf.gz.tbi"), emit: norm_vcf

    script:
    """
    set -euo pipefail

    bcftools norm --multiallelics -both ${mutserve_vcf} | bcftools norm --atomize | bcftools view -f PASS -Oz -o ${mutserve_vcf.getBaseName(2)}.norm.vcf.gz
    bcftools index --tbi ${mutserve_vcf.getBaseName(2)}.norm.vcf.gz
    """

}

process ANNOTATE_MUTSERVE_VCF {
    
    publishDir "${params.outdir}/${sample_id}/variants/mutserve/", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mutserve_norm_vcf), path(mutserve_norm_vcf_index)
    path mitomap_anno_file

    output:
    path("${mutserve_norm_vcf.getBaseName(2)}.annotated.txt")
    path("${mutserve_norm_vcf.getBaseName(2)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail

    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    annotate.py \
    --input ${mutserve_norm_vcf} \
    --annotations ${mitomap_anno_file} \
    --caller mutserve 
    """


}

process VARIANT_CALLS_SNIFFLES {

    //publishDir "${params.outdir}/${sample_id}/variants/sniffles/", mode: 'copy'
    //publishDir "${params.outdir}/${sample_id}/variants/sniffles/assembly_to_ref", pattern: "*assembly.ref.sniffles*", mode: 'copy'
    container params.sniffles
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)
    tuple path(mt_ref), path(mt_ref_index)

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.sniffles.vcf"), emit: vcf
    tuple val(sample_id), path("${input_bam.getBaseName()}.sniffles.snf"), emit: snf

    script:
    """
    set -euo pipefail

    sniffles --qc-output-all \
    --minsvlen 5 \
    --minsupport ${params.min_sv_support} \
    --threads ${task.cpus} \
    --ref ${mt_ref} \
    --input ${input_bam} \
    --snf ${input_bam.getBaseName()}.sniffles.snf \
    --vcf ${input_bam.getBaseName()}.sniffles.vcf

    """
}


process FILTER_SNIFFLES_VCF {
    
    publishDir "${params.outdir}/${sample_id}/variants/", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(sniffles_vcf)

    output:
    tuple val(sample_id), path("${sniffles_vcf.getBaseName()}.ge${params.min_sv_support}.svlen5.vcf")
    tuple val(sample_id), path("${sniffles_vcf.getBaseName()}.ge${params.min_sv_support}.svlen50.vcf")


    script:
    """
    set -euo pipefail

    bcftools filter -i "SUPPORT>=${params.min_sv_support} && (SVLEN>=50 || SVLEN<=-50)" ${sniffles_vcf} > ${sniffles_vcf.getBaseName()}.ge${params.min_sv_support}.svlen50.vcf
    bcftools filter -i "SUPPORT>=${params.min_sv_support} && ((SVLEN>=5 && SVLEN<50) || (SVLEN<=-5 && SVLEN>-50))" ${sniffles_vcf} > ${sniffles_vcf.getBaseName()}.ge${params.min_sv_support}.svlen5.vcf

    """
}

process COMBINE_SV_CALLS {

    publishDir "${params.outdir}/", mode: 'copy'
    container params.sniffles

    input:
    path snf_files

    output:
    path "multisample.sniffles.vcf"

    script:
    """
    sniffles \
    --qc-output-all \
    --minsvlen 5 \
    --combine-low-confidence 0 \
    --combine-low-confidence-abs 0 \
    --combine-null-min-coverage 1 \
    --combine-output-filtered \
    --combine-pair-relabel \
    --input ${snf_files} \
    --vcf multisample.sniffles.vcf
    """

}


process READ_VARIANT_MATRIX {
    
    publishDir "${params.outdir}/${sample_id}/variants/read_variant_matrix", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index), path(baldur_norm_vcf), path(baldur_norm_vcf_index)
    tuple path(mt_ref), path(mt_ref_index)

    output:
    path("${input_bam.getBaseName()}.snv_indel.csv")
    path("${input_bam.getBaseName()}.deletions.csv")

    script:
    """
    set -euo pipefail

    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    make_read_variant_matrix.py \
    --bam ${input_bam} \
    --vcf ${baldur_norm_vcf} \
    --output ${input_bam.getBaseName()} \
    --ref ${mt_ref}

    """
}