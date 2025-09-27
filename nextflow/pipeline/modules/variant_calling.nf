process VARIANT_CALLS_BALDUR {

    publishDir "${params.outdir}/variants/baldur", mode: 'symlink'
    container params.baldur
    tag params.sample_id

    input:
    tuple path(input_bam), path(input_bam_index)
    path mt_ref
    val sample_id

    output:
    path("${input_bam.getBaseName()}.baldur.vcf.gz"), emit: baldur_vcf
    path("${input_bam.getBaseName()}.baldur_del.txt"), emit: baldur_del_file

    script:
    """
    set -euo pipefail

    # baldur variant calls (SNV, small indel, large deletion)
    baldur -n ${sample_id} -l debug --output-deletions -T ${mt_ref} -o "${input_bam.getBaseName()}.baldur" ${input_bam}

    """
}

process NORMALIZE_BALDUR_VCF {

    publishDir "${params.outdir}/variants/baldur", mode: 'symlink'
    container params.bcftools
    tag params.sample_id

    input:
    path baldur_vcf
    path mt_ref

    output:
    tuple path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz.tbi"), emit: baldur_norm_vcf
    tuple path("${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz.tbi"), emit: baldur_norm_indels_vcf
    tuple path("${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz.tbi"), emit: baldur_norm_snvs_vcf

    script:
    """
    set -euo pipefail

    bcftools norm --multiallelics -both ${baldur_vcf} | bcftools norm --atomize | bcftools view -f PASS -Oz -o ${baldur_vcf.getBaseName(2)}.norm.vcf.gz 
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.vcf.gz  

    bcftools view --types indels ${baldur_vcf.getBaseName(2)}.norm.vcf.gz -Oz -o ${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz 
    bcftools view --exclude-types indels ${baldur_vcf.getBaseName(2)}.norm.vcf.gz -Oz -o ${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz 
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz 
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz 

    """
}

process ANNOTATE_BALDUR_SNVS {
    
    publishDir "${params.outdir}/variants/baldur", mode: 'symlink'
    container params.mitoscope
    tag params.sample_id

    input:
    tuple path(baldur_norm_snvs_vcf),path(baldur_norm_snvs_vcf_index)
    path mitomap_anno_file

    output:
    path("${baldur_norm_snvs_vcf.getBaseName(2)}.annotated.txt")
    path("${baldur_norm_snvs_vcf.getBaseName(2)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail

    annotate.py \
    --input ${baldur_norm_snvs_vcf} \
    --annotations ${mitomap_anno_file} \
    --caller baldur 
    """
}

process ANNOTATE_BALDUR_INDELS {
    
    publishDir "${params.outdir}/variants/baldur", mode: 'symlink'
    container params.mitoscope
    tag params.sample_id

    input:
    tuple path(baldur_norm_indels_vcf), path(baldur_norm_indels_vcf_index)
    path mitomap_anno_file

    output:
    path("${baldur_norm_indels_vcf.getBaseName(2)}.annotated.txt")
    path("${baldur_norm_indels_vcf.getBaseName(2)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail

    annotate.py \
    --input ${baldur_norm_indels_vcf} \
    --annotations ${mitomap_anno_file} \
    --caller baldur 
    """
}

process VARIANT_CALLS_MUTSERVE {

    publishDir "${params.outdir}/variants/mutserve", mode: 'symlink'
    container params.mutserve
    tag params.sample_id

    input:
    tuple path(input_bam), path(input_bam_index)
    path mt_ref
    val sample_id

    output:
    path("${input_bam.getBaseName()}.mutserve.vcf.gz"), emit: mutserve_vcf

    script:
    """
    set -euo pipefail

    mutserve call ${input_bam} \
    --output ${input_bam.getBaseName()}.mutserve.vcf.gz \
    --reference ${mt_ref} \
    --threads ${task.cpus} --no-ansi

    """
}

process NORMALIZE_MUTSERVE_VCF {
    
    publishDir "${params.outdir}/variants/mutserve", mode: 'symlink'
    container params.bcftools
    tag params.sample_id

    input:
    path mutserve_vcf
    path mt_ref

    output:
    tuple path("${mutserve_vcf.getBaseName(2)}.norm.vcf.gz"),path("${mutserve_vcf.getBaseName(2)}.norm.vcf.gz.tbi"), emit: mutserve_norm_vcf

    script:
    """
    set -euo pipefail

    bcftools norm --multiallelics -both ${mutserve_vcf} | bcftools norm --atomize | bcftools view -f PASS -Oz -o ${mutserve_vcf.getBaseName(2)}.norm.vcf.gz
    bcftools index --tbi ${mutserve_vcf.getBaseName(2)}.norm.vcf.gz
    """

}

process ANNOTATE_MUTSERVE_VCF {
    
    publishDir "${params.outdir}/variants/mutserve", mode: 'symlink'
    container params.mitoscope
    tag params.sample_id

    input:
    tuple path(mutserve_norm_vcf), path(mutserve_norm_vcf_index)
    path mitomap_anno_file

    output:
    path("${mutserve_norm_vcf.getBaseName(2)}.annotated.txt")
    path("${mutserve_norm_vcf.getBaseName(2)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail

    annotate.py \
    --input ${mutserve_norm_vcf} \
    --annotations ${mitomap_anno_file} \
    --caller mutserve 
    """


}

process VARIANT_CALLS_SNIFFLES {

    publishDir "${params.outdir}/variants/sniffles", mode: 'symlink'
    container params.sniffles
    tag params.sample_id

    input:
    tuple path(input_bam), path(input_bam_index)

    output:
    path("${input_bam.getBaseName()}.sniffles.vcf"), emit: sniffles_vcf

    script:
    """
    set -euo pipefail

    sniffles --qc-output-all --allow-overwrite \
    --minsupport 2 \
    --input ${input_bam} \
    --vcf ${input_bam.getBaseName()}.sniffles.vcf
    """
}


process FILTER_SNIFFLES_VCF_MINSUPPORT {

    publishDir "${params.outdir}/variants/sniffles", mode: 'symlink'
    container params.bcftools
    tag params.sample_id

    input:
    path sniffles_vcf

    output:
    path("${sniffles_vcf.getBaseName()}.ge2.vcf")

    script:
    """
    set -euo pipefail

    bcftools filter -i "SUPPORT>=2" ${sniffles_vcf} > ${sniffles_vcf.getBaseName()}.ge2.vcf

    """
}