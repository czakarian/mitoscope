process VARIANT_CALLS_BALDUR {

    publishDir "${params.outdir}/${sample_id}/variants/baldur/", pattern: "*filtered.baldur*", mode: 'copy'
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
    --small-deletion-limit 20 \
    --large-deletion-limit 20 \
    -q 20 -Q 20 \
    -T ${mt_ref} \
    -n ${sample_id} \
    -o "${input_bam.getBaseName()}.baldur" \
    ${input_bam} 
    """
}

process NORMALIZE_BALDUR_VCF {

    publishDir "${params.outdir}/${sample_id}/variants/baldur/", pattern: "*filtered.baldur*", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_vcf)

    output:
    tuple val(sample_id), path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz.tbi"), emit: norm_vcf
    tuple val(sample_id), path("${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz.tbi"), emit: norm_indels_vcf
    tuple val(sample_id), path("${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz.tbi"), emit: norm_snvs_vcf

    script:
    """
    set -euo pipefail

    bcftools norm --multiallelics -both ${baldur_vcf} | bcftools norm --atomize --atom-overlaps . | bcftools view -f PASS -Oz -o ${baldur_vcf.getBaseName(2)}.norm.vcf.gz 
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.vcf.gz  

    bcftools view --types indels ${baldur_vcf.getBaseName(2)}.norm.vcf.gz -Oz -o ${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz 
    bcftools view --exclude-types indels ${baldur_vcf.getBaseName(2)}.norm.vcf.gz -Oz -o ${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz 
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.indels.vcf.gz 
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.snvs.vcf.gz 

    """
}

process ANNOTATE_BALDUR_SNVS {
    
    publishDir "${params.outdir}/${sample_id}/variants/baldur/", pattern: "*filtered.baldur*", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_norm_snvs_vcf), path(baldur_norm_snvs_vcf_index)
    path mitomap_anno_file

    output:
    path("${baldur_norm_snvs_vcf.getBaseName(2)}.annotated.txt")
    path("${baldur_norm_snvs_vcf.getBaseName(2)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail

    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    annotate.py \
    --input ${baldur_norm_snvs_vcf} \
    --annotations ${mitomap_anno_file} \
    --caller baldur 
    """
}

process ANNOTATE_BALDUR_INDELS {
    
    publishDir "${params.outdir}/${sample_id}/variants/baldur/", pattern: "*filtered.baldur*", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_norm_indels_vcf), path(baldur_norm_indels_vcf_index)
    path mitomap_anno_file

    output:
    path("${baldur_norm_indels_vcf.getBaseName(2)}.annotated.txt")
    path("${baldur_norm_indels_vcf.getBaseName(2)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail
    
    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    annotate.py \
    --input ${baldur_norm_indels_vcf} \
    --annotations ${mitomap_anno_file} \
    --caller baldur 
    """
}

process VARIANT_CALLS_MUTSERVE {

    publishDir "${params.outdir}/${sample_id}/variants/mutserve/", pattern: "*filtered.mutserve*", mode: 'copy'
    container params.mutserve
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)
    tuple path(mt_ref), path(mt_ref_index)
    val contig_name

    output:
    tuple val(sample_id), path("${input_bam.getBaseName()}.mutserve.vcf.gz"), emit: vcf

    script:
    """
    set -euo pipefail

    mutserve call ${input_bam} \
    --output ${input_bam.getBaseName()}.mutserve.vcf.gz \
    --reference ${mt_ref} \
    --contig-name ${contig_name} \
    --threads ${task.cpus} --no-ansi \
    --alignQ 0 --mapQ 20 --baseQ 20 --level 0.005

    """
}

process NORMALIZE_MUTSERVE_VCF {
    
    publishDir "${params.outdir}/${sample_id}/variants/mutserve/", pattern: "*filtered.mutserve*", mode: 'copy'
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
    
    publishDir "${params.outdir}/${sample_id}/variants/mutserve/", pattern: "*filtered.mutserve*", mode: 'copy'
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

    publishDir "${params.outdir}/${sample_id}/variants/sniffles/", pattern: "*filtered.sniffles*", mode: 'copy'
    //publishDir "${params.outdir}/${sample_id}/variants/sniffles/assembly_to_ref", pattern: "*assembly.ref.sniffles*", mode: 'copy'
    container params.sniffles
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam), path(input_bam_index)

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
    --input ${input_bam} \
    --snf ${input_bam.getBaseName()}.sniffles.snf \
    --vcf ${input_bam.getBaseName()}.sniffles.vcf

    """
}


process FILTER_SNIFFLES_VCF {
    
    publishDir "${params.outdir}/${sample_id}/variants/sniffles/", pattern: "*filtered.sniffles*", mode: 'copy'
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