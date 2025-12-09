process VARIANT_CALLS_BALDUR {

    publishDir "${params.outdir}/${sample_id}/variants/baldur/", pattern: "*.baldur_del.txt", mode: 'copy'
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

    // publishDir "${params.outdir}/${sample_id}/variants/baldur/", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_vcf)

    output:
    tuple val(sample_id), path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz"),path("${baldur_vcf.getBaseName(2)}.norm.vcf.gz.tbi"), emit: norm_vcf

    script:
    """
    set -euo pipefail

    bcftools norm --multiallelics -both ${baldur_vcf} | bcftools norm --atomize --atom-overlaps . | bcftools view -f PASS -Oz -o ${baldur_vcf.getBaseName(2)}.norm.vcf.gz 
    bcftools index --tbi ${baldur_vcf.getBaseName(2)}.norm.vcf.gz  

    """
}


process ANNOTATE_BALDUR {
    
    //publishDir "${params.outdir}/${sample_id}/variants/baldur/", pattern: "*.mitomap.txt", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/variants/baldur/", pattern: "*.heteroplasmy.png", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_norm_vcf), path(baldur_norm_vcf_index)
    path mitomap_anno_file

    output:
    tuple val(sample_id), path("${baldur_norm_vcf.getBaseName(2)}.mitomap.txt"), emit: mitomap_txt
    path("${baldur_norm_vcf.getBaseName(2)}.heteroplasmy.png") 

    script:
    """
    set -euo pipefail

    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    annotate.py \
    --input ${baldur_norm_vcf} \
    --annotations ${mitomap_anno_file} \
    --caller baldur 
    """
}

process VEP_BALDUR_VCF {
    //publishDir "${params.outdir}/${sample_id}/variants/baldur/", mode: 'copy'
    container params.vep
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(baldur_norm_vcf), path(baldur_norm_vcf_index)
    tuple path(mt_ref), path(mt_ref_index)
    path vep_cache_dir

    output:
    tuple val(sample_id), path("${baldur_norm_vcf.getBaseName(2)}.vep.vcf.gz"), path("${baldur_norm_vcf.getBaseName(2)}.vep.vcf.gz.tbi"), emit: vep_vcf
    // path("${baldur_norm_vcf.getBaseName(2)}.vep.tsv")

    script:
    """
    set -euo pipefail

    vep --input_file ${baldur_norm_vcf} \
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
    --distance 0

    tabix -p vcf ${baldur_norm_vcf.getBaseName(2)}.vep.vcf.gz

    vep --input_file ${baldur_norm_vcf} \
    --cache \
    --dir_cache ${vep_cache_dir} \
    --offline \
    --fork ${task.cpus} \
    --fasta ${mt_ref} \
    --tab \
    --output_file ${baldur_norm_vcf.getBaseName(2)}.vep.tsv \
    --hgvs \
    --protein \
    --symbol \
    --biotype \
    --no_stats \
    --allow_non_variant \
    --distance 0

    """

}

process ADD_MITOMAP_TO_BALDUR_VCF {
    publishDir "${params.outdir}/${sample_id}/variants/baldur/", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_vcf), path(input_vcf_index), path(mitomap_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.mt.baldur.annotated.vcf.gz"), path("${sample_id}.mt.baldur.annotated.vcf.gz.tbi"), emit: anno_vcf

    script:
    """
    
    cut -f 1-5,35 ${mitomap_tsv} | tail +2 > subset_mitomap.tab
    bgzip subset_mitomap.tab
    tabix -s1 -b2 -e2 subset_mitomap.tab.gz

    echo '##INFO=<ID=MITOMAP,Number=1,Type=String,Description="Annotations from MITOMAP. Format: Gene.Name|Gene.Type|Amino.Acid.Change|GB.Freq.FL|GB.Freq.CR|GB.Seqs.FL|GB.Seqs.CR|Homoplasmy|Heteroplasmy|Disease.Status|References|Additional.Annotations|MitoTIP">' > header.txt

    bcftools annotate ${input_vcf} \
        -c CHROM,POS,-,REF,ALT,INFO/MITOMAP \
        --annotations subset_mitomap.tab.gz \
        --header-lines header.txt \
        --threads ${task.cpus} \
        --output-type z \
        --write-index=tbi \
        --output ${sample_id}.mt.baldur.annotated.vcf.gz 

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
    
    publishDir "${params.outdir}/${sample_id}/variants/sniffles/", mode: 'copy'
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

