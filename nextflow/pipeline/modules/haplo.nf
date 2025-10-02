process HAPLOGREP {
    publishDir "${params.outdir}/${sample_id}/qc/haplogroup", mode: 'symlink'
    container params.haplogrep
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mutserve_vcf)

    output:
    path("${mutserve_vcf.getBaseName(2)}.haplogrep.*")

    script:
    """
    set -euo pipefail

    haplogrep3 classify --extend-report --write-qc \
    --tree ${params.haplogrep_tree} \
    --input ${mutserve_vcf} \
    --output ${mutserve_vcf.getBaseName(2)}.haplogrep.txt
    """

}

process HAPLOCHECK {
    publishDir "${params.outdir}/${sample_id}/qc/haplogroup", mode: 'symlink'
    container params.haplocheck
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mutserve_vcf)

    output:
    path("${mutserve_vcf.getBaseName(2)}.haplocheck.*")

    script:
    """
    set -euo pipefail

    haplocheck --raw --out ${mutserve_vcf.getBaseName(2)}.haplocheck.txt ${mutserve_vcf} 
    """
}



