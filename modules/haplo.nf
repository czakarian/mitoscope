process HAPLOGREP {
    publishDir "${params.outdir}/${sample_id}/qc/haplogroup", mode: 'copy'
    container params.haplogrep
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(mutserve_vcf)

    output:
    tuple val(sample_id), path("${mutserve_vcf.getBaseName(2)}.haplogrep.txt"), emit: haplogrep_txt
    tuple val(sample_id), path("${mutserve_vcf.getBaseName(2)}.haplogrep.qc.txt"), emit: haplogrep_qc


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
    publishDir path:"${params.outdir}/${sample_id}/qc/haplogroup", pattern: "*.txt", mode: 'copy'

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



