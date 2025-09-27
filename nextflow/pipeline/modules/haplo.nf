process HAPLOGREP {
    publishDir "${params.outdir}/qc/haplogroup", mode: 'symlink'
    container params.haplogrep
    tag params.sample_id

    input:
    path mutserve_vcf

    output:
    path("${mutserve_vcf.getBaseName(2)}.haplogrep.*")

    script:
    """
    set -euo pipefail

    haplogrep3 classify --extend-report --write-qc \
    --tree phylotree-rcrs@17.2 \
    --input ${mutserve_vcf} \
    --output ${mutserve_vcf.getBaseName(2)}.haplogrep.txt
    """

}

process HAPLOCHECK {
    publishDir "${params.outdir}/qc/haplogroup", mode: 'symlink'
    container params.haplocheck
    tag params.sample_id

    input:
    path mutserve_vcf

    output:
    path("${mutserve_vcf.getBaseName(2)}.haplocheck.*")

    script:
    """
    set -euo pipefail

    haplocheck --raw --out ${mutserve_vcf.getBaseName(2)}.haplocheck.txt ${mutserve_vcf} 
    """
}



