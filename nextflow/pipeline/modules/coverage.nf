process NUCLEAR_COVERAGE {

    publishDir "${params.outdir}", mode: 'symlink'
    container params.mosdepth

    tag "${params.sample_id}"

    input:
    path bam_file
    path regions_bed
    val sample_id

    output:
    path "${sample_id}.mosdepth_summary.txt"
    path "${sample_id}.mosdepth.global.dist.txt"
    path "${sample_id}.regions.bed.gz"

    script:
    """
    set -euo pipefail

    mosdepth \
    --by ${regions_bed} \
    --no-per-base \
    --threads ${task.cpus} \
    ${sample_id} ${bam_file}
    """
}


