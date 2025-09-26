process FILTER_NUMTS {

    publishDir "${params.outdir}/alignments", mode: 'symlink'
    // container params.samtools

    tag "${params.sample_id}"

    input:
    tuple path(bam_file), path(bam_file_index)

    output:
    path("${bam_file.baseName}.filtered.bam"), emit: mt_filtered_bam
    path("${bam_file.baseName}.discardReads.bam"), emit: mt_numt_bam


    script:
    """
    set -euo pipefail
    
    filter_bam.py --max_sc_threshold 100 -i ${bam_file}
    """
}