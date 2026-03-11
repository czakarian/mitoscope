process KMER_SELECTION {

    // publishDir "${params.outdir}/${sample_id}", mode: 'copy'
    container params.kmctools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fastq_file)
    path kmc_pre
    path kmc_suf

    output:
    tuple val(sample_id), path("${sample_id}.kmer_selection.fastq.gz"), emit:fastq
    tuple val(sample_id), path("${sample_id}.kmer_read_count.txt"), emit:kmer_read_count

    script:
    """
    set -euo pipefail
    
    # Reconstruct KMC DB from staged components
    cp ${kmc_pre} MT.k29.kmc_pre
    cp ${kmc_suf} MT.k29.kmc_suf

    kmc_tools -t${task.cpus} filter MT.k29 -ci1 ${fastq_file} -fq -ci2500 /dev/stdout | tr ' ' '\t' | gzip > ${sample_id}.kmer_selection.fastq.gz
    echo \$((\$(zcat ${sample_id}.kmer_selection.fastq.gz | wc -l) / 4)) > ${sample_id}.kmer_read_count.txt

    """
}