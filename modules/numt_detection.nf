process NUMT_DETECTION_SNIFFLES {

    // publishDir "${params.outdir}/${sample_id}/numts/", pattern: "*.{vcf,snf}", mode: 'copy'
    publishDir "${params.outdir}/${sample_id}/logs/", pattern: "*.log", mode: 'copy'
    container params.sniffles
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_cram), path(input_cram_index)
    path ref

    output:
    tuple val(sample_id), path("${sample_id}.numts.INS.sniffles.vcf"), emit: vcf
    path("numt_detection.sniffles.log"), emit:log

    script:
    """
    set -euo pipefail

    # Log all stdout/stderr to log file and console
    exec > >(tee -a "numt_detection.sniffles.log") 2>&1


    sniffles \
    --reference ${ref} \
    --minsvlen 20 \
    --threads ${task.cpus} \
    --input ${input_cram} \
    --vcf ${sample_id}.numts.INS.sniffles.vcf
    """
}

process NUMT_DETECTION_INSERTIONS_TO_FASTA {

    //publishDir "${params.outdir}/${sample_id}/numts/", pattern: "*.{fasta}", mode: 'copy'
    container params.bcftools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.numts.INS.fasta"), emit: fasta

    script:
    """
    set -euo pipefail

    bcftools view -i 'SVTYPE=="INS"' ${input_vcf} | grep -v '^#' | awk '\$5 != "<INS>" {print ">" \$1 "-" \$2 "\\n" \$5}' > ${sample_id}.numts.INS.fasta

    """
}

process NUMT_DETECTION_INSERTIONS_BLAST {

    //publishDir "${params.outdir}/${sample_id}/numts/", pattern: "*.{txt}", mode: 'copy'
    container params.blast
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_fasta)
    path blast_db

    output:
    tuple val(sample_id), path("${sample_id}.numts.INS.blast.txt")

    script:
    """
    set -euo pipefail

    blastn -query ${input_fasta} -out "${sample_id}.numts.INS.blast.txt" -db ${blast_db}/mt.db -word_size 25 -outfmt 6 -num_threads ${task.cpus} -task blastn -perc_identity 95

    """
}

process NUMT_DETECTION_MTDNA_INSERTIONS_TO_FASTA {

    publishDir "${params.outdir}/${sample_id}/numts/", pattern: "*.{fasta}", mode: 'copy'
    container params.flye
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_fasta), path(blast_output)

    output:
    tuple val(sample_id), path("${sample_id}.numts.INS.blast.mtDNA.fasta")

    script:
    """
    set -euo pipefail

    cut -f 1 ${blast_output} | sort | uniq > blast_record_names.txt
    seqtk subseq ${input_fasta} blast_record_names.txt > ${sample_id}.numts.INS.blast.mtDNA.fasta

    """
}

process NUMT_DETECTION_MTDNA_INSERTIONS_BLAST_CHECK {

    publishDir "${params.outdir}/${sample_id}/numts/", pattern: "*.{txt}", mode: 'copy'
    container params.blast
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_fasta)
    path blast_db

    output:
    tuple val(sample_id), path("${sample_id}.numts.INS.blast.mtDNA.txt")

    script:
    """
    set -euo pipefail

    blastn -query ${input_fasta} -db ${blast_db}/hg38 -word_size 25 -outfmt 6 -num_threads ${task.cpus} -task blastn -perc_identity 95 -max_target_seqs 1 | awk '\$2 == "chrM"' > "${sample_id}.numts.INS.blast.mtDNA.txt"
    """
}


process NUMT_DETECTION_SUPPLEMENTARY {

    publishDir "${params.outdir}/${sample_id}/numts/", pattern: "*.{bam,bai,txt,csv}", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_cram), path(input_cram_index)
    path ref

    output:
    // tuple val(sample_id), path("${sample_id}.numts.SA.bam"), emit:bam
    // path('examine.txt')
    path("${sample_id}.numts.SA.csv")

    script:
    """
    set -euo pipefail

    numt_detection_sa.py -i ${input_cram} -o ${sample_id} -r ${ref} > examine.txt

    """
}


process NUMT_DETECTION_PLOT {

    publishDir "${params.outdir}/${sample_id}/numts/", pattern: "*.{png}", mode: 'copy'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(blast_file)
    path circos_bed

    output:
    tuple val(sample_id), path("${sample_id}.numts.INS.circos.png")

    script:
    """
    set -euo pipefail
    
    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=${params.mplconfigdir}
    mkdir -p \$MPLCONFIGDIR

    numt_circos.py -i ${blast_file} -o "${sample_id}.numts.INS.circos.png" -b ${circos_bed}

    """
}