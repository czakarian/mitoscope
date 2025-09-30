


process MT_COVERAGE {

    publishDir "${params.outdir}/qc/coverage", mode: 'symlink'
    container params.mosdepth
    tag params.sample_id

    input:
    tuple path(input_bam), path(input_bam_index)

    output:
    path("${input_bam.getBaseName()}.per-base.bed.gz"), emit: per_base_bed
    path("${input_bam.getBaseName()}.mosdepth.summary.txt"), emit: mosdepth_summary

    script:
    """
    set -euo pipefail

    mosdepth ${input_bam.getBaseName()} ${input_bam}
    """

}

process MT_READ_LENGTH {

    publishDir "${params.outdir}/qc/read_length", mode: 'symlink'
    container params.samtools
    tag params.sample_id

    input:
    tuple path(input_bam), path(input_bam_index)

    output:
    path("${input_bam.getBaseName()}.read_lengths.txt"), emit: read_length_file

    script:
    """
    set -euo pipefail

    samtools view ${input_bam} | awk '{print length(\$10)}' > "${input_bam.getBaseName()}.read_lengths.txt"
    """

}

process COVERAGE_PLOT {

    publishDir "${params.outdir}/qc/coverage", mode: 'symlink'
    container params.python
    tag params.sample_id

    input:
    path per_base_bed

    output:
    path("${per_base_bed.getBaseName(3)}.mitochondrial_coverage.png")

    script:
    """
    set -euo pipefail

    qc_plots.py --plot coverage --input ${per_base_bed} --outprefix ${per_base_bed.getBaseName(3)}

    """
}


process READ_LENGTH_PLOT {

    publishDir "${params.outdir}/qc/read_length", mode: 'symlink'
    container params.python
    tag params.sample_id

    input:
    path read_length_file

    output:
    path("${read_length_file.getBaseName(2)}.read_length_distribution.png")

    script:
    """
    set -euo pipefail

    qc_plots.py --plot read_length --input ${read_length_file} --outprefix ${read_length_file.getBaseName(2)}

    """
}

// process QC_SUMMARY {
    
    
//     script:
//     """
//     set -euo pipefail
//     ## get overall qc stats and to qc_summary.txt
//     mean_cov=$(awk '$1 == "MT" {print $4}' ${QCDIR}/${FASTQPREFIX}.mosdepth.summary.txt)
//     read_count=$(cat ${QCDIR}/${FASTQPREFIX}.read_lengths.txt | wc -l)
//     avg_length=$(awk '{sum+=$1} END {print sum/NR}' ${QCDIR}/${FASTQPREFIX}.read_lengths.txt)

//     ## calculate N50
//     n50=$(awk '{sum+=$1; arr[NR]=$1} END {
//         half = sum/2;
//         n = asort(arr);
//         total = 0;
//         for (i = n; i >= 1; i--) {
//             total += arr[i];
//             if (total >= half) {
//                 print arr[i];
//                 break;
//             }
//         }
//     }' ${QCDIR}/${FASTQPREFIX}.read_lengths.txt)

//     # get methylation stats
//     avg_meth=$(awk '{sum+=$7} END {print sum/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)
//     meth_gt_1=$(awk 'NR>1 && $7 > 0.01 {count++} END {print count/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)
//     meth_gt_5=$(awk 'NR>1 && $7 > 0.05 {count++} END {print count/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)
//     meth_gt_10=$(awk 'NR>1 && $7 > 0.1 {count++} END {print count/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)

//     echo -e "Sample\tRead_Count\tMean_Coverage\tAverage_Read_Length\tN50\tAverage_Methylation_Percent\tNum_Meth_Sites_GT_1\tNum_Meth_Sites_GT_5\tNum_Meth_Sites_GT_10" > ${QCDIR}/${FASTQPREFIX}.qc_summary.txt
//     echo -e "${FASTQPREFIX}\t${read_count}\t${mean_cov}\t${avg_length}\t${n50}\t${avg_meth}\t${meth_gt_1}\t${meth_gt_5}\t${meth_gt_10}" >> ${QCDIR}/${FASTQPREFIX}.qc_summary.txt

//     echo '==' $(date) '==' Get basic QC stats -- methylation, mean coverage, read count, avg read length, n50 COMPLETED
//     ##  
//     """   
// }

