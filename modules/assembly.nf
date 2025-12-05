process MT_ASSEMBLY {

    publishDir path: "${params.outdir}/${sample_id}/", pattern:"assembly",  mode: 'copy'
    publishDir path: "${params.outdir}/${sample_id}/logs", pattern:"*.log",  mode: 'copy'
    container params.flye
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(filtered_fastq)
    val platform

    output:
    tuple val(sample_id), path("assembly"), emit: assembly_dir
    path("assembly_iterations.log"), emit: log

    script:
    def flye_preset = platform == 'ont' ? '--nano-hq' : 
                      platform == 'pb'  ? '--pacbio-hifi':
                      null
    """
    set -euo pipefail

    iterate="true"
    iter_count=0
    max_iter=${params.flye_iterations}
    mito_length_found="false"

    # Log all stdout/stderr to log file and console
    exec > >(tee -a "assembly_iterations.log") 2>&1

    echo "[\$(date)] Starting Flye iterative assembly for ${sample_id}"
    echo "[\$(date)] Platform preset: ${flye_preset}"
    echo "[\$(date)] Using up to \${max_iter} iterations"

    while [[ \${iterate} == "true" && \${iter_count} -lt \${max_iter} ]]; do
        iter_count=\$((\${iter_count} + 1))
        echo "[\$(date)] Iteration \${iter_count} for ${sample_id}"

        iter_dir="assembly_iter\${iter_count}"

        # Downsample reads
        seqtk sample -s \${iter_count} ${filtered_fastq} ${params.num_downsampled_reads} | gzip > downsampled.iter\${iter_count}.fastq.gz

        # Run Flye assembly
        flye \
        --threads ${task.cpus} \
        --out-dir "\${iter_dir}" \
        -m ${params.flye_min_overlap} \
        --meta \
        ${flye_preset} \
        downsampled.iter\${iter_count}.fastq.gz
        
        ## Confirm circularity and length of assembly
        mito_length_found=\$(
            awk -v target=16569 -v tol=${params.flye_bp_tolerance} '
            NR > 1 {
                len = \$2
                circ = \$4
                if (circ == "Y" && len >= target - tol && len <= target + tol) {
                    found=1
                    exit
                }
            }
            END {
                print found ? "true" : "false"
            }' "\${iter_dir}/assembly_info.txt" 2>/dev/null || echo "false"
        )

        # If mitochondrial-length circular contig found, stop iterating
        if [[ \${mito_length_found} == "true" ]]; then
            echo "[\$(date)] Circular ~16.6 kb contig found in iteration \${iter_count}."
            cp -r "\${iter_dir}" assembly
            cp downsampled.iter\${iter_count}.fastq.gz assembly/downsampled.fastq.gz
            iterate="false"
        else
            echo "[\$(date)] No circular mito contig found — continuing."
        fi
    done

    if [[ \${mito_length_found} != "true" ]]; then
        echo "[\$(date)] WARNING: Reached max iterations (\${max_iter}) without finding a circular mitochondrial contig." >&2
        # Copy the last attempt anyway for inspection
        cp -r "assembly_iter\${iter_count}" assembly
    fi

    """
}

process INDEX_ASSEMBLY {

    publishDir "${params.outdir}/${sample_id}/assembly", mode: 'copy'
    container params.samtools
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly_fasta)

    output:
    tuple val(sample_id), path(assembly_fasta), path("*.fai")

    script:
    """
    samtools faidx ${assembly_fasta}
    """
}


process ROTATE_ASSEMBLY {

    publishDir "${params.outdir}/${sample_id}/assembly", pattern: "*.fasta", mode: 'symlink'
    publishDir path: "${params.outdir}/logs", pattern: "*.log", mode: 'symlink'
    container params.python
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(assembly_fasta), path(assembly_fasta_fai)
    tuple val(sample_id), path(assembly_align_ref_bam), path(assembly_align_ref_bam_index)
    tuple val(sample_id), path(assembly_dir)

    output:
    tuple val(sample_id), path("assembly_rotated.fasta"), optional: true
    // path("rotate_assembly.log"), emit: log

    script:
    """
    rotate.py --bam ${assembly_align_ref_bam} --assembly ${assembly_fasta}
    """
}



