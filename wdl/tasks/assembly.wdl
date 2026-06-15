version 1.0

task mt_assembly {

  input {
    String sample_id
    File filtered_fastq
    String platform
    Int flye_iterations=10
    Int num_downsampled_reads=100
    Int flye_min_overlap=2500
    Int flye_bp_tolerance=100
  }

  Int threads = 2
  Int mem_gb = 8
  
  Map[String, String] presets = {
    "ont": "--nano-hq",
    "pb": "--pacbio-hifi"
  }
  String flye_preset = presets[platform]

  command <<<
    set -euo pipefail

    iterate="true"
    iter_count=0
    max_iter=~{flye_iterations}
    mito_length_found="false"

    # Log all stdout/stderr to log file and console
    exec > >(tee -a "assembly_iterations.log") 2>&1

    echo "[$(date)] Starting Flye iterative assembly for ~{sample_id}"
    echo "[$(date)] Platform preset: ~{flye_preset}"
    echo "[$(date)] Using up to ${max_iter} iterations"

    while [[ ${iterate} == "true" && ${iter_count} -lt ${max_iter} ]]; do
        iter_count=$((${iter_count} + 1))
        echo "[$(date)] Iteration ${iter_count} for ~{sample_id}"

        iter_dir="assembly_iter${iter_count}"

        # Downsample reads
        seqtk sample -s ${iter_count} ~{filtered_fastq} ~{num_downsampled_reads} | gzip > downsampled.iter${iter_count}.fastq.gz

        # Run Flye assembly
        flye \
        --threads ~{threads} \
        --out-dir "${iter_dir}" \
        -m ~{flye_min_overlap} \
        --meta \
        ~{flye_preset} \
        downsampled.iter${iter_count}.fastq.gz
        
        ## Confirm circularity and length of assembly
        mito_length_found=$(
            awk -v target=16569 -v tol=~{flye_bp_tolerance} '
            NR > 1 {
                len = $2
                circ = $4
                if (circ == "Y" && len >= target - tol && len <= target + tol) {
                    found=1
                    exit
                }
            }
            END {
                print found ? "true" : "false"
            }' "${iter_dir}/assembly_info.txt" 2>/dev/null || echo "false"
        )

        # If mitochondrial-length circular contig found, stop iterating
        if [[ ${mito_length_found} == "true" ]]; then
            echo "[$(date)] Circular ~16.6 kb contig found in iteration ${iter_count}."
            cp -r "${iter_dir}" assembly
            cp downsampled.iter${iter_count}.fastq.gz assembly/downsampled.fastq.gz
            iterate="false"
        else
            echo "[$(date)] No circular mito contig found — continuing."
        fi
    done

    if [[ ${mito_length_found} != "true" ]]; then
        echo "[$(date)] WARNING: Reached max iterations (${max_iter}) without finding a circular mitochondrial contig." >&2
        # Copy the last attempt anyway for inspection
        cp -r "assembly_iter${iter_count}" assembly
    fi

  >>>
  
  output {
    Array[File] mt_assembly_dir = glob("assembly/*")
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "czakarian/mitoscope-flye_seqtk:1.0"
  }

}