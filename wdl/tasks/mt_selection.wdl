version 1.0

task kmer_selection {

  input {
    String sample_id
    File fastq_file 
    File kmc_pre
    File kmc_suf
  }

  Int threads = 4
  Int mem_gb = 8

  command <<<
    set -euo pipefail
    
    # Reconstruct KMC DB from staged components
    cp ~{kmc_pre} MT.k29.kmc_pre
    cp ~{kmc_suf} MT.k29.kmc_suf

    kmc_tools -t~{threads} filter MT.k29 -ci1 ~{fastq_file} -fq -ci2500 /dev/stdout | tr ' ' '\t' | gzip > ~{sample_id}.kmer_selection.fastq.gz
    echo $(($(zcat ~{sample_id}.kmer_selection.fastq.gz | wc -l) / 4)) > ~{sample_id}.kmer_read_count.txt

  >>>

  output {
    File kmer_selected_fastq = "~{sample_id}.kmer_selection.fastq.gz"
    File kmer_read_count = "~{sample_id}.kmer_read_count.txt"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    docker: "quay.io/biocontainers/kmc:3.2.1--hf1761c0_2"
  }

}

task align_to_ref {

  input {
    String sample_id
    File fastq_file 
    String platform
    File minimap_index
  }

  Int threads = 4
  Int mem_gb = 8

  Map[String, String] presets = {
    "ont": "map-ont",
    "pb": "map-hifi"
  }
  String preset = presets[platform]

  command <<<
    set -euo pipefail

    minimap2 -ax ~{preset} -Y -y -t ~{threads} ~{minimap_index} ~{fastq_file} > "~{sample_id}.sam"

  >>>

  output {
    File sam = "~{sample_id}.sam"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    docker: "quay.io/biocontainers/minimap2:2.24--h7132678_1"
  }

}

task sam_to_bam {

  input {
    String sample_id
    File sam 
  }

  Int threads = 2
  Int mem_gb = 4

  command <<<
    set -euo pipefail

    samtools sort -@~{threads} -o "~{sample_id}.bam" ~{sam}
    samtools index -@~{threads} "~{sample_id}.bam"

  >>>

  output {
    File bam = "~{sample_id}.bam"
    File bam_index = "~{sample_id}.bam.bai"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    docker: "quay.io/biocontainers/samtools:1.21--h96c455f_1"
  }

}


task filter_numts {

  input {
    File py_script
    String sample_id
    File bam_file
    File bam_file_index
  }

  Int threads = 1
  Int mem_gb = 4

  command <<<
    set -euo pipefail
    
    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR="/tmp/mplconfig"
    mkdir -p $MPLCONFIGDIR

    python3 ~{py_script} -i ~{bam_file} -o ~{sample_id} \
    --max_sc_threshold 200 \
    --max_meth_threshold 0.5

  >>>

  output {
    File mt_bam = "~{sample_id}.mt.bam"
    File mt_bam_index = "~{sample_id}.mt.bam.bai"
    File discard_bam = "~{sample_id}.discard.bam"
    File discard_bam_index = "~{sample_id}.discard.bam.bai"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    docker: "czakarian/mitoscope-python:1.0"
  }

}

task filtered_bam_to_fastq {

  input {
    String sample_id
    File bam_file
    File bam_file_index
  }

  Int threads = 2
  Int mem_gb = 8

  command <<<
    set -euo pipefail

    samtools fastq -T MM,ML -@ ~{threads} -0 ~{sample_id}.fastq.gz ~{bam_file}
  >>>

  output {
    File filtered_fastq = "~{sample_id}.fastq.gz"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    docker: "quay.io/biocontainers/samtools:1.21--h96c455f_1"
  }

}