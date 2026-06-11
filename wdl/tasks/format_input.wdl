version 1.0

task aligned_bam_to_fastq {

  input {
    String sample_id
    File bam_file 
    File bam_file_index
  }

  Int threads = 2
  Int mem_gb = 8
  Int disk_size = ceil(size(bam_file, "GB") * 4)

  command <<<
    set -euo pipefail

    samtools view -h -@ ~{threads} ~{bam_file} chrM \
    | samtools fastq -T MM,ML -@ ~{threads} \
    | tr '\t' ' ' > ~{sample_id}.fastq
    
    samtools view -@ ~{threads} -c -F 2308 ~{bam_file} chrM > ~{sample_id}.chrM_read_count.txt
  >>>

  output {
    File fastq = "~{sample_id}.fastq"
    File chrM_read_count = "~{sample_id}.chrM_read_count.txt"
  }

  runtime {
    cpu: threads
    memory: mem_gb + "G"
    disks: "local-disk " + disk_size + " SSD"
    docker: "quay.io/biocontainers/samtools:1.21--h96c455f_1"
  }

}

task compress_fastq {

  input {
    String sample_id
    File fastq_file 
  }

  Int threads = 2
  Int mem_gb = 8
  Int disk_size = ceil(size(fastq_file, "GB") * 4)

  command <<<
    set -euo pipefail

    pigz -p ~{threads} -c ~{fastq_file} > "~{sample_id}.fastq.gz"
  >>>

  output {
    File fastq_gz = "~{sample_id}.fastq.gz"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    disks: "local-disk " + disk_size + " SSD"
    docker: "quay.io/biocontainers/pigz:2.8"
  }

}