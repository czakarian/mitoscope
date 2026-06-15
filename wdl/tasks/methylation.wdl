version 1.0

task meth_freq_ont {

  input {
    String sample_id
    File input_bam
    File input_bam_index
    File mt_ref
    File mt_ref_index
  }

  Int threads = 2
  Int mem_gb = 8

  command <<<
    set -euo pipefail

    modkit pileup ~{input_bam} ~{sample_id}.modkit.bedmethyl --ref ~{mt_ref} --threads ~{threads} \
    --motif CG 0 --ignore a --log-filepath modkit.log --header --no-filtering

    modkit pileup ~{input_bam} ~{sample_id}.modkit.CH.bedmethyl --ref ~{mt_ref} --threads ~{threads} \
    --motif CG 0 --motif CH 0 --ignore a --log-filepath modkit.CH.log --header --no-filtering

    modkit pileup ~{input_bam} ~{sample_id}.modkit.A.bedmethyl --ref ~{mt_ref} --threads ~{threads} \
    --motif A 0 --log-filepath modkit.A.log --header --no-filtering
    
  >>>

  output {
    File bedmethyl_CG = "~{sample_id}.modkit.bedmethyl"
    File bedmethyl_CH = "~{sample_id}.modkit.CH.bedmethyl"
    File bedmethyl_A = "~{sample_id}.modkit.A.bedmethyl"

  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "quay.io/biocontainers/ont-modkit:0.5.0--hcdda2d0_2"
  }

}

task meth_freq_pb {

  input {
    String sample_id
    File input_bam
    File input_bam_index
    File mt_ref
    File mt_ref_index
  }

  Int threads = 2
  Int mem_gb = 8

  command <<<
    set -euo pipefail

    aligned_bam_to_cpg_scores --bam ~{input_bam} \
    --output-prefix ~{sample_id}.pbcpgtools \
    --pileup-mode count \
    --threads ~{threads}
    
  >>>

  output {
    File bed = "~{sample_id}.pbcpgtools.combined.bed.gz"
    File bed_index = "~{sample_id}.pbcpgtools.combined.bed.gz.tbi"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "quay.io/biocontainers/pb-cpg-tools:3.0.0--h9ee0642_0"
  }

}
