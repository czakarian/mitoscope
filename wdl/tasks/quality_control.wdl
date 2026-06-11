version 1.0

task mt_coverage {
    input {
        String sample_id
        File input_bam
        File input_bam_index
    }

    Int threads = 2
    Int mem_gb = 8

    command <<<
        set -euo pipefail

        mosdepth ~{sample_id} ~{input_bam} --threads ~{threads}
    >>>

    output {
        File per_base_bed = "~{sample_id}.per-base.bed.gz"
        File mosdepth_summary = "~{sample_id}.mosdepth.summary.txt"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        docker: "quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0" 
    }

}

task nuclear_coverage {
    input {
        String sample_id
        File input_bam
        File input_bam_index
        File ref
        File nuc_intervals_bed
    }

    Int threads = 4
    Int mem_gb = 8
    Int disk_size = ceil(size(input_bam, "GB") * 4)

    command <<<
        set -euo pipefail

        mosdepth ~{sample_id} ~{input_bam} --fasta ~{ref} --by ~{nuc_intervals_bed} --threads ~{threads} --no-per-base --fast-mode
    >>>

    output {
        File mosdepth_summary = "~{sample_id}.mosdepth.summary.txt"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/biocontainers/mosdepth:0.3.8--hd299d5a_0" 
    }

}

task mt_read_length {
    input {
        String sample_id
        File input_fastq
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        zcat ~{input_fastq} | sed -n '2~4p'  | awk '{print length($0)}' > "~{sample_id}.read_lengths.txt"
    >>>

    output {
        File read_length_txt = "~{sample_id}.read_lengths.txt"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        docker: "czakarian/mitoscope-python:1.0" 
    }

}

task coverage_plot {
    input {
        File py_script
        String sample_id
        File per_base_bed
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        python3 ~{py_script} --plot coverage --input ~{per_base_bed} --outprefix ~{sample_id}
    >>>

    output {
        File qc_plots = "~{sample_id}.mitochondrial_coverage.png"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        docker: "czakarian/mitoscope-python:1.0" 
    }
}


task read_length_plot {
    input {
        File py_script
        String sample_id
        File read_length_file
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        python3 ~{py_script} --plot read_length --input ~{read_length_file} --outprefix ~{sample_id}
    >>>

    output {
        File qc_plots = "~{sample_id}.read_length_distribution.png"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        docker: "czakarian/mitoscope-python:1.0" 
    }
}

task haplogrep {
    input {
        String sample_id
        File mutserve_vcf
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        haplogrep3 classify --extend-report --write-qc \
        --tree phylotree-rcrs@17.2 \
        --input ~{mutserve_vcf} \
        --output ~{sample_id}.haplogrep.txt
    >>>

    output {
        File haplogrep_txt = "~{sample_id}.haplogrep.txt"
        File haplogrep_qc = "~{sample_id}.haplogrep.qc.txt"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        docker: "quay.io/biocontainers/haplogrep3:3.2.2--hdfd78af_1" 
    }
}

task haplocheck {
    input {
        String sample_id
        File mutserve_vcf
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        haplocheck --raw --out ~{sample_id}.haplocheck.txt ~{mutserve_vcf} 
    >>>

    output {
        File haplocheck_txt = "~{sample_id}.haplocheck.txt"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        docker: "quay.io/biocontainers/haplocheck:1.3.3--h2a3209d_2" 
    }
}


task qc_summary {
    input {
        File py_script
        String sample_id
        File mito_coverage_file
        File nuclear_coverage_file 
        File read_length_file
        File haplogrep_file
        File kmer_read_count 
        File chrM_read_count
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        chrM_rc=$(cat ~{chrM_read_count})
        kmer_rc=$(cat ~{kmer_read_count})

        python3 ~{py_script} \
        -c ~{mito_coverage_file} \
        -n ~{nuclear_coverage_file} \
        -r ~{read_length_file} \
        -g ~{haplogrep_file} \
        -m "${chrM_rc}" \
        -k "${kmer_rc}" \
        -s ~{sample_id}
    >>>

    output {
        File tsv = "~{sample_id}.qc_summary.tsv"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        docker: "czakarian/mitoscope-python:1.0" 
    }
}

