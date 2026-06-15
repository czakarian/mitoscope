version 1.0

task make_rotated_mt_ref {

    input {
        File py_script
        File mt_ref
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        python3 ~{py_script} -i ~{mt_ref} -o MT_rotated.fasta
    >>>

    output {
        File mt_ref_rotated = "MT_rotated.fasta"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "czakarian/mitoscope-python:1.0"
    }

}

task make_blast_dbs {
    input {
        File hg38_ref
        File mt_ref
        File mt_ref_rotated
    }

    Int threads = 2
    Int mem_gb = 4

    command <<<
        set -euo pipefail

        makeblastdb -in ~{hg38_ref} -dbtype nucl -out blast/hg38
        makeblastdb -in ~{mt_ref} -dbtype nucl -out blast/mt
        makeblastdb -in ~{mt_ref_rotated} -dbtype nucl -out blast/mt_rotated

        tar -czvf blast.tar.gz blast/

    >>>

    output {
        File blast_db = "blast.tar.gz"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "quay.io/biocontainers/blast:2.16.0--h66d330f_5"
    }

}

task numt_detection_sniffles {
    input {
        String sample_id
        File input_bam
        File input_bam_index
        File ref
    }

    Int threads = 4
    Int mem_gb = 8
    Int disk_size = ceil(size(input_bam, "GB") * 4)

    command <<<
        set -euo pipefail

        sniffles \
        --reference ~{ref} \
        --minsvlen 20 \
        --threads ~{threads} \
        --input ~{input_bam} \
        --vcf ~{sample_id}.numts.INS.sniffles.vcf
    >>>

    output {
        File vcf = "~{sample_id}.numts.INS.sniffles.vcf"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 1
        maxRetries: 1
        docker: "quay.io/biocontainers/sniffles:2.6.2--pyhdfd78af_0"
    }
}

task numt_detection_insertions_to_fasta {

    input {
        String sample_id
        File input_vcf
    }

    Int threads = 2
    Int mem_gb = 4

    command <<<
        set -euo pipefail

        bcftools view -i 'SVTYPE=="INS"' ~{input_vcf} | grep -v '^#' | awk '$5 != "<INS>" {split($10, fields, ":"); print ">" $1 "-" $2 "-" fields[3] "-" fields[4] "\n" $5}' > ~{sample_id}.numts.INS.fasta
    >>>

    output {
        File fasta = "~{sample_id}.numts.INS.fasta"

    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
    }
}


task numt_detection_insertions_blast {
    input {
        String sample_id
        File input_fasta
        File blast_db_tar
    }

    Int threads = 2
    Int mem_gb = 4

    command <<<
        set -euo pipefail

        mkdir "blast_dir"
        tar -xvf ~{blast_db_tar} -C "blast_dir" 

        blastn -query ~{input_fasta} -out "~{sample_id}.numts.INS.blast.txt" -db "blast_dir/blast/mt" -word_size 25 -outfmt 6 -num_threads ~{threads} -task blastn -perc_identity 95
        blastn -query ~{input_fasta} -out "~{sample_id}.numts.INS.blast.rotated.txt" -db "blast_dir/blast/mt_rotated" -word_size 25 -outfmt 6 -num_threads ~{threads} -task blastn -perc_identity 95

    >>>

    output {
        File blast_output =  "~{sample_id}.numts.INS.blast.txt"
        File blast_rotated_output =  "~{sample_id}.numts.INS.blast.rotated.txt"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "quay.io/biocontainers/blast:2.16.0--h66d330f_5"
    }

}


task numt_detection_mtdna_insertions_to_fasta {

    input {
        String sample_id
        File input_fasta
        File blast_output
    }

    Int threads = 2
    Int mem_gb = 4

    command <<<
        set -euo pipefail

        cut -f 1 ~{blast_output} | sort | uniq > blast_record_names.txt
        seqtk subseq ~{input_fasta} blast_record_names.txt > ~{sample_id}.numts.INS.blast.mtDNA.fasta
    >>>

    output {
        File fasta = "~{sample_id}.numts.INS.blast.mtDNA.fasta"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "czakarian/mitoscope-flye_seqtk:1.0"
    }
}

task numt_detection_mtdna_insertions_blast_check {
    input {
        String sample_id
        File input_fasta
        File blast_db_tar
    }

    Int threads = 2
    Int mem_gb = 4

    command <<<
        set -euo pipefail

        mkdir "blast_dir"
        tar -xvf ~{blast_db_tar} -C "blast_dir" 

        blastn -query ~{input_fasta} -db "blast_dir/blast/hg38" -word_size 25 -outfmt 6 -num_threads ~{threads} -task blastn -perc_identity 95 -max_target_seqs 1 | awk '$2 == "chrM"' > "~{sample_id}.numts.INS.blast.mtDNA.txt"
        
        echo -e "nuc_chrom\tnuc_position\tmt_start\tmt_end\tDR\tDV\tAF" > "~{sample_id}.numts.INS.tsv"
        awk '{split($1, id_parts, "-"); print id_parts[1]"\t"id_parts[2]"\t"$9"\t"$10"\t"id_parts[3]"\t"id_parts[4]"\t"id_parts[4]/(id_parts[3]+id_parts[4])}' "~{sample_id}.numts.INS.blast.mtDNA.txt" >> "~{sample_id}.numts.INS.tsv"
        
    >>>

    output {
        File blast_output = "~{sample_id}.numts.INS.blast.mtDNA.txt"
        File tsv_output = "~{sample_id}.numts.INS.tsv"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "quay.io/biocontainers/blast:2.16.0--h66d330f_5"
    }
}

task numt_detection_plot {
    input {
        File py_script
        String sample_id
        File blast_file
        File circos_bed
    }

    Int threads = 1
    Int mem_gb = 2

    command <<<
        set -euo pipefail
        
        # set up temp cache directory for matplotlib
        export MPLCONFIGDIR="/tmp/mplconfig"
        mkdir -p $MPLCONFIGDIR

        python3 ~{py_script} -i ~{blast_file} -o "~{sample_id}.numts.INS.circos.png" -b ~{circos_bed}
    >>>

    output {
        File circos_plot = "~{sample_id}.numts.INS.circos.png"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "czakarian/mitoscope-python:1.0"
    }
}

task numt_detection_supplementary {
    input {
        File py_script
        String sample_id
        File input_bam
        File input_bam_index
        File hg38_ref
    }

    Int threads = 1
    Int mem_gb = 4
    Int disk_size = ceil(size(input_bam, "GB") * 4)

    command <<<
        set -euo pipefail

        python3 ~{py_script} -i ~{input_bam} -o ~{sample_id} -r ~{hg38_ref}
    >>>

    output {
        File supplementary_numt_tsv = "~{sample_id}.numts.SA.tsv"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        disks: "local-disk " + disk_size + " SSD"
        preemptible: 1
        maxRetries: 1
        docker: "czakarian/mitoscope-python:1.0"
    }

}

