version 1.0

task baldur {

  input {
    String sample_id
    File input_bam
    File input_bam_index
    File mt_ref
    File mt_ref_index
    ## baldur-specific parameters
    Int deletion_size_limit=25
    Float indel_threshold=0.1
    Float snv_hard_threshold=0.0005
    Float snv_soft_threshold=0.0025
    Int homopolymer_limit=5
    Int mapq_threshold=20
    Int base_qual_threshold=20
    Int max_indel_base_qual=30
  }

  Int threads = 4
  Int mem_gb = 8

  command <<<
    set -euo pipefail

    baldur -l debug \
    --output-deletions \
    --small-deletion-limit ~{deletion_size_limit} \
    --large-deletion-limit ~{deletion_size_limit} \
    --indel-thresholds ~{indel_threshold} ~{indel_threshold} \
    --snv-thresholds ~{snv_hard_threshold} ~{snv_soft_threshold}  \
    --homopolymer-limit ~{homopolymer_limit} \
    -q~{mapq_threshold} -Q~{base_qual_threshold} -I~{max_indel_base_qual} \
    -T ~{mt_ref} \
    -n ~{sample_id} \
    -o ~{sample_id}.baldur \
    ~{input_bam} 

  >>>

  output {
    File baldur_vcf = "~{sample_id}.baldur.vcf.gz"
    File? baldur_dels = "~{sample_id}.baldur_del.txt"

  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "czakarian/mitoscope-baldur:1.0"
  }

}

task normalize_vcf {

  input {
    String sample_id
    File baldur_vcf
  }

  Int threads = 1
  Int mem_gb = 2

  command <<<
    set -euo pipefail

    bcftools norm --multiallelics -both ~{baldur_vcf} \
    | bcftools norm --atomize --atom-overlaps . \
    | bcftools view -f PASS -Oz -o ~{sample_id}.baldur.norm.vcf.gz 

    bcftools index --tbi ~{sample_id}.baldur.norm.vcf.gz  

  >>>

  output {
    File baldur_norm_vcf = "~{sample_id}.baldur.norm.vcf.gz"
    File baldur_norm_vcf_index = "~{sample_id}.baldur.norm.vcf.gz.tbi"

  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
  }

}

task pull_mitomap_annos {

  input {
    File py_script
    String sample_id
    File baldur_norm_vcf
    File baldur_norm_vcf_index
    File mitomap_anno_file
  }

  Int threads = 1
  Int mem_gb = 2

  command <<<
    set -euo pipefail

    # set up temp cache directory for matplotlib
    export MPLCONFIGDIR=/tmp/mplconfig
    mkdir -p $MPLCONFIGDIR

    python3 ~{py_script} \
    --input ~{baldur_norm_vcf} \
    --output ~{sample_id} \
    --annotations ~{mitomap_anno_file} \
    --caller baldur 
  >>>

  output {
    File mitomap_txt= "~{sample_id}.mitomap.txt"
    File heteroplasmy_plot = "~{sample_id}.heteroplasmy.png"

  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "czakarian/mitoscope-python:1.0"
  }

}

task run_vep {

    input {
        String sample_id
        File baldur_norm_vcf
        File baldur_norm_vcf_index
        File mt_ref
        File vep_cache_dir_tar
    }

    Int threads = 2
    Int mem_gb = 2

    command <<<
        set -euo pipefail

        mkdir "vep_dir"
        tar -xvf ~{vep_cache_dir_tar} -C "vep_dir" 

        vep --input_file ~{baldur_norm_vcf} \
        --cache \
        --dir_cache "vep_dir/vep" \
        --offline \
        --fork ~{threads} \
        --fasta ~{mt_ref} \
        --vcf \
        --compress_output bgzip \
        --output_file ~{sample_id}.vep.vcf.gz \
        --hgvs \
        --protein \
        --symbol \
        --biotype \
        --no_stats \
        --allow_non_variant \
        --distance 0

        tabix -p vcf ~{sample_id}.vep.vcf.gz

    >>>

    output {
        File vep_vcf = "~{sample_id}.vep.vcf.gz"
        File vep_vcf_index = "~{sample_id}.vep.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "quay.io/biocontainers/ensembl-vep:115--pl5321h2a3209d_0"
    }
}

task add_mitomap_to_vcf {
    input {
        String sample_id
        File input_vcf
        File input_vcf_index
        File mitomap_txt
    }

    Int threads = 2
    Int mem_gb = 2

    command <<<
    
        cut -f 1-5,35 ~{mitomap_txt} | tail +2 > subset_mitomap.tab
        bgzip subset_mitomap.tab
        tabix -s1 -b2 -e2 subset_mitomap.tab.gz

        echo '##INFO=<ID=MITOMAP,Number=1,Type=String,Description="Annotations from MITOMAP. Format: Gene.Name|Gene.Type|Amino.Acid.Change|GB.Freq.FL|GB.Freq.CR|GB.Seqs.FL|GB.Seqs.CR|Homoplasmy|Heteroplasmy|Disease.Status|References|Additional.Annotations|MitoTIP">' > header.txt

        bcftools annotate ~{input_vcf} \
            -c CHROM,POS,-,REF,ALT,INFO/MITOMAP \
            --annotations subset_mitomap.tab.gz \
            --header-lines header.txt \
            --threads ~{threads} \
            --output-type z \
            --write-index=tbi \
            --output ~{sample_id}.mt.baldur.annotated.vcf.gz 
    >>>

    output {
        File out_vcf = "~{sample_id}.mt.baldur.annotated.vcf.gz"
        File out_vcf_index = "~{sample_id}.mt.baldur.annotated.vcf.gz.tbi"
    }

    runtime {
        cpu: threads
        memory: mem_gb + " GB"
        preemptible: 1
        maxRetries: 1
        docker: "quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
    }
}

task sniffles {

  input {
    String sample_id
    File input_bam
    File input_bam_index
    File mt_ref
    File mt_ref_index
    ## sniffles specific parameters
    Int min_sv_support=4

  }

  Int threads = 2
  Int mem_gb = 8

  command <<<
    set -euo pipefail

    sniffles --qc-output-all \
    --minsvlen 5 \
    --minsupport ~{min_sv_support} \
    --threads ~{threads} \
    --ref ~{mt_ref} \
    --input ~{input_bam} \
    --snf ~{sample_id}.sniffles.snf \
    --vcf ~{sample_id}.sniffles.vcf

  >>>

  output {
    File sniffles_vcf = "~{sample_id}.sniffles.vcf"
    File sniffles_snf = "~{sample_id}.sniffles.snf"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "quay.io/biocontainers/sniffles:2.6.2--pyhdfd78af_0"
  }

}

task filter_sniffles {

  input {
    String sample_id
    File sniffles_vcf
    ## sniffles specific parameters
    Int min_sv_support=4

  }

  Int threads = 1
  Int mem_gb = 4

  command <<<
    set -euo pipefail

    bcftools filter -i "SUPPORT>=~{min_sv_support} && (SVLEN>=50 || SVLEN<=-50)" ~{sniffles_vcf} > ~{sample_id}.sniffles.svlen50.vcf
    bcftools filter -i "SUPPORT>=~{min_sv_support} && ((SVLEN>=5 && SVLEN<50) || (SVLEN<=-5 && SVLEN>-50))" ~{sniffles_vcf} > ~{sample_id}.sniffles.svlen5.vcf

  >>>

  output {
    File sniffles_svlen5_vcf = "~{sample_id}.sniffles.svlen5.vcf"
    File sniffles_svlen50_vcf = "~{sample_id}.sniffles.svlen50.vcf"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "quay.io/biocontainers/bcftools:1.21--h3a4d415_1"
  }

}


task mutserve {

  input {
    String sample_id
    File input_bam
    File input_bam_index
    File mt_ref
    File mt_ref_index
  }

  Int threads = 2
  Int mem_gb = 4

  command <<<
    set -euo pipefail

    mutserve call ~{input_bam} \
    --output ~{sample_id}.mutserve.vcf.gz \
    --reference ~{mt_ref} \
    --contig-name MT \
    --threads ~{threads} --no-ansi \
    --alignQ 0 --mapQ 20 --baseQ 20 --level 0.005
  >>>

  output {
    File vcf = "~{sample_id}.mutserve.vcf.gz"
  }

  runtime {
    cpu: threads
    memory: mem_gb + " GB"
    preemptible: 1
    maxRetries: 1
    docker: "quay.io/biocontainers/mutserve:2.0.3--hdfd78af_0"
  }

}