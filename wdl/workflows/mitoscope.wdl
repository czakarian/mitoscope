version 1.0

import "../tasks/format_input.wdl" as format_input
import "../tasks/mt_selection.wdl" as mt_selection
import "../tasks/assembly.wdl" as assembly
import "../tasks/variant_calling.wdl" as variant_calling
import "../tasks/numt_pipeline.wdl" as numt_pipeline
import "../tasks/quality_control.wdl" as quality_control
import "../tasks/methylation.wdl" as methylation


workflow mitoscope {
    input {
        ## user input
        String platform
        String sample_id
        File bam_file
        File bam_file_index
        File hg38_ref
        File? numt_sniffles_vcf
        Boolean run_numt_pipeline

        ### resources 
        File kmc_pre
        File kmc_suf
        File minimap_index_pb
        File minimap_index_ont
        File mt_ref
        File mt_ref_index
        File mitomap_anno_file
        File vep_cache_dir_tar
        File circos_bed
        File nuc_intervals_bed
        File script_filter_bam_py
        File script_annotate_py
        File script_rotate_py
        File script_qc_plots_py
        File script_qc_summary_py
        File script_numt_circos_py
        File script_numt_detection_sa_py
    }

    Map[String, File] minimap_index_map = {
      "ont": minimap_index_ont,
      "pb": minimap_index_pb
    }
    File minimap_index = minimap_index_map[platform]

    call format_input.aligned_bam_to_fastq {
        input:
            sample_id = sample_id,
            bam_file = bam_file,
            bam_file_index = bam_file_index
    }

    call format_input.compress_fastq {
        input:
            sample_id = sample_id,
            fastq_file = aligned_bam_to_fastq.fastq
    }

    call mt_selection.kmer_selection {
        input:
            sample_id = sample_id,
            fastq_file = compress_fastq.fastq_gz,
            kmc_pre = kmc_pre,
            kmc_suf = kmc_suf,
    }

    call mt_selection.align_to_ref {
        input:
          sample_id = sample_id,
          fastq_file = kmer_selection.kmer_selected_fastq,
          platform = platform,
          minimap_index = minimap_index

    }

    call mt_selection.sam_to_bam {
        input:
          sample_id = sample_id,
          sam = align_to_ref.sam,
    }

    call mt_selection.filter_numts {
        input:
          py_script = script_filter_bam_py,
          sample_id = sample_id,
          bam_file = sam_to_bam.bam,
          bam_file_index = sam_to_bam.bam_index
    }

    call mt_selection.filtered_bam_to_fastq {
        input:
          sample_id = sample_id,
          bam_file = filter_numts.mt_bam,
          bam_file_index = filter_numts.mt_bam_index
    }

    call assembly.mt_assembly {
        input:
          sample_id = sample_id,
          filtered_fastq = filtered_bam_to_fastq.filtered_fastq,
          platform = platform
    }


    call variant_calling.baldur {
        input:
          sample_id = sample_id,
          input_bam = filter_numts.mt_bam,
          input_bam_index = filter_numts.mt_bam_index,
          mt_ref = mt_ref,
          mt_ref_index = mt_ref_index
    }

    call variant_calling.normalize_vcf {
        input:
          sample_id = sample_id,
          baldur_vcf = baldur.baldur_vcf
    }

    call variant_calling.pull_mitomap_annos {
        input:
          py_script = script_annotate_py,
          sample_id = sample_id,
          baldur_norm_vcf = normalize_vcf.baldur_norm_vcf,
          baldur_norm_vcf_index = normalize_vcf.baldur_norm_vcf_index,
          mitomap_anno_file = mitomap_anno_file
    }

    call variant_calling.run_vep {
        input:
          sample_id = sample_id,
          baldur_norm_vcf = normalize_vcf.baldur_norm_vcf,
          baldur_norm_vcf_index = normalize_vcf.baldur_norm_vcf_index,
          mt_ref = mt_ref,
          vep_cache_dir_tar = vep_cache_dir_tar
    }

    call variant_calling.add_mitomap_to_vcf {
        input:
          sample_id = sample_id,
          input_vcf = run_vep.vep_vcf,
          input_vcf_index = run_vep.vep_vcf_index,
          mitomap_txt = pull_mitomap_annos.mitomap_txt
    }


    call variant_calling.sniffles {
        input:
          sample_id = sample_id,
          input_bam = filter_numts.mt_bam,
          input_bam_index = filter_numts.mt_bam_index,
          mt_ref = mt_ref,
          mt_ref_index = mt_ref_index
    }

    call variant_calling.filter_sniffles {
        input:
          sample_id = sample_id,
          sniffles_vcf = sniffles.sniffles_vcf,
    }

    call variant_calling.mutserve {
        input:
          sample_id = sample_id,
          input_bam = filter_numts.mt_bam,
          input_bam_index = filter_numts.mt_bam_index,
          mt_ref = mt_ref,
          mt_ref_index = mt_ref_index
    }

    call quality_control.mt_coverage {
        input:
          sample_id = sample_id,
          input_bam = filter_numts.mt_bam,
          input_bam_index = filter_numts.mt_bam_index
    }

    call quality_control.nuclear_coverage {
        input:
          sample_id = sample_id,
          input_bam = bam_file,
          input_bam_index = bam_file_index,
          ref = hg38_ref,
          nuc_intervals_bed = nuc_intervals_bed
    }

    call quality_control.mt_read_length {
        input:
          sample_id = sample_id,
          input_fastq = filtered_bam_to_fastq.filtered_fastq
    }

    call quality_control.coverage_plot {
        input:
          py_script = script_qc_plots_py,
          sample_id = sample_id,
          per_base_bed = mt_coverage.per_base_bed
    }

    call quality_control.read_length_plot {
        input:
          py_script = script_qc_plots_py,
          sample_id = sample_id,
          read_length_file = mt_read_length.read_length_txt
    }

    call quality_control.haplogrep {
        input:
          sample_id = sample_id,
          mutserve_vcf = mutserve.vcf
    }

    call quality_control.haplocheck {
        input:
          sample_id = sample_id,
          mutserve_vcf = mutserve.vcf
    }

    call quality_control.qc_summary {
        input:
          py_script = script_qc_summary_py,
          sample_id = sample_id,
          mito_coverage_file = mt_coverage.mosdepth_summary,
          nuclear_coverage_file = nuclear_coverage.mosdepth_summary,
          read_length_file = mt_read_length.read_length_txt,
          haplogrep_file = haplogrep.haplogrep_txt,
          kmer_read_count = kmer_selection.kmer_read_count,
          chrM_read_count = aligned_bam_to_fastq.chrM_read_count
    }

    if (run_numt_pipeline == true) {

        call numt_pipeline.make_rotated_mt_ref {
            input:
              py_script = script_rotate_py,
              mt_ref = mt_ref
        }

        call numt_pipeline.make_blast_dbs {
            input:
              hg38_ref = hg38_ref,
              mt_ref = mt_ref,
              mt_ref_rotated = make_rotated_mt_ref.mt_ref_rotated
        }

        if (!defined(numt_sniffles_vcf)) {
          call numt_pipeline.numt_detection_sniffles {
              input:
                sample_id = sample_id,
                input_bam = bam_file,
                input_bam_index = bam_file_index,
                ref = hg38_ref
          }
        }

        call numt_pipeline.numt_detection_insertions_to_fasta {
            input:
              sample_id = sample_id,
              input_vcf = select_first([numt_sniffles_vcf, numt_detection_sniffles.vcf])
        }

        call numt_pipeline.numt_detection_insertions_blast {
            input:
              sample_id = sample_id,
              input_fasta = numt_detection_insertions_to_fasta.fasta,
              blast_db_tar = make_blast_dbs.blast_db
        }

        call numt_pipeline.numt_detection_mtdna_insertions_to_fasta {
            input:
              sample_id = sample_id,
              input_fasta = numt_detection_insertions_to_fasta.fasta,
              blast_output = numt_detection_insertions_blast.blast_output
        }

        call numt_pipeline.numt_detection_mtdna_insertions_blast_check {
            input:
              sample_id = sample_id,
              input_fasta = numt_detection_mtdna_insertions_to_fasta.fasta,
              blast_db_tar = make_blast_dbs.blast_db
        }

        call numt_pipeline.numt_detection_plot {
            input:
              py_script = script_numt_circos_py,
              sample_id = sample_id,
              blast_file = numt_detection_mtdna_insertions_blast_check.blast_output,
              circos_bed = circos_bed,
        }

        call numt_pipeline.numt_detection_supplementary {
            input:
              py_script = script_numt_detection_sa_py,
              sample_id = sample_id,
              input_bam = bam_file,
              input_bam_index = bam_file_index,
              hg38_ref = hg38_ref
        }
    }
    if (platform == 'pb') {
      call methylation.meth_freq_pb {
          input:
            sample_id = sample_id,
            input_bam = filter_numts.mt_bam,
            input_bam_index = filter_numts.mt_bam_index,
            mt_ref = mt_ref,
            mt_ref_index = mt_ref_index
      }
    }
    if (platform == 'ont') {
        call methylation.meth_freq_ont {
          input:
            sample_id = sample_id,
            input_bam = filter_numts.mt_bam,
            input_bam_index = filter_numts.mt_bam_index,
            mt_ref = mt_ref,
            mt_ref_index = mt_ref_index
      }
    }

    output {
      File filtered_bam = filter_numts.mt_bam
      File filtered_bam_index = filter_numts.mt_bam_index
      File numt_bam = filter_numts.discard_bam
      File numt_bam_index = filter_numts.discard_bam_index
      Array[File] assembly_output = mt_assembly.mt_assembly_dir
      File anno_baldur_vcf = add_mitomap_to_vcf.out_vcf
      File anno_baldur_vcf_index = add_mitomap_to_vcf.out_vcf_index
      File? baldur_dels = baldur.baldur_dels
      File sniffles_svlen5_vcf = filter_sniffles.sniffles_svlen5_vcf
      File sniffles_svlen50_vcf = filter_sniffles.sniffles_svlen50_vcf
      File qc_table = qc_summary.tsv
      File? methylation_ont = meth_freq_ont.bedmethyl_CG
      File? methylation_pb= meth_freq_pb.bed
      File? numt_insertion_output = numt_detection_mtdna_insertions_blast_check.tsv_output
      File? numt_insertion_plot = numt_detection_plot.circos_plot
      File? numt_supplementary_output = numt_detection_supplementary.supplementary_numt_tsv

    }

}