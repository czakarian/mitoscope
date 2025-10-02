   // !!! do sep workflow for validation !!!
    // if (params.inputfile.endsWith('.bam')) {
    //     bam_ch = Channel.fromPath(params.inputfile)
    // } else if (params.inputfile.endsWith('.cram')) {
    //     cram_ch = Channel.fromPath(params.inputfile)
    // } else if (params.inputfile.endsWith('.fastq') || params.inputfile.endsWith('.fastq.gz')) {
    //     fastq_ch = Channel.fromPath(params.inputfile)
    // } else {
    //     error "Unsupported input format: ${params.inputfile}, must be .bam, .cram, .fastq, or .fastq.gz"
    // }

    // // check that aligned bam or cram has existing index file
    // if (params.is_aligned) {
    //     if (fastq_ch) {
    //         error "Fastq provided but `is_aligned = true`. Please set to false if not aligned bam or cram."
    //     } else if (bam_ch && file(params.inputfile + ".bai").exists()) {
    //         bam_index_ch = Channel.fromPath(params.inputfile + ".bai")
    //     } else if (cram_ch && file(params.inputfile + ".crai").exists()) {
    //         cram_index_ch = Channel.fromPath(params.inputfile + ".crai")
    //     } else {
    //         error "Bam or cram is missing index file. Make sure .bai or .crai file exists."
    //     }
    // }

    // // check if cram input file that ref file is provided 
    // if (cram_ch) {
    //     if (!params.reference) {
    //         error "Missing genome reference fasta (--reference) to accompany .cram input"
    //     } else if (!file(params.reference).exists()) {
    //         error "Genome reference fasta file does not exist: ${params.reference}"
    //     } else {
    //         ref_ch = Channel.fromPath(params.reference)
    //     }
    // } else {
    //     ref_ch = null
    // }