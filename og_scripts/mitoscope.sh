#!/bin/bash
######################################################################
# Copyright (C) 2019-2024 Dr Chia-Lin Wei Laboratory, NWGC
# Adapted from ecLegoV3 for human mitochondria
#
# PENDING:
# C2. <INV>, <DEL>, <INS> required different minimum read support
# C3. address <DUP> at sub-genomic level
# C4. visualization:
#     circos or CGView (https://js.cgview.ca/examples/index.html)
# C5. population based analysis
# C6. visualization at population-level (inter-samples)
# C8. modifications analysis (on non-PCR library construction)
#
# PENDING (technical):
# T1. Convert to nextflow for parallel execution of non-dependent 
#     subjobs
#
# NOTE:
# N1. .gv visualization with https://dreampuf.github.io/GraphvizOnline/
#
#
######################################################################


set -eu -o pipefail

usage() {
    echo "Usage: mitoscope.sh -i <input_fastq> -p <ont|pb> [-t <threads>] [-m <minreadsupport>] [-c <cncov>]"
    echo
    echo "Required arguments:"
    echo "  -i <input>            Input FASTQ or BAM file."
    echo "  -o <outputdir>        Output directory."
    echo "  -p <ont|pb>           Platform: 'ont' for Oxford Nanopore or 'pb' for PacBio."
    echo
    echo "Optional arguments:"
    echo "  -t <threads>          Number of threads to use (default: 4)."
    echo "  -m <minreadsupport>   Minimum read support for SV calling (default: 2)."
    echo "  -c <cncov>            CN Coverage? (default: 12)."
}

# Check if no arguments were passed
if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

# Default values
export THREADS=4  
export MINREADSUPPORT=2
export CNCOV=""

while getopts "i:o:p:t:m:c:h" FLAG; do
    case ${FLAG} in
        i) INPUTFILE=${OPTARG};;
        o) OUTDIR=${OPTARG};;
        p) PLATFORM=${OPTARG};;
        t) THREADS=${OPTARG};;
        m) MINREADSUPPORT=${OPTARG};;
        c) CNCOV=${OPTARG};;
        h) 
            usage
            exit 0
            ;;
        *) 
            echo "Invalid or missing arguments."
            usage
            exit 1
            ;;
    esac
done

# Check if required arguments are set
if [[ -z "${INPUTFILE}"  ||  -z "${OUTDIR}"  ||  -z "${PLATFORM}" ]]; then
    echo "Error: -i <input>, -o <outputdir>, and -p <ont/pb> options are all required."
    exit 1
fi

if [[ "${PLATFORM}" != "ont"  &&  "${PLATFORM}" != "pb" ]]; then
    echo "Invalid option for -p: ${PLATFORM}. Must be 'ont' or 'pb'."
    exit 1
fi

if [[ ! -f "${INPUTFILE}" ]]; then
    echo "${INPUTFILE} does not exist"
    exit 2
fi

if [[ "${INPUTFILE}" != *.bam  && "${INPUTFILE}" != *.fastq  && "${INPUTFILE}" != *.fastq.gz ]]; then
    echo "File ${INPUTFILE} does not have .bam or .fastq/.fastq.gz extension."
    exit 2
fi

## set ont/hifi platform variables
if [ ${PLATFORM} == "ont" ]; then
    export FLYEPLATFORM="--nano-hq"
    export MINIMAPPLATFORM="map-ont"
    export MINIMAPINDEX="MT-ont.mmi"
elif [ ${PLATFORM} == "pb" ]; then
    export FLYEPLATFORM="--pacbio-hifi "
    export MINIMAPPLATFORM="map-hifi"
    export MINIMAPINDEX="MT-hifi.mmi"
fi

echo "# Input file = ${INPUTFILE}"
echo "# Output directory = ${OUTDIR}"
echo "# Threads = ${THREADS}"
echo "# Minimum read support = ${MINREADSUPPORT}"
covCN=12
if [[ -n "${CNCOV}" ]]; then
    covcn_value=$(echo "${CNCOV}" | bc)
    if [[ $covcn_value -gt 0 ]]; then
        echo "# Adjusting covCN ${covCN} to user-defined value ${covcn_value}"
        covCN="${covcn_value}"
    else
        echo "#####"
        echo "# WARNING: Ignoring unusable coverage threshold provided ${CNCOV}"
        echo "#          Using covCN ${covCN} instead"
        echo "#####"
    fi
else
    echo "# Assumed 25x coverage, using covCN ${covCN}"
fi


export MITOSCOPE_ROOT="$(dirname "$(readlink -f "${BASH_SOURCE}")")"
export MITOSCOPE_RESOURCES="${MITOSCOPE_ROOT}/resources"
export MITOSCOPE_SINGULARITY="${MITOSCOPE_ROOT}/singularity"
export MITOSCOPE_TOOLS="${MITOSCOPE_ROOT}/tools"
export KMCTOOLSCMD="${MITOSCOPE_SINGULARITY}/kmc_3.2.1.sif kmc_tools"
export FLYECMD="${MITOSCOPE_SINGULARITY}/flye_2.9.6.sif flye"
export MINIMAP2CMD="${MITOSCOPE_SINGULARITY}/minimap2_2.24.sif minimap2"
export SAMTOOLSCMD="${MITOSCOPE_SINGULARITY}/samtools_1.21.sif samtools"
export BCFTOOLSCMD="${MITOSCOPE_SINGULARITY}/bcftools_1.21.sif bcftools"
export GENOMECOVERAGEBEDCMD="${MITOSCOPE_SINGULARITY}/bedtools_2.31.0.sif genomeCoverageBed"
export SORTBEDCMD="${MITOSCOPE_SINGULARITY}/bedtools_2.31.0.sif sortBed"
export BG2BWCMD="${MITOSCOPE_SINGULARITY}/ucsc-bedgraphtobigwig_445.sif bedGraphToBigWig"
export SNIFFLESCMD="${MITOSCOPE_SINGULARITY}/sniffles_2.6.2.sif sniffles"
#export GATKCMD="${MITOSCOPE_SINGULARITY}/gatk_4.5.0.0.sif gatk"
export MUTSERVECMD="${MITOSCOPE_SINGULARITY}/mutserve_2.0.3.sif mutserve"
export HAPLOGREPCMD="${MITOSCOPE_SINGULARITY}/haplogrep_3.2.2.sif haplogrep3"
export HAPLOCHECKCMD="${MITOSCOPE_SINGULARITY}/haplocheck_1.3.3.sif haplocheck"
export BALDURCMD="${MITOSCOPE_SINGULARITY}/baldur_1.2.2.sif baldur"
export MOSDEPTHCMD="${MITOSCOPE_SINGULARITY}/mosdepth_0.3.8.sif mosdepth"
export MINIMODCMD="${MITOSCOPE_SINGULARITY}/minimod_0.4.0.sif minimod"
#export MINIMODCMD="${MITOSCOPE_TOOLS}/minimod/minimod-v0.4.0/minimod"


INPUTFILE=$(readlink -f "${INPUTFILE}") ## resolve the symbolic link to actual file path
export INPUTDIR=$(dirname "${INPUTFILE}")
OUTDIR="${OUTDIR%/}"
[ ! -d "${OUTDIR}" ] && mkdir "${OUTDIR}"
export RESULTDIR=${OUTDIR}/mitoscope
export ALIGNDIR=${OUTDIR}/mitoscope/alignments
export VARIANTDIR=${OUTDIR}/mitoscope/variants
export QCDIR=${OUTDIR}/mitoscope/qc
mkdir -p "${RESULTDIR}"
mkdir -p "${ALIGNDIR}"
mkdir -p "${VARIANTDIR}"
mkdir -p "${QCDIR}"

export FASTQPREFIX="$(basename "${INPUTFILE}")"
export FASTQPREFIX="${FASTQPREFIX%.gz}"
export FASTQPREFIX="${FASTQPREFIX%.*}"

export SINGULARITY_BINDPATH="${MITOSCOPE_ROOT},${INPUTDIR}"
if [[ "${OUTDIR}" != "${INPUTDIR}" ]]; then
    export SINGULARITY_BINDPATH="${MITOSCOPE_ROOT},${INPUTDIR},${OUTDIR}"
fi

export KMCTOOLSTHREADS=${THREADS}
export PIGZTHREADS=${THREADS}
export MINIMAP2THREADS=${THREADS}
export SAMTOOLSTHREADS=${THREADS}
export FLYETHREADS=${THREADS}
export MUTSERVETHREADS=${THREADS}
export MUTECTTHREADS=${THREADS}
export MINIMODTHREADS=${THREADS}

export FLYEMINOVERLAP=2500

# if input file is a bam convert to fastq
if [[ "${INPUTFILE}" == *.fastq || "${INPUTFILE}" == *.fastq.gz ]]; then
    echo '==' $(date) '==' Fastq input provided, bam to fastq conversion SKIPPED
    export FASTQ=${INPUTFILE}
elif [[ "${INPUTFILE}" == *.bam ]]; then
    echo '==' $(date) '==' Bam to fastq conversion STARTED
    export FASTQ=${OUTDIR}/${FASTQPREFIX}.fastq.gz
    ${SAMTOOLSCMD} fastq -T MM,ML -@ ${SAMTOOLSTHREADS} ${INPUTFILE} | tr '\t' ' ' | pigz -p ${PIGZTHREADS} - > ${FASTQ}
    echo '==' $(date) '==' Bam to fastq conversion COMPLETE
fi
#

## select long-reads which are likely from MT
echo '==' $(date) '==' MT candidate fastq generation STARTED
${KMCTOOLSCMD} -t${KMCTOOLSTHREADS} filter ${MITOSCOPE_RESOURCES}/MT.k29 -ci1 ${FASTQ} -fq -ci2500 ${RESULTDIR}/${FASTQPREFIX}.MT.fastq 
cat ${RESULTDIR}/${FASTQPREFIX}.MT.fastq | tr ' ' '\t' | pigz -p ${PIGZTHREADS} - > ${RESULTDIR}/${FASTQPREFIX}.MT.fastq.gz
rm ${RESULTDIR}/${FASTQPREFIX}.MT.fastq 
echo '==' $(date) '==' MT candidate fastq generation ENDED
#

# align MT reads to ref for read filtration + debugging
echo '==' $(date) '==' MT candidates fastq reference mapping STARTED
${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -Y -y -t ${MINIMAP2THREADS} ${MITOSCOPE_RESOURCES}/${MINIMAPINDEX} ${RESULTDIR}/${FASTQPREFIX}.MT.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.bam 
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.bam 
echo '==' $(date) '==' MT candidates fastq reference mapping COMPLETED
#

# remove foldback reads and NUMTs 
echo '==' $(date) '==' Removal of foldback MT candidates STARTED
python ${MITOSCOPE_ROOT}/filter_bam.py -i ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.bam 
echo '==' $(date) '==' Removal of foldback MT candidates COMPLETED
# 

### some qc steps // coverage / read lengths / n50s / methylation
echo '==' $(date) '==' Get basic QC stats -- methylation, mean coverage, read count, avg read length, n50 STARTED

# plot average methylation levels across MT genome 
${MINIMODCMD} freq -t ${MINIMODTHREADS} -m 0.5 -o ${QCDIR}/${FASTQPREFIX}.minimod.tsv ${MITOSCOPE_RESOURCES}/MT.fasta ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam 
sort -k2 -n ${QCDIR}/${FASTQPREFIX}.minimod.tsv -o ${QCDIR}/${FASTQPREFIX}.minimod.tsv
${MITOSCOPE_ROOT}/qc_plots.py --plot methylation --input ${QCDIR}/${FASTQPREFIX}.minimod.tsv --outprefix ${QCDIR}/${FASTQPREFIX}
tail +2 ${QCDIR}/${FASTQPREFIX}.minimod.tsv | cut -f 1-3,7  > ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.minimod.bedGraph

## plot mean coverage across MT genome
${MOSDEPTHCMD} ${QCDIR}/${FASTQPREFIX} ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam
${MITOSCOPE_ROOT}/qc_plots.py --plot coverage --input ${QCDIR}/${FASTQPREFIX}.per-base.bed.gz --outprefix ${QCDIR}/${FASTQPREFIX}

## plot avg read length
${SAMTOOLSCMD} view ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam | awk '{print length($10)}' > ${QCDIR}/${FASTQPREFIX}.read_lengths.txt
${MITOSCOPE_ROOT}/qc_plots.py --plot read_length --input ${QCDIR}/${FASTQPREFIX}.read_lengths.txt --outprefix ${QCDIR}/${FASTQPREFIX}

## get overall qc stats and to qc_summary.txt
mean_cov=$(awk '$1 == "MT" {print $4}' ${QCDIR}/${FASTQPREFIX}.mosdepth.summary.txt)
read_count=$(cat ${QCDIR}/${FASTQPREFIX}.read_lengths.txt | wc -l)
avg_length=$(awk '{sum+=$1} END {print sum/NR}' ${QCDIR}/${FASTQPREFIX}.read_lengths.txt)

## calculate N50
n50=$(awk '{sum+=$1; arr[NR]=$1} END {
    half = sum/2;
    n = asort(arr);
    total = 0;
    for (i = n; i >= 1; i--) {
        total += arr[i];
        if (total >= half) {
            print arr[i];
            break;
        }
    }
}' ${QCDIR}/${FASTQPREFIX}.read_lengths.txt)

# get methylation stats
avg_meth=$(awk '{sum+=$7} END {print sum/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)
meth_gt_1=$(awk 'NR>1 && $7 > 0.01 {count++} END {print count/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)
meth_gt_5=$(awk 'NR>1 && $7 > 0.05 {count++} END {print count/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)
meth_gt_10=$(awk 'NR>1 && $7 > 0.1 {count++} END {print count/(NR-1)}' ${QCDIR}/${FASTQPREFIX}.minimod.tsv)

echo -e "Sample\tRead_Count\tMean_Coverage\tAverage_Read_Length\tN50\tAverage_Methylation_Percent\tNum_Meth_Sites_GT_1\tNum_Meth_Sites_GT_5\tNum_Meth_Sites_GT_10" > ${QCDIR}/${FASTQPREFIX}.qc_summary.txt
echo -e "${FASTQPREFIX}\t${read_count}\t${mean_cov}\t${avg_length}\t${n50}\t${avg_meth}\t${meth_gt_1}\t${meth_gt_5}\t${meth_gt_10}" >> ${QCDIR}/${FASTQPREFIX}.qc_summary.txt

echo '==' $(date) '==' Get basic QC stats -- methylation, mean coverage, read count, avg read length, n50 COMPLETED
##  

# baldur variant calls (SNV, small indel, large deletion)
echo '==' $(date) '==' Baldur variant calling START
mkdir -p ${VARIANTDIR}/baldur

${BALDURCMD} -n ${FASTQPREFIX} -l debug --output-deletions -T ${MITOSCOPE_RESOURCES}/MT.fasta \
-o ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam

${BCFTOOLSCMD} norm --multiallelics -both ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.vcf.gz | \
${BCFTOOLSCMD} norm --atomize | ${BCFTOOLSCMD} view -f PASS -Oz -o ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.vcf.gz 
${BCFTOOLSCMD} index --tbi ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.vcf.gz 

${BCFTOOLSCMD} view --types indels ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.vcf.gz -Oz -o ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.indels.vcf.gz 
${BCFTOOLSCMD} view --exclude-types indels ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.vcf.gz -Oz -o ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.snvs.vcf.gz 
${BCFTOOLSCMD} index --tbi ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.indels.vcf.gz 
${BCFTOOLSCMD} index --tbi ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.snvs.vcf.gz 

echo '==' $(date) '==' Baldur variant calling COMPLETED
#

# add MITOMAP annotations to baldur output 
echo '==' $(date) '==' MITOMAP annotation of baldur output STARTED
python ${MITOSCOPE_ROOT}/annotate.py \
--input ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.indels.vcf.gz \
--annotations ${MITOSCOPE_ROOT}/annotations/CombinedDiseaseVariantDB.csv \
--caller baldur 

python ${MITOSCOPE_ROOT}/annotate.py \
--input ${VARIANTDIR}/baldur/${FASTQPREFIX}.MT.ref.filtered.baldur.norm.snvs.vcf.gz \
--annotations ${MITOSCOPE_ROOT}/annotations/CombinedDiseaseVariantDB.csv \
--caller baldur 
echo '==' $(date) '==' MITOMAP annotation of baldur output COMPLETED
#

# call SNVs using mutserve 
echo '==' $(date) '==' Mutserve SNV calling STARTED
mkdir -p ${VARIANTDIR}/mutserve

## move into mutserve directory (so .txt file gets outputted in there)
workdir=$(pwd)
cd ${VARIANTDIR}/mutserve

${MUTSERVECMD} call ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam \
--output ${FASTQPREFIX}.MT.ref.filtered.mutserve.vcf.gz \
--reference ${MITOSCOPE_RESOURCES}/MT.fasta \
--threads ${MUTSERVETHREADS} --no-ansi

## rename txt output from mutserve since it always truncates full name
#mv ${VARIANTDIR}/mutserve/${FASTQPREFIX%%.*}.txt ${VARIANTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.txt
cd ${workdir}

# ## normalize multiallelics
${BCFTOOLSCMD} norm --multiallelics -both ${VARIANTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.vcf.gz | \
${BCFTOOLSCMD} norm --atomize | ${BCFTOOLSCMD} view -f PASS -Oz -o ${VARIANTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.norm.vcf.gz
${BCFTOOLSCMD} index --tbi ${VARIANTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.norm.vcf.gz

#
echo '==' $(date) '==' Mutserve SNV calling COMPLETED

# add MITOMAP annotations to mutserve output 
echo '==' $(date) '==' MITOMAP annotation of mutserve output STARTED
python ${MITOSCOPE_ROOT}/annotate.py \
--input ${VARIANTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.norm.vcf.gz \
--annotations ${MITOSCOPE_ROOT}/annotations/CombinedDiseaseVariantDB.csv \
--caller mutserve 
echo '==' $(date) '==' MITOMAP annotation of mutserve output COMPLETED
#

# haplogroup classification using haplogrep3
mkdir -p ${VARIANTDIR}/haplogrep

echo '==' $(date) '==' Haplogroup classification STARTED
TREE="phylotree-rcrs@17.2"
${HAPLOGREPCMD} classify \
--tree ${TREE} \
--input ${VARIANTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.vcf.gz \
--output ${VARIANTDIR}/haplogrep/${FASTQPREFIX}.MT.ref.filtered.haplogrep.txt \
--extend-report --write-qc
echo '==' $(date) '==' Haplogroup classification COMPLETED
#

# haplocheck contamination check
echo '==' $(date) '==' Haplocheck contamination check STARTED
${HAPLOCHECKCMD} --raw --out ${VARIANTDIR}/haplogrep/${FASTQPREFIX}.MT.ref.filtered.haplocheck.txt \
${VARIANTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.vcf.gz 
echo '==' $(date) '==' Haplocheck contamination check COMPLETED
#

# call SVs using sniffles
echo '==' $(date) '==' Filtered MT candidates fastq variation against reference STARTED
${SNIFFLESCMD} --qc-output-all --allow-overwrite \
--minsupport ${MINREADSUPPORT} \
--input ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam \
--vcf ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam.raw.vcf

${BCFTOOLSCMD} filter -i "SUPPORT>=${MINREADSUPPORT}" ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam.raw.vcf > ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam.raw.ge${MINREADSUPPORT}.vcf
echo '==' $(date) '==' Filtered MT candidates fastq variation against reference COMPLETED
#

# assemble the selected long-reads
echo '==' $(date) '==' Bam to fastq conversion for NUMT/foldback filtered bam STARTED
${SAMTOOLSCMD} fastq -T MM,ML -@ ${SAMTOOLSTHREADS} -0 ${RESULTDIR}/${FASTQPREFIX}.MT.filtered.fastq.gz ${ALIGNDIR}/${FASTQPREFIX}.MT.ref.filtered.bam 
echo '==' $(date) '==' Bam to fastq conversion for NUMT/foldback filtered bam COMPLETED

echo '==' $(date) '==' de Novo MT assembly STARTED
${FLYECMD} --threads ${FLYETHREADS} --meta ${FLYEPLATFORM} \
${RESULTDIR}/${FASTQPREFIX}.MT.filtered.fastq.gz --out-dir ${RESULTDIR}/MT_assembly -m ${FLYEMINOVERLAP}
echo '==' $(date) '==' de Novo MT assembly COMPLETED
#

# map selected long-reads to assembled contig(s)
echo '==' $(date) '==' MT candidates fastq assembly mapping STARTED
(${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -Y -y -t ${MINIMAP2THREADS} ${RESULTDIR}/MT_assembly/assembly.fasta ${RESULTDIR}/${FASTQPREFIX}.MT.filtered.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam ; \
export DSNAME=${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} | ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${SAMTOOLSCMD} faidx ${RESULTDIR}/MT_assembly/assembly.fasta ; \
${BG2BWCMD} ${DSNAME}.bg ${RESULTDIR}/MT_assembly/assembly.fasta.fai ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' MT candidates fastq assembly mapping COMPLETED
#

## should i adjust the assembly to match ref coordinates here? then i can do sniffles + other variant calling?

# call variations of selected long-reads w.r.t. assembled contig(s)
echo '==' $(date) '==' MT candidates fastq variation against assembly STARTED
${SNIFFLESCMD} --qc-output-all --allow-overwrite \
--minsupport ${MINREADSUPPORT} \
--input ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam \
--vcf ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf

${BCFTOOLSCMD} filter -i "SUPPORT>=${MINREADSUPPORT}" ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf \
> ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.ge${MINREADSUPPORT}.vcf
echo '==' $(date) '==' MT candidates fastq variation against assembly COMPLETED
#

# for inter-sample anchoring + debugging
echo '==' $(date) '==' Assembly reference mapping STARTED
(${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -Y -t ${MINIMAP2THREADS} ${MITOSCOPE_RESOURCES}/${MINIMAPINDEX} ${RESULTDIR}/MT_assembly/assembly.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.ref.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.ref.bam ; \
export DSNAME=${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.ref.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} | ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg ${MITOSCOPE_RESOURCES}/MT.fasta.fai ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' Assembly reference mapping COMPLETED
#

# for inter-sample anchoring + debugging
echo '==' $(date) '==' Assembly variation against reference STARTED
${SNIFFLESCMD} --qc-output-all --allow-overwrite \
--input ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.ref.bam \
--vcf ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.ref.bam.raw.vcf
echo '==' $(date) '==' Assembly variation against reference COMPLETED
#


# # for debugging
# echo '==' $(date) '==' Graph_before reference mapping STARTED
# (${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -Y -t ${MINIMAP2THREADS} ${MITOSCOPE_RESOURCES}/${MINIMAPINDEX} ${RESULTDIR}/MT_assembly/20-repeat/graph_before_rr.fasta \
# | ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${ALIGNDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam ; \
# ${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${ALIGNDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam ; \
# export DSNAME=${ALIGNDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam; \
# export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
# echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
# ${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} | ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
# echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
# ${BG2BWCMD} ${DSNAME}.bg ${MITOSCOPE_RESOURCES}/MT.fasta.fai ${DSNAME}.bw && rm ${DSNAME}.bg ; \
# echo '==' $(date) '==' Done. ) &
# wait
# echo '==' $(date) '==' Graph_before fastq reference mapping COMPLETED
# #


# generate subpopulation
# TODO: determine at what threshold we claim a new subpopulation
echo '==' $(date) '==' Subpopulations generation STARTED
export SIEVEDGRAPHDIR=${RESULTDIR}/MT_assembly/sieved_graph
[ ! -d "${SIEVEDGRAPHDIR}" ] && mkdir -p "${SIEVEDGRAPHDIR}"
pushd "${SIEVEDGRAPHDIR}"
# 

echo '==' $(date) '==' check assembly graph for circular genome..
perl ${MITOSCOPE_ROOT}/ecLegov2.pl sievegraph \
--gv ../assembly_graph.gv \
--diploidcov $(expr 2 \* ${covCN}) \
--oprefix ${FASTQPREFIX}.MT.assembly \
--bam ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.ref.bam

# FIXME: lower dicncov to 4 for SV threshold as >=2
echo '==' $(date) '==' subpopulation script setup..
perl ${MITOSCOPE_ROOT}/cgSupPop.pl setup \
--assembly ${FASTQPREFIX}.MT.assembly.cn$(expr 2 \* ${covCN}).disjointcyclic.gv.overview.xls \
--vcf ${ALIGNDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf \
--dicncov 4 --subgraph 1 --rsg 1 --sample ${FASTQPREFIX} \
--scriptPath ${MITOSCOPE_ROOT}
#
mv process.sh process.sh.tmp
cat process.sh.tmp | sed s?t2tv2.fasta?MT.fasta? > process.sh
chmod u+x *.sh
echo '==' $(date) '==' writing possible subpopulations..
./setup.sh
echo '==' $(date) '==' variation calling on possible subpopulations..
./process.sh
echo '==' $(date) '==' refining alignments possible subpopulations..
./refine.sh
#
popd
#
NUMOFCONTIGS=$(expr $(cat ${RESULTDIR}/MT_assembly/assembly_info.txt | wc -l) - 1)
if [[ -n "${NUMOFCONTIGS}" ]]; then
    NUMOFCONTIGS_value=$(echo "${NUMOFCONTIGS}" | bc)
    if [[ $NUMOFCONTIGS_value -gt 1 ]]; then
        echo ""
        echo "#####"
        echo "# WARNING: Multiple contigs (${NUMOFCONTIGS}) assembled."
        echo "#          Please check the assembly results for major structural variations."
        echo "#####"
    else
        echo ""
        echo "#####"
        echo "# INFO: Single assembled contig."
        echo "#####"
    fi
else
    echo ""
    echo "#####"
    echo "# WARNING: Cannot determine the number of assembled contigs"
    echo "#          Please check the assembly result."
    echo "#####"
fi
echo '==' $(date) '==' Subpopulations generation COMPLETED
#
