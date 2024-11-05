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
# C7. support .bam as input or propagate MM,ML-tags - for modications
# C8. modifications analysis (on non-PCR library construction)
#
# PENDING (technical):
# T1. Convert to nextflow for parallel execution of non-dependent 
#     subjobs
# T2. Add command-line options with getopts
# T3. add bam to fastq conversion into script?
#
# NOTE:
# N1. .gv visualization with https://dreampuf.github.io/GraphvizOnline/
#
#
######################################################################


set -eu -o pipefail

usage() {
    echo "Usage: mitoscope.sh -f <input_fastq> -p <ont|pb> [-t <threads>] [-m <minreadsupport>] [-c <cncov>]"
    echo
    echo "Required arguments:"
    echo "  -f <input_fastq>      Input FASTQ file."
    echo "  -p <ont|pb>           Technology type: 'ont' for Oxford Nanopore or 'pb' for PacBio."
    echo
    echo "Optional arguments:"
    echo "  -t <threads>          Number of threads to use (default: 4)."
    echo "  -m <minreadsupport>   Minimum read support for SV calling (default: 2)."
    echo "  -c <cncov>            Coverage? (default: auto)??."
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

while getopts "f:p:t:m:c:h" FLAG; do
    case ${FLAG} in
        f) FASTQ=${OPTARG};;
        p) TECH=${OPTARG};;
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
if [ -z "$FASTQ" ] || [ -z "$TECH" ]; then
    echo "Error: -f and -p options are required."
    exit 1
fi

if [[ "$TECH" != "ont" && "$TECH" != "pb" ]]; then
    echo "Invalid option for -p: $TECH. Must be 'ont' or 'pb'."
    exit 1
fi

[ -z "${FASTQ}" ] && { echo "Fastq file is empty" ; exit 1; }
[ ! -f "${FASTQ}" ] && { echo "${FASTQ} does not exist" ; exit 2 ; }

## set ont/hifi platform variables
if [ "$TECH" == "ont" ]; then
    export FLYEPLATFORM="--nano-hq"
    export MINIMAPPLATFORM="map-ont"
    export MINIMAPINDEX="MT-ont.mmi"
elif [ "$TECH" == "pb" ]; then
    export FLYEPLATFORM="--pacbio-hifi "
    export MINIMAPPLATFORM="map-hifi"
    export MINIMAPINDEX="MT-hifi.mmi"
fi

echo "# FASTQ = ${FASTQ}"
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


# export MITOSCOPE_ROOT="$(dirname "$(readlink -f "${BASH_SOURCE}")")"
# TODO: set path to your copy of mitoscope
export MITOSCOPE_ROOT="/net/nwgc/vol1/home/czaka/tools/mitoscope"

export MITOSCOPE_RESOURCES="${MITOSCOPE_ROOT}/resources"
export MITOSCOPE_SINGULARITY="${MITOSCOPE_ROOT}/singularity"
export MITOSCOPE_TOOLS="${MITOSCOPE_ROOT}/tools"
export KMCTOOLSCMD="${MITOSCOPE_SINGULARITY}/kmc_3.2.1.sif kmc_tools"
export FLYECMD="${MITOSCOPE_SINGULARITY}/flye_2.9.1.sif flye"
export MINIMAP2CMD="${MITOSCOPE_SINGULARITY}/minimap2_2.24.sif minimap2"
export SAMTOOLSCMD="${MITOSCOPE_SINGULARITY}/samtools_v1.15.1.sif samtools"
export GENOMECOVERAGEBEDCMD="${MITOSCOPE_SINGULARITY}/bedtools_2.31.0.sif genomeCoverageBed"
export SORTBEDCMD="${MITOSCOPE_SINGULARITY}/bedtools_2.31.0.sif sortBed"
export BG2BWCMD="${MITOSCOPE_SINGULARITY}/ucsc-bedgraphtobigwig_445.sif bedGraphToBigWig"
export SNIFFLESCMD="${MITOSCOPE_SINGULARITY}/sniffles_2.3.3.sif sniffles"
export BCFTOOLSCMD="${MITOSCOPE_SINGULARITY}/bcftools_1.19.sif bcftools"
export GATKCMD="${MITOSCOPE_SINGULARITY}/gatk_4.5.0.0.sif gatk"
export MUTSERVECMD="${MITOSCOPE_TOOLS}/mutserve_2.0.1/mutserve"
export HAPLOGREPCMD="${MITOSCOPE_TOOLS}/haplogrep_3.2.2/haplogrep3"
export HAPLOCHECKCMD="${MITOSCOPE_TOOLS}/haplocheck_1.3.3/haplocheck"

## temp use R module for annotate.R 
module load modules modules-init modules-gs
module load R/4.3.2

export FASTQDIR="$(dirname "$(readlink -f "${FASTQ}")")"
export RESULTDIR=${FASTQDIR}/mitoscope
[ ! -d "${RESULTDIR}" ] && mkdir "${RESULTDIR}"
export DEBUGDIR=${FASTQDIR}/mitoscope/debug
[ ! -d "${DEBUGDIR}" ] && mkdir "${DEBUGDIR}"

export FASTQPREFIX="$(basename "$(readlink -f "${FASTQ}")")"
export FASTQPREFIX="${FASTQPREFIX%.gz}"
export FASTQPREFIX="${FASTQPREFIX%.*}"

export SINGULARITY_BINDPATH="${MITOSCOPE_ROOT},${FASTQDIR}"

export KMCTOOLSTHREADS=${THREADS}
export PIGZTHREADS=${THREADS}
export MINIMAP2THREADS=${THREADS}
export SAMTOOLSTHREADS=${THREADS}
export FLYETHREADS=${THREADS}
export MUTSERVETHREADS=${THREADS}
export MUTECTTHREADS=${THREADS}

export FLYEMINOVERLAP=2500

# select long-reads which are likely from MT
echo '==' $(date) '==' MT candidate fastq generation STARTED
time (${KMCTOOLSCMD} -t${KMCTOOLSTHREADS} filter \
${MITOSCOPE_RESOURCES}/MT.k29 -ci1 ${FASTQ} \
-fq -ci2500 ${RESULTDIR}/${FASTQPREFIX}.MT.all.fastq) &
wait
echo '==' $(date) '==' MT candidate fastq generation ENDED
echo '==' $(date) '==' MT oversized candidate removal STARTED
(cat ${RESULTDIR}/${FASTQPREFIX}.MT.all.fastq \
| perl -e 'while ($L1=<STDIN>) { $L2=<STDIN>; $L3=<STDIN>; $L4=<STDIN>;
  chomp($L2); if (length($L2)<=16569) {
    print $L1; print $L2,"\n"; print $L3; print $L4;
  }
}' \
| pigz -p ${PIGZTHREADS} - > ${RESULTDIR}/${FASTQPREFIX}.MT.fastq.gz) &
wait
echo '==' $(date) '==' MT oversized candidate removal ENDED
echo '==' $(date) '==' MT all candidate fastq compression STARTED
time (pigz -p ${PIGZTHREADS} ${RESULTDIR}/${FASTQPREFIX}.MT.all.fastq) &
wait
echo '==' $(date) '==' MT all candidate fastq compression ENDED
#

# align MT reads to ref for NUMT and foldback filtration + debugging
echo '==' $(date) '==' MT candidates fastq reference mapping STARTED
(${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -t ${MINIMAP2THREADS} \
"${MITOSCOPE_RESOURCES}/${MINIMAPINDEX}" \
${RESULTDIR}/${FASTQPREFIX}.MT.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam ; \
export DSNAME=${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} | ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg "${MITOSCOPE_RESOURCES}/MT.fasta.fai" ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' MT candidates fastq reference mapping COMPLETED
#

# for debugging
echo '==' $(date) '==' MT candidates fastq variation against reference STARTED
${SNIFFLESCMD} \
--output-rnames \
--qc-output-all --allow-overwrite \
--minsupport ${MINREADSUPPORT} \
--input ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam \
--vcf ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam.raw.vcf

${BCFTOOLSCMD} filter \
-i "SUPPORT>=${MINREADSUPPORT}" \
${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam.raw.vcf \
> ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam.raw.ge${MINREADSUPPORT}.vcf
echo '==' $(date) '==' MT candidates fastq variation against reference COMPLETED
#       

# remove foldback reads before assembly
echo '==' $(date) '==' Removal of foldback MT candidates STARTED
python ${MITOSCOPE_ROOT}/filter_bam.py -i ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam 

${SAMTOOLSCMD} fastq -@ ${SAMTOOLSTHREADS} ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.filtered.bam | pigz -p ${PIGZTHREADS} > ${RESULTDIR}/${FASTQPREFIX}.MT.filtered.fastq.gz
echo '==' $(date) '==' Removal of foldback MT candidates COMPLETED
#

# call indels using Mutect2
echo '==' $(date) '==' Mutect2 variant calling STARTED
mkdir -p ${RESULTDIR}/mutect

${GATKCMD} AddOrReplaceReadGroups \
-I ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.filtered.bam \
-O ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.filtered.RG.bam  \
--RGLB lib1 --RGPL wgs --RGPU unit1 --RGSM ${FASTQPREFIX} --SORT_ORDER coordinate --CREATE_INDEX true 

${GATKCMD} Mutect2 --mitochondria-mode \
-R ${MITOSCOPE_RESOURCES}/MT.fasta \
-I ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.filtered.RG.bam \
-O ${RESULTDIR}/mutect/${FASTQPREFIX}.MT.ref.filtered.mutect2_raw.vcf.gz \
--native-pair-hmm-threads ${MUTECTTHREADS}

${GATKCMD} FilterMutectCalls --mitochondria-mode \
-R ${MITOSCOPE_RESOURCES}/MT.fasta \
-V ${RESULTDIR}/mutect/${FASTQPREFIX}.MT.ref.filtered.mutect2_raw.vcf.gz \
-O ${RESULTDIR}/mutect/${FASTQPREFIX}.MT.ref.filtered.mutect2_filters.vcf.gz

${BCFTOOLSCMD} view -f PASS -Oz -o ${RESULTDIR}/mutect/${FASTQPREFIX}.MT.ref.filtered.mutect2_PASS.vcf.gz \
${RESULTDIR}/mutect/${FASTQPREFIX}.MT.ref.filtered.mutect2_filters.vcf.gz

echo '==' $(date) '==' Mutect2 variant calling COMPLETED
#

# call SNVs using mutserve 
echo '==' $(date) '==' Mutserve SNV calling STARTED
mkdir -p ${RESULTDIR}/mutserve

${MUTSERVECMD} call ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.filtered.bam \
--output ${RESULTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.vcf.gz \
--reference ${MITOSCOPE_RESOURCES}/MT.fasta \
--threads ${MUTSERVETHREADS} --no-ansi
echo '==' $(date) '==' Mutserve SNV calling COMPLETED
## rename txt output from mutserve since it always truncates full name
mv ${RESULTDIR}/mutserve/${FASTQPREFIX%%.*}.txt ${RESULTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.txt
#

# haplogroup classification using haplogrep3
echo '==' $(date) '==' Haplogroup classification STARTED
TREE="phylotree-rcrs@17.2"
${HAPLOGREPCMD} classify \
--tree ${TREE} \
--input ${RESULTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.vcf.gz \
--output ${RESULTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.haplogrep.txt \
--extend-report --write-qc
echo '==' $(date) '==' Haplogroup classification COMPLETED
#

# haplocheck contamination check
echo '==' $(date) '==' Haplocheck contamination check STARTED
${HAPLOCHECKCMD} --raw --out ${RESULTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.haplocheck.txt \
${RESULTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.vcf.gz 
echo '==' $(date) '==' Haplocheck contamination check COMPLETED
#

# add MITOMAP annotations to mutserve output 
## generate static anno file to use
echo '==' $(date) '==' MITOMAP annotation of mutserve output STARTED
Rscript ${MITOSCOPE_ROOT}/annotate.R \
${RESULTDIR}/mutserve/${FASTQPREFIX}.MT.ref.filtered.mutserve.txt \
${MITOSCOPE_ROOT}/annotations/CombinedDiseaseVariantDB.csv
echo '==' $(date) '==' MITOMAP annotation of mutserve output COMPLETED
#

# assemble the selected long-reads
echo '==' $(date) '==' de Novo MT assembly STARTED
${FLYECMD} --threads ${FLYETHREADS} --meta \
${FLYEPLATFORM} ${RESULTDIR}/${FASTQPREFIX}.MT.filtered.fastq.gz \
--out-dir ${RESULTDIR}/MT_assembly \
-m ${FLYEMINOVERLAP} &
wait
echo '==' $(date) '==' de Novo MT assembly COMPLETED
#

# map selected long-reads to assembled contig(s)
echo '==' $(date) '==' MT candidates fastq assembly mapping STARTED
(${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -t ${MINIMAP2THREADS} \
"${RESULTDIR}/MT_assembly/assembly.fasta" \
${RESULTDIR}/${FASTQPREFIX}.MT.filtered.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam ; \
export DSNAME=${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} | ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${SAMTOOLSCMD} faidx "${RESULTDIR}/MT_assembly/assembly.fasta" ; \
${BG2BWCMD} ${DSNAME}.bg "${RESULTDIR}/MT_assembly/assembly.fasta.fai" ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' MT candidates fastq assembly mapping COMPLETED
#

# call variations of selected long-reads w.r.t. assembled contig(s)
echo '==' $(date) '==' MT candidates fastq variation against assembly STARTED
${SNIFFLESCMD} \
--output-rnames \
--qc-output-all --allow-overwrite \
--minsupport ${MINREADSUPPORT} \
--input ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam \
--vcf ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf

${BCFTOOLSCMD} filter \
-i "SUPPORT>=${MINREADSUPPORT}" \
${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf \
> ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.ge${MINREADSUPPORT}.vcf
echo '==' $(date) '==' MT candidates fastq variation against assembly COMPLETED
#

# for inter-sample anchoring + debugging
echo '==' $(date) '==' Assembly reference mapping STARTED
(${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -t ${MINIMAP2THREADS} \
"${MITOSCOPE_RESOURCES}/${MINIMAPINDEX}" \
${RESULTDIR}/MT_assembly/assembly.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam ; \
export DSNAME=${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} | ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg "${MITOSCOPE_RESOURCES}/MT.fasta.fai" ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' Assembly reference mapping COMPLETED
#

# for inter-sample anchoring + debugging
echo '==' $(date) '==' Assembly variation against reference STARTED
${SNIFFLESCMD} \
--output-rnames \
--qc-output-all --allow-overwrite \
--input ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam \
--vcf ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam.raw.vcf
echo '==' $(date) '==' Assembly variation against reference COMPLETED
#


# for debugging
echo '==' $(date) '==' Graph_before reference mapping STARTED
(${MINIMAP2CMD} -ax ${MINIMAPPLATFORM} -t ${MINIMAP2THREADS} \
"${MITOSCOPE_RESOURCES}/${MINIMAPINDEX}" \
${RESULTDIR}/MT_assembly/20-repeat/graph_before_rr.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${DEBUGDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${DEBUGDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam ; \
export DSNAME=${DEBUGDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} | ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg "${MITOSCOPE_RESOURCES}/MT.fasta.fai" ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' Graph_before fastq reference mapping COMPLETED
#


# generate subpopulation
# TODO: determine at what threshold we claim a new subpopulation
echo '==' $(date) '==' Subpopulations generation STARTED
export SIEVEDGRAPHDIR=${FASTQDIR}/mitoscope/MT_assembly/sieved_graph
[ ! -d "${SIEVEDGRAPHDIR}" ] && mkdir -p "${SIEVEDGRAPHDIR}"
pushd "${SIEVEDGRAPHDIR}"
# 
echo '==' $(date) '==' check assembly graph for circular genome..
perl ${MITOSCOPE_ROOT}/ecLegov2.pl sievegraph \
--gv ../assembly_graph.gv \
--diploidcov $(expr 2 \* ${covCN}) \
--oprefix ${FASTQPREFIX}.MT.assembly \
--bam ../../${FASTQPREFIX}.MT.assembly.ref.bam
# FIXME: lower dicncov to 4 for SV threshold as >=2
echo '==' $(date) '==' subpopulation script setup..
perl ${MITOSCOPE_ROOT}/cgSupPop.pl setup \
--assembly ${FASTQPREFIX}.MT.assembly.cn$(expr 2 \* ${covCN}).disjointcyclic.gv.overview.xls \
--vcf ../../${FASTQPREFIX}.MT.assembly.bam.raw.vcf \
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


#
# TODO: <INV>, <DEL>, <INS> required different minimum read support
# TODO: SNV calling
#
# TODO: visualization : circos or CGView (https://js.cgview.ca/examples/index.html)
# NOTE: .gv visualization with https://dreampuf.github.io/GraphvizOnline/
#
# TODO: population based analysis and visualization?
#

