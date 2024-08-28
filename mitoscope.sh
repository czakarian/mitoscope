#!/bin/bash
######################################################################
# Copyright (C) 2019-2024 Dr Chia-Lin Wei Laboratory, NWGC
# Adapted from ecLegoV3 for human mitochondria
#
# PENDING:
# C1. SNV/SNP calling support
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
# T2. containerize sniffles
# T3. use containerized samtools in cgSupPop.pl, cgbam.pl, ecLegov2.pl
#
# NOTE:
# N1. .gv visualization with https://dreampuf.github.io/GraphvizOnline/
#
#
######################################################################


set -eu -o pipefail

export FASTQ=${1}
export MINREADSUPPORT=${2:-"2"}
export CNCOV=${3:-""}
[ -z "${FASTQ}" ] && { echo "Fastq file is empty" ; exit 1; }
[ ! -f "${FASTQ}" ] && { echo "${FASTQ} does not exist" ; exit 2 ; }

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
export KMCTOOLSCMD="${MITOSCOPE_SINGULARITY}/kmc_3.2.1.sif kmc_tools"
export FLYECMD="${MITOSCOPE_SINGULARITY}/flye_2.9.1.sif flye"
export MINIMAP2CMD="${MITOSCOPE_SINGULARITY}/minimap2_2.24.sif minimap2"
export SAMTOOLSCMD="${MITOSCOPE_SINGULARITY}/samtools_v1.15.1.sif samtools"
export GENOMECOVERAGEBEDCMD="${MITOSCOPE_SINGULARITY}/bedtools_2.31.0.sif genomeCoverageBed"
export SORTBEDCMD="${MITOSCOPE_SINGULARITY}/bedtools_2.31.0.sif sortBed"
export BG2BWCMD="${MITOSCOPE_SINGULARITY}/ucsc-bedgraphtobigwig_445.sif bedGraphToBigWig"

export FASTQDIR="$(dirname "$(readlink -f "${FASTQ}")")"
export RESULTDIR=${FASTQDIR}/mitoscope
[ ! -d "${RESULTDIR}" ] && mkdir "${RESULTDIR}"
export DEBUGDIR=${FASTQDIR}/mitoscope/debug
[ ! -d "${DEBUGDIR}" ] && mkdir "${DEBUGDIR}"

export FASTQPREFIX="$(basename "$(readlink -f "${FASTQ}")")"
export FASTQPREFIX="${FASTQPREFIX%.gz}"
export FASTQPREFIX="${FASTQPREFIX%.*}"

export SINGULARITY_BINDPATH="${MITOSCOPE_ROOT},${FASTQDIR}"

export KMCTOOLSTHREADS=4
export PIGZTHREADS=4
export MINIMAP2THREADS=4
export SAMTOOLSTHREADS=4
export FLYETHREADS=4
export FLYEMINOVERLAP=2500

module load modules{,-init,-gs}
module load samtools/1.17
module load bcftools/1.19
# TODO: to containerize
module load python/3.10.14_spec
module load sniffles/2.3.3

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

# assemble the selected long-reads
echo '==' $(date) '==' de Novo MT assembly STARTED
${FLYECMD} --threads ${FLYETHREADS} --meta \
--nano-hq ${RESULTDIR}/${FASTQPREFIX}.MT.fastq.gz \
--out-dir ${RESULTDIR}/MT_assembly \
-m ${FLYEMINOVERLAP} &
wait
echo '==' $(date) '==' de Novo MT assembly COMPLETED
#

# map selected long-reads to assembled contig(s)
echo '==' $(date) '==' MT candidates fastq assembly mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${RESULTDIR}/MT_assembly/assembly.fasta" \
${RESULTDIR}/${FASTQPREFIX}.MT.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam ; \
export DSNAME=${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} \
| ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${SAMTOOLSCMD} faidx "${RESULTDIR}/MT_assembly/assembly.fasta" ; \
${BG2BWCMD} ${DSNAME}.bg "${RESULTDIR}/MT_assembly/assembly.fasta.fai" ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' MT candidates fastq assembly mapping COMPLETED
#

# TODO: to containerize
# call variations of selected long-reads w.r.t. assembled contig(s)
echo '==' $(date) '==' MT candidates fastq variantion against assembly STARTED
sniffles \
--output-rnames \
--qc-output-all --allow-overwrite \
--minsupport ${MINREADSUPPORT} \
--input ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam \
--vcf ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf
bcftools filter -i "SUPPORT>=${MINREADSUPPORT}" ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.vcf \
> ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.bam.raw.ge${MINREADSUPPORT}.vcf
echo '==' $(date) '==' MT candidates fastq variantion against assembly COMPLETED
#

# for inter-sample anchoring + debugging
echo '==' $(date) '==' Assembly reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${MITOSCOPE_RESOURCES}/MT.mmi" \
${RESULTDIR}/MT_assembly/assembly.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam ; \
export DSNAME=${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} \
| ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg "${MITOSCOPE_RESOURCES}/MT.fasta.fai" ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' Assembly reference mapping COMPLETED
#

# TODO: to containerize
# for inter-sample anchoring + debugging
echo '==' $(date) '==' Assembly variation against reference STARTED
sniffles \
--output-rnames \
--qc-output-all --allow-overwrite \
--input ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam \
--vcf ${RESULTDIR}/${FASTQPREFIX}.MT.assembly.ref.bam.raw.vcf
echo '==' $(date) '==' Assembly variation against reference COMPLETED
#

# for debugging
echo '==' $(date) '==' MT candidates fastq reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${MITOSCOPE_RESOURCES}/MT.mmi" \
${RESULTDIR}/${FASTQPREFIX}.MT.fastq.gz \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam ; \
export DSNAME=${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} \
| ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
echo '==' $(date) '==' Generating bigwig for ${DSNAME}.. ; \
${BG2BWCMD} ${DSNAME}.bg "${MITOSCOPE_RESOURCES}/MT.fasta.fai" ${DSNAME}.bw && rm ${DSNAME}.bg ; \
echo '==' $(date) '==' Done. ) &
wait
echo '==' $(date) '==' MT candidates fastq reference mapping COMPLETED
#

# TODO: to containerize
# for debugging
echo '==' $(date) '==' MT candidates fastq variantion against reference STARTED
sniffles \
--output-rnames \
--qc-output-all --allow-overwrite \
--minsupport ${MINREADSUPPORT} \
--input ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam \
--vcf ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam.raw.vcf
bcftools filter -i "SUPPORT>=${MINREADSUPPORT}" ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam.raw.vcf \
> ${DEBUGDIR}/${FASTQPREFIX}.MT.ref.bam.raw.ge${MINREADSUPPORT}.vcf
echo '==' $(date) '==' MT candidates fastq variantion against reference COMPLETED
#

# for debugging
echo '==' $(date) '==' Graph_before reference mapping STARTED
(${MINIMAP2CMD} -ax map-ont -t ${MINIMAP2THREADS} \
"${MITOSCOPE_RESOURCES}/MT.mmi" \
${RESULTDIR}/MT_assembly/20-repeat/graph_before_rr.fasta \
| ${SAMTOOLSCMD} sort -O BAM -@${SAMTOOLSTHREADS} -o ${DEBUGDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam ; \
${SAMTOOLSCMD} index -@${SAMTOOLSTHREADS} ${DEBUGDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam ; \
export DSNAME=${DEBUGDIR}/${FASTQPREFIX}.MT.graph_before_rr.ref.bam; \
export CNSCALE=$(awk -v DICNSCALE=$(expr 2 \* ${covCN}) 'BEGIN{print 2/DICNSCALE}') ; \
echo '==' $(date) '==' Generating bedgraph for ${DSNAME}.. ; \
${GENOMECOVERAGEBEDCMD} -bg -split -scale ${CNSCALE} -ibam ${DSNAME} \
| ${SORTBEDCMD} -i - > ${DSNAME}.bg ; \
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
echo '==' $(date) '==' variations calling on possible subpopulations..
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

