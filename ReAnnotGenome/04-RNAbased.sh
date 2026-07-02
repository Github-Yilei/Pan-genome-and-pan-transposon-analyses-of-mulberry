#!/bin/bash
#PBS -N RNAbased
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node08:ppn=8
#PBS -l mem=45GB




################################
### Project information ###
cd /data_group/wuyilei/Project/ReAnnotGenome/04-RNAbased
SampleList=/data_group/wuyilei/Project/ReAnnotGenome/samples.list

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
RNAID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

source activate vep_env
mkdir ${sampleID}
hisat2-build ../${sampleID}.genome.fasta.new.masked  ${sampleID}
hisat2 --dta -p ${threads} -x ${sampleID}/${sampleID} \
       -1 /data_group/wuyilei/Project/MulberryFruit/CleanReads/${RNAID}_1.fq.gz \
       -2 /data_group/wuyilei/Project/MulberryFruit/CleanReads/${RNAID}_2.fq.gz -S ${sampleID}/${RNAID}.sam

samtools sort -@ ${threads} ${sampleID}/${RNAID}.sam -O BAM -o ${sampleID}/${RNAID}.bam

samtools merge -@ ${threads}  *.bam 

source activate tadlib_env
stringtie -p 8 -o ${sampleID}.gtf ${sampleID}.bam
