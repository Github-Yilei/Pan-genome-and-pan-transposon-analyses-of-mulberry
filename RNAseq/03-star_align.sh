#!/bin/bash
#PBS -N FastpQC
#PBS -o QsubLog/
#PBS -e QsubLog/
#PBS -q wuyileiQueue
#PBS -l nodes=node03:ppn=12
#PBS -l mem=20GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/MulberryFruit

SampleList=${ProjectDir}/bach2rna.idx

CleanReads=${ProjectDir}/CleanReads



#mkdir -p ${RawFastQC}
#mkdir -p ${CleanReads}
#mkdir -p ${CleanFastQC}
#mkdir -p ${FastpDir}
#mkdir -p ${FastpLog}
#mkdir -p ${PBSLog}


################################
### Tools ###

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Main process ###
cd ${ProjectDir}
source activate pairadise_env

STAR \
    --genomeDir Index_Star \
     --runThreadN 12 \
    --readFilesIn  ${CleanReads}/${sampleID}_1.fq.gz ${CleanReads}/${sampleID}_2.fq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix Star/${sampleID} \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 12 \
        --outReadsUnmapped Fastx
