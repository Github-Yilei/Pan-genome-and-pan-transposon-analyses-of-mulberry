#!/bin/bash
#PBS -N FastpQC
#PBS -o QsubLog/FastpLog/
#PBS -e QsubLog/FastpLog/
#PBS -q cypQueue
#PBS -l nodes=node24:ppn=5
#PBS -l mem=20GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/cunyupeng/huling/daiyi/Project/Mulberry

SampleList=${ProjectDir}/samples.tsv

RawReads=/data_group/cunyupeng/huling/daiyi/Project/Mulberry/NewReads

RawFastQC=${ProjectDir}/FastQCDir/RawFastQC
CleanFastQC=${ProjectDir}/FastQCDir/CleanFastQC
CleanReads=${ProjectDir}/CleanReads

FastpDir=${ProjectDir}/FastpDir
FastpLog=${ProjectDir}/FastpDir/FastpLog

PBSLog=${ProjectDir}/QsubLog/FastpLog

mkdir -p ${RawFastQC}
mkdir -p ${CleanReads}
mkdir -p ${CleanFastQC}
mkdir -p ${FastpDir}
mkdir -p ${FastpLog}
mkdir -p ${PBSLog}


################################
### Tools ###
Fastp=fastp
Fastqc=fastqc

################################
### Main parameters ###
RequiredCPU=5

IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Main process ###
cd ${ProjectDir}

### FastQC for raw reads ###

### Quality control ###
${Fastp} -i ${RawReads}/${sampleID}_1.fastq.gz -I ${RawReads}/${sampleID}_2.fastq.gz -o ${CleanReads}/${sampleID}_1.fq.gz \
        -O ${CleanReads}/${sampleID}_2.fq.gz -W 5 -M 20 -5 -3 -l 50 -w ${RequiredCPU} \
        -j ${FastpDir}/${sampleID}.json -h ${FastpDir}/${sampleID}.html > ${FastpLog}/${sampleID}.qc.log 2>&1
