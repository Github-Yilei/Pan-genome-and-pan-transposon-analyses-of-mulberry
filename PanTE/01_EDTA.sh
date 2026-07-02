#!/bin/bash
#PBS -N EDTA
#PBS -o QsubLog/EDTA/
#PBS -e QsubLog/EDTA/
#PBS -q wuyileiQueue
#PBS -l nodes=node09:ppn=24
#PBS -l mem=100GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/PanTE

SampleList=/data_group/wuyilei/Project/PanTE/samples.list




################################
### Tools ###

################################
### Main parameters ###
RequiredCPU=24

IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Main process ###
source activate edta_env


cd ${ProjectDir}/01_EDTA


EDTA.pl -genome ${ProjectDir}/Genome/${sampleID}.genome.fasta -species others -step all -t ${RequiredCPU} -sensitive 1 -anno 1 -evaluate 1
