#!/bin/bash
#PBS -N EDTA
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node04:ppn=8
#PBS -l mem=80GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/PanTE/05_ReEDTA

SampleList=/data_group/wuyilei/Project/PanTE/samples.list

### Tools ###

################################
### Main parameters ###
RequiredCPU=8

IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

################################
### Main process ###
source activate edta_env

cd ${ProjectDir}


EDTA.pl --genome ${ProjectDir}/Genome/PanGenome.fa -t ${RequiredCPU} --anno 1 --curatedlib panEDTA.TElib.FixedId.fa
