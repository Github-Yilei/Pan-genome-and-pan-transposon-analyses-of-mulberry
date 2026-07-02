#!/bin/bash
#PBS -N psvcp
#PBS -o psvcp.log
#PBS -e psvcp.err
#PBS -q wuyileiQueue
#PBS -l nodes=node07:ppn=40
#PBS -l mem=200GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/MulberryGenome/Psvcp


################################
### Tools ###
#Fastp=fastp
#Fastqc=fastqc

################################
### Main parameters ###

#IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
#sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`


################################
### Main process ###
cd ${ProjectDir}
source activate syri_env

python3 /opt/service/miniconda3/envs/syri_env/psvcp_v1.01/1Genome_construct_Pangenome.py -t 40 ./ genome_list
