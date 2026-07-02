#!/bin/bash
#PBS -N GeneMark
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node03:ppn=24
#PBS -l mem=40GB




################################
### Project information ###
cd /data_group/wuyilei/Project/ReAnnotGenome/02-Gmes
SampleList=/data_group/wuyilei/Project/ReAnnotGenome/samples.list

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`


threads=24
source activate GALBA_env

mkdir ${sampleID}
/opt/service/miniconda3/envs/GALBA_env/GeneMark-ETP/bin/gmes/gmes_petap.pl --ES \
        --sequence ../${sampleID}.genome.fasta.new.masked \
        --cores ${threads} --work_dir ${sampleID}
