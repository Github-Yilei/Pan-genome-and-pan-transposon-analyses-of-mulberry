#!/bin/bash
#PBS -N trinity
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node03:ppn=20
#PBS -l mem=220GB




################################
### Project information ###
cd /data_group/wuyilei/Project/ReAnnotGenome/05-PASA
SampleList=/data_group/wuyilei/Project/ReAnnotGenome/RNA.list

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

threads=20
source activate dyw_py36_env
mkdir ${sampleID} && cd ${sampleID}
Trinity \
        --seqType fq --CPU ${threads}  --max_memory 220G \
        --samples_file /data_group/wuyilei/Project/ReAnnotGenome/05-PASA/sample.tsv \
        --output Trinity_out


#basic PASA 
Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g ../${sampleID}.genome.fasta.new.masked -t Trinity.fasta --ALIGNERS blat