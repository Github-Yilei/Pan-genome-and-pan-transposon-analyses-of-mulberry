#!/bin/bash
#PBS -N EDTA
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node02:ppn=12
#PBS -l mem=40GB




################################
### Project information ###
cd /data_group/wuyilei/Project/ReAnnotGenome/01-augustsus
SampleList=/data_group/wuyilei/Project/ReAnnotGenome/samples.list

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

source activate vep_env
seqkit split -i ../${sampleID}.genome.fasta.new.masked -O ${sampleID}_Split

# parallel working
source activate lzf_py36_env
ls ${sampleID}_Split | while read name;
do
echo "augustus --species=arabidopsis --gff3=on --softmasking=1 --UTR=on ${sampleID}_Split/${name} > ${name}.gff" >> ${sampleID}.cmd.list;
done

/opt/service/miniconda3/envs/r4_py37_env/bin/ParaFly -c ${sampleID}.cmd.list -CPU 12

ls ${sampleID}_Split | while read name;
do
cat ${name}.gff | perl /opt/service/miniconda3/envs/lzf_py36_env/bin/join_aug_pred.pl | grep -v '^#' >> ${sampleID}.joined.gff;
done


source activate vep_env
bedtools sort -i ${sampleID}.joined.gff > ${sampleID}.gff3
rm *joined.gff
rm -rf *_Split
