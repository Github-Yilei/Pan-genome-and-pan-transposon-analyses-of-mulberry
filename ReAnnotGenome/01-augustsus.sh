#!/bin/bash
#PBS -N EDTA
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node03:ppn=24
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
seqkit split -i ${sampleID}.genome.fasta.new.masked -O Split

# parallel working
source activate lzf_py36_env
ls Split | while read name;
do
echo "augustus --species=arabidopsis --gff3=on --softmasking=1 --UTR=on Split/${name} > ${name}.gff " >> cmd.list;
done

ParaFly -c cmd.list -CPU 24

ls Split | while read name;
do
cat ${name}.gff | perl /opt/service/miniconda3/envs/lzf_py36_env/bin/join_aug_pred.pl | grep -v '^#' >> temp.joined.gff;
done


source activate vep_env
bedtools sort -i temp.joined.gff > ${sampleID}.gff3

rm -rf Split
