#!/bin/bash
#PBS -N CPC2
#PBS -o CPC2.log
#PBS -e CPC2.err
#PBS -q wuyileiQueue
#PBS -l nodes=node08:ppn=5
#PBS -l mem=100GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/MulberryGenome/CPC2

cd ${ProjectDir}


source activate tadlib_env
sample_id=Morus_yunnanensis
CPC2.py -i ../${sample_id}.mRNA.fa -o ${sample_id}.CPC2
svm-scale -r /opt/service/miniconda3/envs/tadlib_env/data/cpc2.range ${sample_id}.CPC2.tmp.1 > ${sample_id}.CPC2.tmp.2
svm-predict -b 1 -q ${sample_id}.CPC2.tmp.2 /opt/service/miniconda3/envs/tadlib_env/data/cpc2.model ${sample_id}.CPC2.tmp.out
CPC2.py -i ../${sample_id}.mRNA.fa -o ${sample_id}.CPC2
rm *feat *tmp*


grep -f Unpaired_alba_gene.idx Morus_alba.CPC2.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tUnpaired"}' > Unpaired_alba
grep -f paired_alba_gene.idx Morus_alba.CPC2.txt |  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tPaired"}' > paired_alba

grep -f Unpaired_atropurpurea_gene.idx Morus_atropurpurea.CPC2.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tUnpaired"}' > Unpaired_atrop
grep -f paired_atropurpurea_gene.idx Morus_atropurpurea.CPC2.txt | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tPaired"}' > paired_atrop
