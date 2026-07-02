#!/bin/bash
#PBS -N Wget
#PBS -o Wget.log
#PBS -e Wget.err
#PBS -q wuyileiQueue
#PBS -l nodes=node09:ppn=8
#PBS -l mem=50GB

################################
### Variable information ###


ProjectDir=/data_group/wuyilei/Project/MulberryGenome/Syngap

cd ${ProjectDir}
source activate syngap_env

syngap dual \
    --sp1fa=../Morus_atropurpurea.Chrom.fa \
    --sp1gff=../Morus_atropurpurea.Chrom.gff \
    --sp2fa=../Morus_alba.Chrom.fa \
    --sp2gff=../Morus_alba.Chrom.gff \
    --sp1=atropurpurea \
    --sp2=alba

syngap genepair \
    --sp1fa=../Morus_atropurpurea.Chrom.fa \
    --sp1gff=atropurpurea.SynGAP.gff3 \
    --sp2fa=../Morus_alba.Chrom.fa \
    --sp2gff=alba.SynGAP.gff3 \
    --sp1=atropurpurea \
    --sp2=alba


syngap evi \
    --genepair=alba.atropurpurea.final.clean.genepair \
    --sp1exp=alba.S4_S6_S10.transcript.TPM.xls \
    --sp2exp=atrop.S4_S6_S10.transcript.TPM.xls