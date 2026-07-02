#!/bin/bash
#PBS -N minimap
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node09:ppn=24
#PBS -l mem=45GB


threads=24
source activate vep_env

cd /data_group/wuyilei/Project/MulberryGenome/Psvcp/Minimap
source activate vep_env
minimap2 -ax map-hifi -t ${threads} /data_group/wuyilei/Project/MulberryGenome/Morus_mongolica.Chrom.fa CRR778716.fq --MD > notabilis2mongolica.sam
samtools view -@ ${threads} -S -b -o notabilis2mongolica.bam notabilis2mongolica.sam
samtools sort -@ ${threads} notabilis2mongolica.bam -o notabilis2mongolica.sorted_hifi.bam
samtools index -@ ${threads} notabilis2mongolica.sorted_hifi.bam

source activate alphamissense_env
sniffles -m notabilis2mongolica.sorted_hifi.bam -v reads_sv.vcf