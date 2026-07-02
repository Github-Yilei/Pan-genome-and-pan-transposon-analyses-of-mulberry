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
minimap2 -ax map-hifi -t ${threads} /data_group/wuyilei/Project/MulberryGenome/Morus_notabilis.Chrom.fa CRR778716.fq --MD > notabilis2notabilis.sam
samtools view -@ ${threads} -S -b -o notabilis2notabilis.bam notabilis2notabilis.sam
samtools sort -@ ${threads} notabilis2notabilis.bam -o notabilis2notabilis.sorted_hifi.bam
samtools index -@ ${threads} notabilis2notabilis.sorted_hifi.bam
sambamba depth window -t ${threads} -w 1000 notabilis2notabilis.sorted_hifi.bam > notabilis2notabilis.sorted_hifi_1000_depth.txt

minimap2 -ax map-hifi -t ${threads} /data_group/wuyilei/Project/MulberryGenome/Morus_notabilis.Chrom.fa SRR31188304.fastq.gz  > mongolica2notabilis.sam
samtools view -@ ${threads} -S -b -o mongolica2notabilis.bam mongolica2notabilis.sam
samtools sort -@ ${threads} mongolica2notabilis.bam -o mongolica2notabilis.sorted_hifi.bam
samtools index -@ ${threads} mongolica2notabilis.sorted_hifi.bam
sambamba depth window -t ${threads} -w 1000 mongolica2notabilis.sorted_hifi.bam > mongolica2notabilis.sorted_hifi_1000_depth.txt

