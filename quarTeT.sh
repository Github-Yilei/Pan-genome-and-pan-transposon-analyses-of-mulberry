#!/bin/bash
#PBS -N CentroMiner
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node07:ppn=5
#PBS -l mem=20GB
cd /data_group/wuyilei/Project/MulberryGenome/Psvcp
source activate syri_env

python3 /opt/service/miniconda3/envs/syri_env/quarTeT/quartet.py CentroMiner -t 20 -i Morus_mongolica.Chrom.fa --gene Morus_mongolica.Chrom.gff  \
    --TE /data_group/wuyilei/Project/PanTE/05_ReEDTA/Morus_mongolica.genome.fasta.mod.EDTA.TEanno.gff3 -p mongolica_quartet

python3 /opt/service/miniconda3/envs/syri_env/quarTeT/quartet.py TeloExplorer -i Morus_mongolica.Chrom.fa -c plant -p mongolica_TeloExplorer
