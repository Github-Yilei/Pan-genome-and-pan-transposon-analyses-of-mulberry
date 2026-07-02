#!/bin/bash
#PBS -N FastpQC
#PBS -o QsubLog/
#PBS -e QsubLog/
#PBS -q cypQueue
#PBS -l nodes=node23:ppn=12
#PBS -l mem=40GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/cunyupeng/huling/daiyi/Project/Mulberry


#mkdir -p ${RawFastQC}
#mkdir -p ${CleanReads}
#mkdir -p ${CleanFastQC}
#mkdir -p ${FastpDir}
#mkdir -p ${FastpLog}
#mkdir -p ${PBSLog}


################################
### Tools ###
Fastp=fastp
Fastqc=fastqc

################################
### Main parameters ###
RequiredCPU=5

################################
### Main process ###
cd ${ProjectDir}
source activate pairadise_env

STAR --runMode genomeGenerate \
     --runThreadN 12 \
     --genomeDir Index_Star \
     --genomeFastaFiles Morusnotabilis.genome.fasta \
     --sjdbGTFfile Morusnotabilis.gff \
     --sjdbOverhang 100 --sjdbGTFfeatureExon CDS
