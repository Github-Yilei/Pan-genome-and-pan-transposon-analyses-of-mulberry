#!/bin/bash
#PBS -N gatk
#PBS -o QsubLog2/
#PBS -e QsubLog2/
#PBS -q wuyileiQueue
#PBS -l nodes=node05:ppn=2
#PBS -l mem=10GB



################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/MulberryFruit

SampleList=${ProjectDir}/sample2RNAseq.tsv



IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`
################################
### Main process ###
cd ${ProjectDir}
source activate vep_env

stringtie  -e -B -p 4 -G Morusnotabilis.gff Star/${sampleID}Aligned.sortedByCoord.out.rg_added.bam -o \
    Result/${sampleID}/${sampleID}.transcripts.gtf -A Result/${sampleID}/${sampleID}_gene_abundances.tsv
