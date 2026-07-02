#!/bin/bash
#PBS -N Orthofinder
#PBS -o Orthofinder.log
#PBS -e Orthofinder.err
#PBS -q wuyileiQueue
#PBS -l nodes=node01:ppn=10
#PBS -l mem=20GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/MulberryGenome/JCVI
SampleList=/data_group/wuyilei/Project/MulberryGenome/SampleList.tsv


## prep data
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

python -m jcvi.formats.gff bed --type=mRNA --primary_only ${sampleID}.Chrom.gff > ${sampleID}.Chrom.bed
ln -s /data_group/wuyilei/Project/MulberryGenome/${sampleID}.Protein.longest.fa ./


## run JCVI
python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names Artocarpus_heterophyllus Morus_alba
python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names Artocarpus_heterophyllus Morus_atropurpurea
python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names Artocarpus_heterophyllus Morus_mongolica
python -m jcvi.compara.catalog ortholog --dbtype prot --no_strip_names Artocarpus_heterophyllus Morus_yunnanensis