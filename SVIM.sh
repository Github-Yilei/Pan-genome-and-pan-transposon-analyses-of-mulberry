#!/bin/bash
#PBS -N callSV
#PBS -o callSV.log
#PBS -e callSV.err
#PBS -q wuyileiQueue
#PBS -l nodes=node03:ppn=24
#PBS -l mem=40GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/MulberryGenome/SVIM

SampleList=/data_group/wuyilei/Project/MulberryGenome/SVIM/samples.list

################################
### Tools ###

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

### Main process ###
cd ${ProjectDir}
source activate vep_env
minimap2 -a -x asm5 --cs -r2k -t 24  Morus_mongolica.Chrom.fa  ../${sampleID}.Chrom.fa > ${sampleID}.sam
samtools sort -m 40G -@ 24 -o ${sampleID}.sorted.bam ${sampleID}.sam
samtools index ${sampleID}.sorted.bam


source activate zpy_assembly_env
svim-asm haploid --min_sv_size 40 --query_names --sample ${sampleID}  ${sampleID} ${sampleID}.sorted.bam PanGenome.fa


source activate pms_py27_env
ls output.vcfs/*.vcf > vcf.list

SURVIVOR merge vcf.list 0.1 0 1 0 0 100 svim_asm_variants.vcf


awk '!/^#/{print $1"\t"$2"\t"$3"\t"$8}' svim_asm_variants.vcf | \
    sed 's#SUPP.*;SVLEN=##' | sed 's#;SVTYPE.*##' | awk '{if($4>0) {print $1"\t"$2"\t"$2+$4"\t"$3} else {print $1"\t"$2+$4"\t"$2"\t"$3}}'> svim_asm_variants.bed