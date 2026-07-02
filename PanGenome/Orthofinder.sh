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
ProjectDir=/data_group/wuyilei/Project/MulberryGenome

################################
### Main process ###
cd ${ProjectDir}


source activate PanSyn_env

orthofinder -f ./ -M msa -S diamond -T fasttree -t 20 -a 20
