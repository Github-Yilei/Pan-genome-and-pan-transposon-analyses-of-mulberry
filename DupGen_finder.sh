#!/bin/bash
#PBS -N blast
#PBS -o blast.log
#PBS -e blast.err
#PBS -q wuyileiQueue
#PBS -l nodes=node08:ppn=3
#PBS -l mem=50GB

################################
### Variable information ###

################################
### Project information ###
ProjectDir=/data_group/wuyilei/Project/MulberryGenome/DupGen_finder


################################
### Tools ###
source activate vep_env
################################
### Main parameters ###
RequiredCPU=3



################################
### Main process ###
cd ${ProjectDir}
blastp -num_threads ${RequiredCPU} -query Morus_yunnanensis.Protein.longest.fa \
    -db Morus_yunnanensis.db -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -out ${ProjectDir}/Yunnanensis/yunnanensis.blast

blastp -num_threads ${RequiredCPU} -query Morus_yunnanensis.Protein.longest.fa \
    -db Csativa.db  -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -out ${ProjectDir}/Yunnanensis/yunnanensis_Csativa.blast




blastp -num_threads ${RequiredCPU} -query Artocarpus_heterophyllus.Protein.longest.fa \
    -db Artocarpus_heterophyllus.db -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -out ${ProjectDir}/Heterophyllus/heterophyllus.blast

blastp -num_threads ${RequiredCPU} -query Artocarpus_heterophyllus.Protein.longest.fa \
    -db Csativa.db  -evalue 1e-10 -max_target_seqs 5 -outfmt 6 -out ${ProjectDir}/Heterophyllus/heterophyllus_Csativa.blast

