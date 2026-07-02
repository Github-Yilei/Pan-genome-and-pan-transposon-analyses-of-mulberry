#!/bin/bash
#PBS -N ParaAT
#PBS -o ParaAT.log
#PBS -e ParaAT.err
#PBS -q wuyileiQueue
#PBS -l nodes=node07:ppn=22
#PBS -l mem=180GB


ProjectDir=/data_group/wuyilei/Project/MulberryGenome/DupGen_finder/Heterophyllus
CDS_path=/data_group/wuyilei/Project/MulberryGenome/DupGen_finder/Artocarpus_heterophyllus.CDS.longest.fa
pep_path=/data_group/wuyilei/Project/MulberryGenome/DupGen_finder/Artocarpus_heterophyllus.Protein.longest.fa

source activate vep_env

# run info
export PATH=/opt/service/ParaAT2.0/Epal2nal.pl:$PATH
export PATH=/opt/service/ParaAT2.0:$PATH
export PATH=/opt/service/KaKs_Calculator2.0/src/AXTConvertor:$PATH


cd ${ProjectDir}
# dispersed
awk 'NR>1{print $1"\t"$3}' results/*.dispersed.pairs > tmp

/opt/service/ParaAT2.0/ParaAT.pl  -h tmp \
    -n ${CDS_path} \
    -a ${pep_path} -p proc \
    -m mafft -f axt -g -k -o Dispersed_dir

# proximal
awk 'NR>1{print $1"\t"$3}' results/*.proximal.pairs > tmp

/opt/service/ParaAT2.0/ParaAT.pl  -h tmp \
    -n ${CDS_path} \
    -a ${pep_path} -p proc \
    -m mafft -f axt -g -k -o Proximal_dir

# tandem
awk 'NR>1{print $1"\t"$3}' results/*.tandem.pairs > tmp

/opt/service/ParaAT2.0/ParaAT.pl  -h tmp \
    -n ${CDS_path} \
    -a ${pep_path} -p proc \
    -m mafft -f axt -g -k -o Tandem_dir

# transposed
awk 'NR>1{print $1"\t"$3}' results/*.transposed.pairs > tmp

/opt/service/ParaAT2.0/ParaAT.pl  -h tmp \
    -n ${CDS_path} \
    -a ${pep_path} -p proc \
    -m mafft -f axt -g -k -o Transposed_dir

# wgd
awk 'NR>1{print $1"\t"$3}' results/*.wgd.pairs > tmp

/opt/service/ParaAT2.0/ParaAT.pl  -h tmp \
    -n ${CDS_path} \
    -a ${pep_path} -p proc \
    -m mafft -f axt -g -k -o Wgd_dir
