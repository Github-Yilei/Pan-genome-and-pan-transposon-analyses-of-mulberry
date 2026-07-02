#!/bin/bash
#PBS -N EDTA
#PBS -o EDTA.log
#PBS -e EDTA.err
#PBS -q wuyileiQueue
#PBS -l nodes=node09:ppn=8
#PBS -l mem=100GB


#### main parameters####
workdir=/data_group/wuyilei/Project/PanTE/04_DeepTE

########################
source activate edta_env

#wget https://de.cyverse.org/dl/d/89D2FE7A-41BA-4F64-80E2-B9C26D49E99F/Plants_model.tar.gz


deepte=/opt/service/miniconda3/envs/edta_env/DeepTE/DeepTE.py
modeldir=${workdir}/Plants_model
cd ${workdir}


## reannot LTR
#python $deepte -i LTR_unknown.fa -sp P -m_dir $modeldir -fam LTR
python3 /opt/service/miniconda3/envs/edta_env/DeepTE/DeepTE.py -i LTR_unknown.fa -sp P -o output_dir -d ./ -m P

sed 's/LTR\/unknown__ClassI_LTR_Copia/LTR\/Copia/' opt_DeepTE.fasta \
    | sed 's/LTR\/unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/' \
    | sed 's/LTR\/unknown__ClassI_LTR/LTR\/unknown/' \
    > LTR_unknown_DeepTE.fa

## reannot TEs
python $deepte -i TE_unknown.fa -sp P -m_dir $modeldir
sed 's/Unknown__ClassI_LTR_Copia/LTR\/Copia/' opt_DeepTE.fasta \
    | sed 's/Unknown__ClassI_LTR_Gypsy/LTR\/Gypsy/' \
    | sed 's/Unknown__ClassI_LTR/LTR\/unknown/' \
    | sed 's/Unknown__ClassI_nLTR_DIRS/DIRS/' \
    | sed 's/Unknown__ClassIII_Helitron/DNA\/Helitron/' \
    | sed 's/Unknown__unknown/Unknown/' \
    | sed 's/Unknown__ClassII_nMITE/DNA\/unknown/' \
    | sed 's/Unknown__ClassII_MITE/MITE\/unknown/' \
    | sed 's/Unknown__ClassI_nLTR_LINE_I/LINE\/RII/' \
    | sed 's/Unknown__ClassI_nLTR_LINE_L1/LINE\/RIL/' \
    | sed 's/Unknown__ClassI_nLTR_LINE/LINE\/unknown/' \
    | sed 's/Unknown__ClassI_nLTR_PLE/nonLTR\/RPP/' \
    | sed 's/Unknown__ClassI_nLTR_SINE_tRNA/SINE\/RST/' \
    | sed 's/Unknown__ClassI_nLTR_SINE/SINE\/unknown/' \
    | sed 's/Unknown__ClassI_nLTR/nonLTR\/unknown/' \
    | sed 's/Unknown__ClassII_DNA_CACTA_MITE/MITE\/DTC/' \
    | sed 's/Unknown__ClassII_DNA_CACTA_nMITE/DNA\/DTC/' \
    | sed 's/Unknown__ClassII_DNA_CACTA_unknown/Unknown/' \
    | sed 's/Unknown__ClassII_DNA_Harbinger_MITE/MITE\/DTH/' \
    | sed 's/Unknown__ClassII_DNA_Harbinger_nMITE/DNA\/DTH/' \
    | sed 's/Unknown__ClassII_DNA_Harbinger_unknown/Unknown/' \
    | sed 's/Unknown__ClassII_DNA_hAT_MITE/MITE\/DTA/' \
    | sed 's/Unknown__ClassII_DNA_hAT_nMITE/DNA\/DTA/' \
    | sed 's/Unknown__ClassII_DNA_hAT_unknown/Unknown/' \
    | sed 's/Unknown__ClassII_DNA_Mutator_MITE/MITE\/DTM/' \
    | sed 's/Unknown__ClassII_DNA_Mutator_nMITE/DNA\/DTM/'\
    | sed 's/Unknown__ClassII_DNA_Mutator_unknown/Unknown/' \
    | sed 's/Unknown__ClassII_DNA_P_nMITE/DNA\/DTP/' \
    | sed 's/Unknown__ClassII_DNA_TcMar_MITE/MITE\/DTT/' \
    | sed 's/Unknown__ClassII_DNA_TcMar_nMITE/DNA\/DTT/' \
    | sed 's/Unknown__ClassII_DNA_TcMar_unknown/Unknown/' \
    | sed 's/Unknown__ClassII/Unknown/' \
    | sed 's/Unknown__ClassI/Unknown/' \
    > TE_unknown_DeepTE.fa



cat LTR_unknown_DeepTE.fa TE_unknown_DeepTE.fa TE_known.fa  > TElib_All_DeepTE.fa
python DeepTE2Fixid.py > panEDTA.TElib.FixedId.fa

echo "`date` All done!"

