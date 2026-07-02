#!/bin/bash
#PBS -N EDTA
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node02:ppn=24
#PBS -l mem=100GB



#### main parameters####
path=/opt/service/miniconda3/envs/edta_env/share/EDTA/bin/
dir=/data_group/wuyilei/Project/PanTE/02_PanTE
genome_list=${dir}/samples.list
threads=10


#######################



source activate edta_env
# Set paths
cd ${dir}

# aggregate TE libs
i=0
for j in `cat $genome_list`; do
        genome=`basename $(echo $j|awk '{print $1}') 2>/dev/null`

        # skip empty lines
        if [ $genome == '' ]; then
                break
        fi

        i=$(($i+5000));
        perl $path/rename_TE.pl $genome.mod.EDTA.TElib.fa.keep.ori $i;
done | perl $path/rename_TE.pl - > panEDTA.TElib.fa.raw

# remove redundant
echo "Generate the panEDTA library"
perl $path/cleanup_nested.pl -in panEDTA.TElib.fa.raw -cov 0.80 -minlen 80 -miniden 80 -t $threads
cp panEDTA.TElib.fa.raw.cln panEDTA.TElib.fa

# Extra step if --curatedlib is provided:
if [ -s "$curatedlib" ]; then
        cat $curatedlib >> panEDTA.TElib.fa
fi
echo `date`
echo -e "\tpanEDTA library of is generated!"

### prep for DeepTE
source activate detectEVE_env
grep "LTR/unknown" panEDTA.TElib.fa | sed 's/>//' > tmp.idx
seqkit grep -f tmp.idx panEDTA.TElib.fa -o LTR_unknown.fa

grep "\#unknown" panEDTA.TElib.fa | grep -v "LTR/unknown" | sed 's/>//' > tmp.idx
seqkit grep -f tmp.idx panEDTA.TElib.fa -o TE_unknown.fa

grep '>' panEDTA.TElib.fa | grep -v "LTR/unknown" | grep -v "\#unknown" | sed 's/>//' > tmp.idx
seqkit grep -f tmp.idx -n panEDTA.TElib.fa -o TE_known.fa
