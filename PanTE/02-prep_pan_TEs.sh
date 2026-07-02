#!/usr/bin/env bash
path=/opt/service/miniconda3/envs/edta_env/share/EDTA/bin/
dir=/data_group/wuyilei/Project/PanTE
genome_list=${dir}/samples.list
threads=16

# Set paths
cd ${dir}

# ln -s ../Genome/*/*mod.EDTA.anno/*.mod.EDTA.RM.out  ./
# ln -s ../Genome/*/*mod.EDTA.TElib.novel.fa ./

# Idenfity full-length TEs
for i in `cat $genome_list`; do
        genome=`basename $(echo $i|awk '{print $1}') 2>/dev/null`

        # skip empty lines
        if [ $genome == '' ]; then
                break
        fi

        echo "Idenfity full-length TEs for genome $genome"
        perl $path/find_flTE.pl $genome.*mod.EDTA.RM.out | \
                awk '{print $10}'| \
                sort| \
                uniq -c |\
                perl -snale 'my ($count, $id) = (split); next if $count < $fl_copy; print $_' -- -fl_copy=$fl_copy | \
                awk '{print $2"#"}' > $genome.mod.EDTA.TElib.fa.keep.list
done

# extract pan-TE library candidate sequences
for i in `cat $genome_list`; do
        genome=`basename $(echo $i|awk '{print $1}') 2>/dev/null`

        # skip empty lines
        if [ $genome == '' ]; then
                break
        fi

        if [ -s "$curatedlib" ]; then

                # a) if --curatedlib is provided
                for j in `cat $genome.mod.EDTA.TElib.fa.keep.list`; do
                        grep $j $genome.*mod.EDTA.TElib.novel.fa;
                done | \
                        perl $path/output_by_list.pl 1 $genome.*mod.EDTA.TElib.novel.fa 1 - -FA > $genome.mod.EDTA.TElib.fa.keep.ori
        else
                # b) if --curatedlib is not provided
                for j in `cat $genome.mod.EDTA.TElib.fa.keep.list`; do
                        grep $j $genome.*mod.EDTA.TElib.fa;
                done | \
                        perl $path/output_by_list.pl 1 $genome.*mod.EDTA.TElib.fa 1 - -FA > $genome.mod.EDTA.TElib.fa.keep.ori
        fi
done
