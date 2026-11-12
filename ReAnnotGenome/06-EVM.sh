#!/bin/bash
#PBS -N EVM
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node22:ppn=12
#PBS -l mem=100GB




################################
### Project information ###
#SampleList=/data_group/wuyilei/Project/ReAnnotGenome/RNA.list
cd /data_group/wuyilei/Project/ReAnnotGenome/06-EVM/Morus_atropurpurea


################################
### Main parameters ###

source activate py27

/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/partition_EVM_inputs.pl \
        --genome /data_group/wuyilei/Project/ReAnnotGenome/Morus_atropurpurea.genome.fasta.new.masked \
        --gene_predictions ../Morus_atropurpurea.gene_predictions.gff3 \
        --protein_alignments ../Morus_atropurpurea.genomethreader.gff3 \
        --transcript_alignments ../Morus_atropurpurea.transcript_alignments.gff3 \
        --segmentSize 100000 \
        --overlapSize 10000 \
        --partition_listing partitions_list.out



/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/write_EVM_commands.pl \
    --weights /data_group/wuyilei/Project/ReAnnotGenome/06-EVM/Morus_atropurpurea/weights.txt \
    --genome /data_group/wuyilei/Project/ReAnnotGenome/Morus_atropurpurea.genome.fasta.new.masked \
    --gene_predictions ../Morus_atropurpurea.gene_predictions.gff3 \
    --protein_alignments ../Morus_atropurpurea.genomethreader.gff3 \
        --transcript_alignments ../Morus_atropurpurea.transcript_alignments.gff3 \
        --output_file_name evm.out \
    --partitions partitions_list.out >commands.list

/opt/service/miniconda3/envs/r4_py37_env/bin/ParaFly -c commands.list -CPU 12

/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl \
        --partitions partitions_list.out --output_file_name evm.out
/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl \
        --partitions partitions_list.out --output evm.out \
        --genome /data_group/wuyilei/Project/ReAnnotGenome/Morus_atropurpurea.genome.fasta.new.masked
