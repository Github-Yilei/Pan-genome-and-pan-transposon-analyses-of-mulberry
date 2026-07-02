#!/bin/bash
#PBS -N Updata
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node22:ppn=6
#PBS -l mem=100GB




################################
### Project information ###
cd /data_group/wuyilei/Project/ReAnnotGenome/06-EVM

################################
### Main parameters ###

source activate py27

/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/partition_EVM_inputs.pl \
        --genome ../Morus_indica.genome.fasta.new.masked \
        --gene_predictions Morus_indica.augustsus.gff3 \
        --gene_predictions Morus_indica.genemark.gff3 \
        --protein_alignments Morus_indica.genomethreader.gff3 \
        --transcript_alignments Morus_indica.gtf \
        --transcript_alignments  Morus_indica.sqlite.pasa_assemblies.gff3 \
        --segmentSize 100000 \
        --overlapSize 10000 \
        --partition_listing partitions_list.out



/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/write_EVM_commands.pl \
    --weights /data_group/wuyilei/Project/ReAnnotGenome/06-EVM/weights.txt \
    --genome /data_group/wuyilei/Project/ReAnnotGenome/Morus_indica.genome.fasta.new.masked \
    --gene_predictions Morus_indica.augustsus.gff3 \
    --gene_predictions Morus_indica.genemark.gff3 \
    --protein_alignments Morus_indica.genomethreader.gff3 \
    --transcript_alignments Morus_indica.gtf \
    --transcript_alignments  Morus_indica.sqlite.pasa_assemblies.gff3 \
    --output_file_name evm.out \
    --partitions partitions_list.out >commands.list

ParaFly -c commands.list -CPU 20

/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
/opt/service/miniconda3/envs/py27/opt/evidencemodeler-1.1.1/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome ../Morus_indica.genome.fasta.new.masked

find . -regex ".*evm.out.gff3" -exec cat {} \; > Morus_indica.EVidenceModeler_combined.gff3