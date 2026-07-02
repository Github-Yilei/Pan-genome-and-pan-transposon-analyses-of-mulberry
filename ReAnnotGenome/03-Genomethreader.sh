#!/bin/bash
#PBS -N gth
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node22:ppn=10
#PBS -l mem=60GB




################################
### Project information ###
cd /data_group/wuyilei/Project/ReAnnotGenome/03-Genomethreader
SampleList=/data_group/wuyilei/Project/ReAnnotGenome/samples.list

################################
### Main parameters ###
IFS=$'\r\n' GLOBIGNORE='*' command eval  'ARRAY=($(cat $SampleList))'
sampleID=`echo ${ARRAY[(($PBS_ARRAY_INDEX-1))]}`

source activate zyp_gemoma_env
# tblastn
makeblastdb -in ../${sampleID}.genome.fasta.new.masked -parse_seqids -dbtype nucl -out ${sampleID}.db
tblastn -query all.pep.fa -out all_pep.blast -db ${sampleID}.db -outfmt 6 -evalue 1e-5 -num_threads 20 -qcov_hsp_perc 50.0 -num_alignments 5

# extract pprotein id
awk '{print $1}' all_pep.blast | sort | uniq >unique_pep.list

# extact uniq protein 
seqkit grep -f unique_pep.list all.pep.fa -o ${sampleID}_unique_pep.fa


gth -genomic ../${sampleID}.genome.fasta.new.masked  -protein ${sampleID}_unique_pep.fa -intermediate -gff3out -o ${sampleID}.gth_homology.gff -force
python genomeThreader_to_evm_gff3.py --input_file ${sampleID}.gth_homology.gff --output ${sampleID}.genomethreader.gff3

