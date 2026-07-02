#!/bin/bash
#PBS -N PASA
#PBS -o QsubLog.log
#PBS -e QsubLog.err
#PBS -q wuyileiQueue
#PBS -l nodes=node22:ppn=10
#PBS -l mem=100GB




################################
### Project information ###
cd /data_group/wuyilei/Project/ReAnnotGenome/07-Update/Morus_alba

################################
### Main parameters ###


## genome_guided_trinity
source activate dyw_py36_env

Trinity --genome_guided_bam Morus_alba.bam --genome_guided_max_intron 10000 --max_memory 100G --CPU 6 --output Morus_alba_trinity


## PASA update
source activate py27
/opt/service/miniconda3/envs/py27/opt/pasa-2.5.2/misc_utilities/accession_extractor.pl  <Trinity-GG.fasta > tdn.accs
/opt/service/miniconda3/envs/py27/opt/pasa-2.5.2/Launch_PASA_pipeline.pl \
    -c alignAssembly.config -t Trinity-GG.fasta \
    -C -R -g  /data_group/wuyilei/Project/ReAnnotGenome/Morus_alba.genome.fasta.new.masked \
    --ALIGNERS blat --CPU 20  \
   --TDN tdn.accs

/opt/service/miniconda3/envs/py27/opt/pasa-2.5.2/scripts/build_comprehensive_transcriptome.dbi \
    -c alignAssembly.config -t Trinity-GG.fasta \
    --min_per_ID 95 \
    --min_per_aligned 30

/opt/service/miniconda3/envs/py27/opt/pasa-2.5.2/scripts/Load_Current_Gene_Annotations.dbi -c alignAssembly.config \
    -g /data_group/wuyilei/Project/ReAnnotGenome/Morus_alba.genome.fasta.new.masked \
    -P /data_group/wuyilei/Project/ReAnnotGenome/06-EVM/Morus_alba.EVidenceModeler_combined.gff3

/opt/service/miniconda3/envs/py27/opt/pasa-2.5.2/Launch_PASA_pipeline.pl \
    -c annotCompare.config -A -g /data_group/wuyilei/Project/ReAnnotGenome/Morus_alba.genome.fasta.new.masked \
    -t compreh_init_build/compreh_init_build.fasta

