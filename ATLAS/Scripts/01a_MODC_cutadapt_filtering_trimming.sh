#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 05/07/2021 15:20
##################################
#v1 - modified path to tools script

# Creates submission script to use cutadapt to remove adapters from demultiplexed Illumina libraries. 
# Poly-A and -T filters have been removed as we decided that these would be dropped due to mapping quality by bwa mem. 
# Additional filters included in Lymantria monacha dataset
# -fwad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -fwad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA   
# -rvad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -rvad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SCRIPTS=VelocityPipeline/tools
SPECIES=E3_Aphantopus_hyperantus
INFILE=00_raw_reads_modern
OUTFILE=01a_MODC_cutadapt_reads

$SHAREDFOLDER/$SCRIPTS/01a_parallel_cutadapt_UCL.sh \
-i $SHAREDFOLDER/$SPECIES/$INFILE \
-o $OUTFILE -t 1 -m 8 -ph 33 \
-fwad1 AGATCGGAAGAGCACACGTCTGAACTCCAGTC -rvad1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-minl 20 -phredq 20;
