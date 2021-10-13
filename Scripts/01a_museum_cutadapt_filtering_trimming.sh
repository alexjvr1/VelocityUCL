#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 06/07/2021 11:22
##################################

# Creates submission script to use cutadapt to remove adapters from demultiplexed Illumina libraries. 
# Poly-A and -T filters have been removed as we decided that these would be dropped due to mapping quality by bwa mem. 
#
# Change INFILE and OUTFILE as needed. Change SPECIES.
# Specify the number of individuals after -t
#
# -fwad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -fwad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA   
# -rvad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -rvad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

SHAREDFOLDER=/SAN/ugi/LepGenomics
SCRIPTS=VelocityPipeline/tools
SPECIES=D3_Pararge_aegeria
INFILE=00_raw_reads_museum1
OUTFILE=01a_filtered_reads_museum1

$SHAREDFOLDER/$SCRIPTS/01a_parallel_cutadapt_UCL.sh \
-i $SHAREDFOLDER/$SPECIES/$INFILE \
-o $SHAREDFOLDER/$SPECIES/$OUTFILE -t 48 -m 8 -ph 33 \
-fwad1 AGATCGGAAGAGCACACGTCTGAACTCCAGTC -rvad1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-minl 20 -phredq 20;
