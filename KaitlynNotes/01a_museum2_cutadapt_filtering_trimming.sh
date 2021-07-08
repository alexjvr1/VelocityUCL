#!/bin/bash
#################################
# Kaitlyn Louth under supervision of Alexandra Jansen Van Rensburg
# alexjvr@gmail.com 
# Last modified 08/07/2021 13:47
#################################

# Creates submission script to use cutadapt to remove adapters from demultiplexed Illumina libraries.
# Additional filters included in Lymantria monacha dataset
# -fwad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -fwad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# -rvad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -rvad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

SHAREDFOLDER=/SAN/ugi/LepGenomics
SCRIPTS=VelocityPipeline/tools
SPECIES=C3_Aricia_agestis
INFILE=00_raw_reads_museum2
OUTFILE=local/01a_museum2_cutadapt_reads

$SHAREDFOLDER/$SCRIPTS/01a_parallel_cutadapt_UCL.sh \
-i $SHAREDFOLDER/$SPECIES/$INFILE \
-o $OUTFILE -t 1 -m 8 -ph 33 \
-fwad1 AGATCGGAAGAGCACACGTCTGAACTCCAGTC -rvad1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-minl 20 -phredq 20;
