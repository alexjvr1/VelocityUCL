#!/bin/bash
##################################
# Kaitlyn Louth under supervision of Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 08/07/2021 14:24
##################################

# Creates submission script to use cutadapt to remove adapters from demultiplexed Illumina li$
# Poly-A and -T filters have been removed as we decided that these would be dropped due to ma$
# Additional filters included in Lymantria monacha dataset
# -fwad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -fwad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# -rvad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -rvad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

SHAREDFOLDER=/SAN/ugi/LepGenomics
SCRIPTS=VelocityPipeline/tools
SPECIES=C3_Aricia_agestis
INFILE=00_raw_reads_modern_exp
OUTFILE=local/01a_filtered_reads_modern_exp

$SHAREDFOLDER/$SCRIPTS/01a_parallel_cutadapt_UCL.sh \
-i $SHAREDFOLDER/$SPECIES/$INFILE \
-o $OUTFILE -t 1 -m 8 -ph 33 \
-fwad1 AGATCGGAAGAGCACACGTCTGAACTCCAGTC -rvad1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-minl 20 -phredq 20;
