#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 24/10/2018 14:28
##################################

# Creates submission script to use cutadapt to remove adapters from demultiplexed Illumina libraries. 
# Poly-A and -T filters have been removed as we decided that these would be dropped due to mapping quality by bwa mem. 
# Additional filters included in Lymantria monacha dataset
# -fwad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -fwad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA   
# -rvad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -rvad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA



/panfs/panasas01/bisc/aj18951/bristol-velocity/AJvR_VelocityPipeline/tools/01a_parallel_cutadapt_bluecp3.sh \
-i 00_raw_reads_museum1 \
-o 01a_museum_cutadapt_reads -n 1 -t 8 -m 8 -ph 33 \
-fwad1 AGATCGGAAGAGCACACGTCTGAACTCCAGTC -rvad1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-minl 20 -phredq 20;
