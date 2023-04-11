#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# Last Modified 15/11/2021
##################################
#v1 -modified path to wrapper script
#v1 -Change Trimmomatic to be the first adapter removal step. Changed output directory

# Generates submission script for second part of the trimming process using Trimmomatic. 
# Romain found that the Cutadapt step (01a) does not remove all of the adapter sequence. 


/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/parallel_trimmomatic_UCL.sh \
-i 00_raw_reads_museum_FINAL \
-o 01a_Trimmomatic_museum \
-n 1 -t 8 -m 16 -ph 33 -c 150 -hc 0 \
-ad /SAN/ugi/LepGenomics/Software/Trimmomatic-0.39/adapters/Velocity.fa \
-illclip '2:30:8:1:true' -le 20 -tr 20 -sw '4:20' -minl 20 -avg 20;
