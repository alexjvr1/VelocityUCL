#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# Last Modified 05/07/2021
##################################
#v1 -modified path to wrapper script

# Generates submission script for second part of the trimming process using Trimmomatic. 
# Romain found that the Cutadapt step (01a) does not remove all of the adapter sequence. 


/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/parallel_trimmomatic_bluecp3.sh \
-i 01a_modern_cutadapt_filtered_reads \
-o 01b_modern_trimmomatic_filtered_reads \
-n 1 -t 8 -m 4 -ph 33 -c 150 -hc 0 \
-ad /fastdata/bo1revx/Liverpool/useful_files/adaptor_files/NEBNext_adaptors.fa \
-illclip '2:30:8:1:true' -le 20 -tr 20 -sw '4:20' -minl 20 -avg 20;
