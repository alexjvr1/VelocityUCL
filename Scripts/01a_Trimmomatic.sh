#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# Last Modified 05/07/2021
##################################
#v1 -modified path to wrapper script

# Generates submission script for first part of the trimming process using Trimmomatic. 
# Change the job name (-N), and the input and output directories

/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/parallel_trimmomatic.sh \
-i 00_raw_reads_museum/TEST \
-o TrimmomaticTest \
-n 1 -t 8 -m 4 -l 16 -ph 33 -c 150 -hc 0 -N Trim.mus \
-ad /SAN/ugi/LepGenomics/Software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa \
-illclip '2:30:8:1:true' -le 20 -tr 20 -sw '4:20' -minl 20 -avg 20;
