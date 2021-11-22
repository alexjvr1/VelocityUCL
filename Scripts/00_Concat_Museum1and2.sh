#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# Last Modified 22/11/2021
##################################

# Generates submission script to concatenate museum1 and museum2 samples together. 


/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/concat.sh \
-I 00_raw_reads_museum/ALLSAMPLES \
-i 00_raw_reads_museum2/ALLSAMPLES \
-S E3_Aphantopus_hyperantus \
-o 00_raw_reads_museum_FINAL \
-N E3.concat.mus;
