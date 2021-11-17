#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# Last Modified 17/11/2021
##################################
#v1

# Generates submission script for second part of the trimming process using AdapterRemoval
# Romain found that the Cutadapt step (01a) does not remove all of the adapter sequence. We've changed the pipeline to use Trimmomatic followed by AdapterRemoval. 
# AdapterRemoval will merge any overlapping reads as well, which is needed for the museum samples. 


/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/parallel_adapterremoval_UCL.sh \
-i 01a_Trimmomatic_museum \
-o 01b_AdapterRemoval_museum \
-n 1 -t 1 -m 8 -N AdaptRem.mus;
