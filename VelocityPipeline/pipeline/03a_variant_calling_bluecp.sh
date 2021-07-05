#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 23/10/2018 09:35
##################################

# Creates submission script for variant calling using mpileup
# using arrays on BlueCrystal p3

#run job in working directory
cd $PBS_O_WORKDIR

#load your program if it is installed globally or the modules you used to install your program locally (compilers, etc) 
#Specify modules

~/bristol-velocity/AJvR_VelocityPipeline/wrapper/03a_call_SNVs_bluecp3.sh -i ~/2a_WPA/variant_calling_input_folder/ -r ~/2a_WPA/GCF_000002315.5_GRCg6a_genomic.fna -o 2a_WPA/03_variants -c c -v 1 -d 0 -s 20 -p 0.05 -module1 apps/samtools-1.8 -job newtest;
