#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 05/07/2021 15:14
##################################
#v1 - modified path to wrapper script

# Creates submission script to check quality of demultiplexed raw reads using fastqc

/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/00_parallel_fastqc_bcp3.sh -i 00_raw_reads_museum/ -o /00_raw_reads_museum/ -n 1 -t 8 -m 4;
