#!/bin/bash
##################################
# Kaitlyn Louth supervised by Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 08/07/2021 14:40
##################################

# Creates submission script to check quality of demultiplexed raw reads using fastqc

/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/00_parallel_fastqc_UCL.sh
-i 00_raw_reads_museum2/
-o /00_raw_reads_museum2/
-t 8 -m 4 -v 4;
