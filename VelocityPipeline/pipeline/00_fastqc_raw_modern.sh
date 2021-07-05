#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 24/10/2018 14:28
##################################

# Creates submission script to check quality of demultiplexed raw reads using fastqc

/panfs/panasas01/bisc/aj18951/bristol-velocity/AJvR_VelocityPipeline/wrapper/00_parallel_fastqc_bcp3.sh -i 00_raw_reads_museum/ -o /00_raw_reads_museum/ -n 1 -t 8 -m 4;
