#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 06/07/2021 11:20
##################################

# Creates submission script to check quality of demultiplexed raw reads using fastqc


SHAREDFOLDER=/SAN/ugi/LepGenomics
SCRIPTS=VelocityPipeline/wrapper
SPECIES=C3_Aricia_agestis
INFILE=00_raw_reads_museum
OUTFILE=00_raw_reads_museum


$SHAREDFOLDER/$SCRIPTS/00_parallel_fastqc_UCL.sh \
-i $SHAREDFOLDER/$SPECIES/$INFILE \
-o $SHAREDFOLDER/$SPECIES/$OUTFILE \
-t 1 -m 32 -v 32;
