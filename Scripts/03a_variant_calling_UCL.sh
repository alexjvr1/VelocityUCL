#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 05/07/2021 15:56
##################################
#v1 - modified path to wrapper script

# Creates submission script for variant calling using mpileup
# using arrays on BlueCrystal p3

#run job in working directory
cd $SGE_O_WORKDIR

#load your program if it is installed globally or the modules you used to install your program locally (compilers, etc) 
#Specify modules
export PATH=/share/apps/perl-5.30.0/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/perl-5.30.0/lib:$LD_LIBRARY_PATH


SHAREDFOLDER="/SAN/ugi/LepGenomics"
SPECIES="C3_Aricia_agestis"
INPUT="02a_mapped_modern"
OUTPUT="03.2_Variant_calling_modern"
REF="RefGenome/GCA_905147365.1_ilAriAges1.1_genomic.fna"
SAMTOOLS="/share/apps/genomics/samtools-1.9/bin/samtools"
JOBNAME="C3_mod_mpileup"
CALLER="/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/03a_call_SNVs_UCL.sh"

$CALLER -i $SHAREDFOLDER/$SPECIES/$INPUT \
-r $SHAREDFOLDER/$SPECIES/$REF -o $SHAREDFOLDER/$SPECIES/$OUTPUT \
-c c -v 1 -d 0 -s 20 -p 0.05 \
-job $JOBNAME;
