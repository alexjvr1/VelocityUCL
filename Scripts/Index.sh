#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.Index  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)

#Script to index a list of bamfiles specified in an input textfile

#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
samtools=/share/apps/genomics/samtools-1.9/bin/samtools

#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
FOLDER=E3_SubsetTests
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern_Downsampled/modc.Downsampled.merged.bamlist
OUTPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern_Downsampled
TASK=theta


#Run analysis
while IFS=  read -r line; do $samtools index $line; done < $INPUT
