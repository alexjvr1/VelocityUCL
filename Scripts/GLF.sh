#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.GLF  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=3:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-38


#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
FOLDER=E3_SubsetTests
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS/MODC/modc.bamlist
OUTPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS/MODC
TASK=GLF

NAME=$(sed "${SGE_TASK_ID}q;d" modc.bamlist)

#Run analysis
$ATLAS task=$TASK bam=${NAME} 
