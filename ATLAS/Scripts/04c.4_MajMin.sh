#!/bin/bash
#$ -S /bin/bash
#$ -N Esp.MajMin  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=20:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)



#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=Spanish_Pararge_aegeria
INPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS
GLFLIST=glflist
TASK=majorMinor


#Run analysis
time $ATLAS task=$TASK glf=$INPUT/$GLFLIST
