#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODE.ATLAS.recal  ##job name
#$ -l tmem=64G #RAM
#$ -l h_vmem=64G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)


#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.19/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS/MODE/recal.bamlist
OUTPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS/MODE
TASK=recal


#Run analysis
while IFS=  read -r line; do time $ATLAS task=$TASK bam=$line region=../A.hyperantus_LR76only.bam minDepth=2  maxDepth=4 limitWindows=2 verbose > $line.recal -l avx2=yes; done < $INPUT
