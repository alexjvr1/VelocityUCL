#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.SplitMerge  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
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
FOLDER=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_MODC.unmerged/modc.bamlist
OUTPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_MODC.unmerged/SplitMerge
TASK=splitMerge

NAME=$(sed "${SGE_TASK_ID}q;d" modc.bamlist)


#Run analysis
#while IFS= read -r line; do $ATLAS task=$TASK bam=$line readGroupSettings=RG.txt; done < $INPUT
$ATLAS task=$TASK bam=${NAME} readGroupSettings=RG.txt
