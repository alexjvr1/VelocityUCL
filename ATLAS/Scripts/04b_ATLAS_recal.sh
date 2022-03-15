#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MUS.RG1recal.CHRSonly  ##job name
#$ -l tmem=64G #RAM
#$ -l h_vmem=64G #enforced limit on shell memory usage
#$ -l h_rt=24:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)


#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
TASK=recal
BAM=RG1.merged.bam
#PMDfile only required for Museum samples
PMD=MUS_Empiric.txt   
INVARIANT=A.hyperantus_LR76only.bed
POOLRGS=mergeTheseRGs.txt
RG=RG1

#Script

time $ATLAS task=$TASK bam=$BAM pmdFile=$PMD regions=$INVARIANT minDepth=2 poolReadGroups=$POOLRGS >> $RG_recal_full.txt
