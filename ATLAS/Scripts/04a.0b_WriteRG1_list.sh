#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MUS.findRG1  ##job name
#$ -l tmem=8G #RAM
#$ -l h_vmem=8G #enforced limit on memory shell usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

#run job in working directory
cd $SGE_O_WORKDIR 


##Software
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH


#Variables
SHAREDPATH=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus
POP=02a_mapped_museum
RG1="M_K00210:144"

#ARRAY
NAME=$(sed "${SGE_TASK_ID}q;d" MUS.bamlist)


#Script
samtools view $SHAREDPATH/$POP/${NAME} | grep $RG1 |awk '{print $1}' >> $SHAREDPATH/$POP/${NAME}.RG1.list
