#!/bin/bash
#$ -S /bin/bash
#$ -N E3_C.merge  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-38

#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH

#ARRAY
NAME=$(sed "${SGE_TASK_ID}q;d" bamlist)

#script 
samtools merge ${NAME}.FULL.bam ${NAME}.RG*bam
