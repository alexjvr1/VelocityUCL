#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MUS.RG1split  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on memory shell usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

#run job in working directory
cd $SGE_O_WORKDIR 

#software
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar

#variables
SHAREDPATH=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus
POP=02a_mapped_museum


#ARRAY
NAME=$(sed "${SGE_TASK_ID}q;d" MUS.RG2.bamlist)


#script
java -Xmx6g -Xms6g -jar $PICARD FilterSamReads I=${NAME} O=${NAME}.RG2.bam READ_LIST_FILE=${NAME}.RG2.list FILTER=includeReadList
