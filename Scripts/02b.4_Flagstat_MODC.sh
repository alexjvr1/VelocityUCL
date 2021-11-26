#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.flagstat  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-35


#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH
SAMTOOLS=/share/apps/genomics/samtools-1.9/bin/samtools


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged
TAIL="realn.bam"


#Set up ARRAY job
#ls *bam | awk -F "." '{print $1}' >> modc.names
NAME=$(sed "${SGE_TASK_ID}q;d" modc.names)


#Flagstat for each bam
time $SAMTOOLS flagstat	$INPUT/${NAME}.$TAIL > $OUTPUT/${NAME}.$TAIL.flagstat
