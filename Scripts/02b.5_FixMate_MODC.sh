#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.Picard  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-5


#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged
TAIL="realn.bam"

#Set up ARRAY job
#ls *bam | awk -F "." '{print $1}' >> modc.names
NAME=$(sed "${SGE_TASK_ID}q;d" modc.names.tofix)

#Generic script
#java -jar picard.jar FixMateInformation \ I=input.bam \ O=fixed_mate.bam \ ADD_MATE_CIGAR=true

echo "java -jar $PICARD FixMateInformation \
INPUT=$INPUT/${NAME}.realn.bam \
OUTPUT=$OUTPUT/${NAME}.fixed.bam \
ADD_MATE_CIGAR=true" >> 02b.4_FixMate.log


time java -jar $PICARD FixMateInformation \
INPUT=$INPUT/${NAME}.realn.bam \
OUTPUT=$OUTPUT/${NAME}.fixed.bam \
ADD_MATE_CIGAR=true 
