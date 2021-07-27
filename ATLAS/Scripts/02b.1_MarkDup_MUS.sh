#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MUS.MarkDups  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48


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
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum
TAIL="RG.bam"

#Set up ARRAY job
#ls *bam | awk -F "." '{print $1}' >> mus.names
NAME=$(sed "${SGE_TASK_ID}q;d" mus.names)

echo "java -jar $PICARD MarkDuplicates \
INPUT=$INPUT/${NAME}.$TAIL \
OUTPUT=$OUTPUT/${NAME}.rmdup.bam \
METRICS_FILE=$OUTPUT/${NAME}.dup.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true" >> 02b.1_MarkDup.log


time java -jar $PICARD MarkDuplicates \
INPUT=$INPUT/${NAME}.$TAIL \
OUTPUT=$OUTPUT/${NAME}.rmdup.bam \
METRICS_FILE=$OUTPUT/${NAME}.dup.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true
