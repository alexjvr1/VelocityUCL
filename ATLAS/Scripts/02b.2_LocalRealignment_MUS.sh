#!/bin/bash
#$ -S /bin/bash
#$ -N E3.LocalRealignment  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/jdk1.8.0_131/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/jdk1.8.0_131/lib:$LD_LIBRARY_PATH


# Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum

GenomeAnalysisTK=/share/apps/genomics/GenomeAnalysisTK-3.8.1.0/GenomeAnalysisTK.jar


#Set up ARRAY job
ls $INPUT/*rmdup.bam | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' > mus.names
NAME=$(sed "${SGE_TASK_ID}q;d" mus.names)


# Identify targets to realign
java -jar $GenomeAnalysisTK -T RealignerTargetCreator \
-R $REF \
-o $OUTPUT/${NAME}.intervals \
-I $INPUT/${NAME}.rmdup.bam


# use IndelRealigner to realign the regions found in the RealignerTargetCreator step
java -jar $GenomeAnalysisTK -T IndelRealigner \
-R $REF \
-targetIntervals $INPUT/${NAME}.intervals \
-I $INPUT/${NAME}.rmdup.bam \
-o $OUTPUT/${NAME}.realn.bam
