#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.LocalRealignment  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
#$ -l tscratch=20G
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-38

#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/jdk1.8.0_131/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/jdk1.8.0_131/lib:$LD_LIBRARY_PATH


# Define variables
USERNAME=ajansen
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged

GenomeAnalysisTK=/share/apps/genomics/GenomeAnalysisTK-3.8.1.0/GenomeAnalysisTK.jar


#Set up ARRAY job
#ls $INPUT/*rmdup.bam |awk -F "/" '{print $NF}' | awk -F "." '{print $1}' > modc.names 
NAME=$(sed "${SGE_TASK_ID}q;d" modc.names)


#Set up	scratch	space
mkdir -p /scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID
TMP_DIR=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID
TMPDIR=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID


# Identify targets to realign
java -Xmx4g -Xms4g -Djava.io.tmpdir=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID -jar $GenomeAnalysisTK -T RealignerTargetCreator \
-R $REF \
-o $OUTPUT/${NAME}.intervals \
-I $INPUT/${NAME}.rmdup.bam


# use IndelRealigner to realign the regions found in the RealignerTargetCreator step
java -Xmx4g -Xms4g -Djava.io.tmpdir=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID -jar $GenomeAnalysisTK -T IndelRealigner \
-R $REF \
-targetIntervals $INPUT/${NAME}.intervals \
-I $INPUT/${NAME}.rmdup.bam \
-o $OUTPUT/${NAME}.realn.bam


function finish {
    rm -rf /scratch0/ajansen/$JOB_ID.$SGE_TASK_ID
}

trap finish EXIT ERR INT TERM
