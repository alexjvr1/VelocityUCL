#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODE_AddRG  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -l tscratch=20G
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-35

#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar


#Define variables
USERNAME=ajansen
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_MODC.unmerged

#Set up ARRAY job
#ls 02a_mapped_modern/*bam | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}' > modc.names 
NAME=$(sed "${SGE_TASK_ID}q;d" modc.names)


#Set up scratch space
mkdir -p /scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID
TMP_DIR=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID
TMPDIR=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID



##Add readgroups

echo "java -Xmx4g -Xms4g -Djava.io.tmpdir=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID -jar $PICARD AddOrReplaceReadGroups \
       I=$INPUT/${NAME}.bam \
       O=$OUTPUT/${NAME}.RG.bam \
       RGID=E3modc \
       RGLB=mod03 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=${NAME}" >> 02b.0_AddRG.log


time java -Xmx4g -Xms4g -Djava.io.tmpdir=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID -jar $PICARD AddOrReplaceReadGroups \
       I=$INPUT/${NAME}.bam \
       O=$OUTPUT/${NAME}.RG.bam \
       RGID=E3modc \
       RGLB=mod03 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=${NAME}

function finish {
    rm -rf /scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID
}

trap finish EXIT ERR INT TERM
