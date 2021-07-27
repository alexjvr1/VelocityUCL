#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MUS_AddRG  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-33

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
TAIL="R1.repaired.fastq.gz.repaired.merged.fastq.gz.bam"

#Set up ARRAY job
ls 02a_mapped_museum/*bam | awk -F "/" '{print $NF}' | awk -F "." '{print $1}' > mus.names 
NAME=$(sed "${SGE_TASK_ID}q;d" mus.names)


##Add readgroups

echo "java -jar $PICARD AddOrReplaceReadGroups \
       I=$INPUT/${NAME}.$TAIL \
       O=$OUTPUT/${NAME}.RG.bam \
       RGID=E3mus \
       RGLB=museum0204 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=${NAME}" >> 02b.0_AddRG.log


time java -jar $PICARD AddOrReplaceReadGroups \
       I=$INPUT/${NAME}.$TAIL \
       O=$OUTPUT/${NAME}.RG.bam \
       RGID=E3mus \
       RGLB=museum0204 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=${NAME}
