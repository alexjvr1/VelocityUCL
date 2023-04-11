#!/bin/bash
#$ -S /bin/bash
#$ -N C1.fix.modc  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-2


#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=C1_Aricia_artaxerxes
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_modern.unmerged
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_modern.unmerged
TAIL="realn.bam"

#Create a list of all the samples to be corrected:
#grep ERROR *validatesam | awk -F "." '{print $1}' |sort |uniq > modc.names.tofix

#Set up ARRAY job
NAME=$(sed "${SGE_TASK_ID}q;d" tofix.names)

#Generic script
#java -jar picard.jar FixMateInformation \ I=input.bam \ O=fixed_mate.bam \ ADD_MATE_CIGAR=true

echo "java -Xmx4g -Xms4g -Djava.io.tmpdir=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID -jar $PICARD FixMateInformation \
INPUT=$INPUT/${NAME}.realn.bam \
OUTPUT=$OUTPUT/${NAME}.fixed.bam \
ADD_MATE_CIGAR=true" >> 02b.4_FixMate.log


time java -Xmx4g -Xms4g -Djava.io.tmpdir=/scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID -jar $PICARD FixMateInformation \
INPUT=$INPUT/${NAME}.realn.bam \
OUTPUT=$OUTPUT/${NAME}.fixed.bam \
ADD_MATE_CIGAR=true 

function finish {
    rm -rf /scratch0/$USERNAME/$JOB_ID.$SGE_TASK_ID
}
