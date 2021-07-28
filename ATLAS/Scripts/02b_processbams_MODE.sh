#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODE.processbams  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)

###Script to process mapped bam files in series
##1. Add RG information
##2. Remove duplicates
##3. Local realignment with GATK3.8
##4. Validate sam file
####################################

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
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_modern_exp
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_modern_exp
TAIL1=`ls 02a_mapped_modern_exp/*bam | awk -F "/" '{print $NF}' | awk '{print substr($0,14)}' | head -n 1`
TAIL2="RG.bam"


############################
#1. Add ReadGroups
############################

#Create input file
ls 02a_mapped_modern_exp/*bam | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}' > mode.names 

#Write command to log file

echo "while read NAME1 <&1; do java -jar $PICARD AddOrReplaceReadGroups \
       I=$INPUT/$NAME1.$TAIL1 \
       O=$OUTPUT/$NAME1.RG.bam \
       RGID=E3mode \
       RGLB=modern4 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=$NAME1; done 1<mode.names" >> 02b.0_AddRG.log


time while read NAME1 <&1; do java -jar $PICARD AddOrReplaceReadGroups \
       I=$INPUT/$NAME1$TAIL1 \
       O=$OUTPUT/$NAME1.RG.bam \
       RGID=E3mode \
       RGLB=modern4 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=$NAME1; done 1<mode.names 



#############################
#2. Mark Duplicates
#############################

#Create list of bam files to process
ls $INPUT/*RG.bam | awk -F "." '{print $1}' >> mode.names

#Write command to log file
echo "while read NAME2 <&2; do 
java -jar $PICARD MarkDuplicates \
INPUT=$INPUT/$NAME2.$TAIL2 \
OUTPUT=$OUTPUT/$NAME2.rmdup.bam \
METRICS_FILE=$OUTPUT/$NAME2.dup.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true; done 2<mode.names" >> 02b.ProcessBams.log

#Run command
time while read NAME2 <&2; do java -jar $PICARD MarkDuplicates \
INPUT=$INPUT/$NAME2.$TAIL2 \
OUTPUT=$OUTPUT/$NAME2.rmdup.bam \
METRICS_FILE=$OUTPUT/$NAME2.dup.txt \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=SILENT \
CREATE_INDEX=true; done 2<mode.names



#############################
#3. Local realignment
#############################

#We import an old version of GATK and an older version of java for this step
export PATH=/share/apps/jdk1.8.0_131/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/jdk1.8.0_131/lib:$LD_LIBRARY_PATH
GenomeAnalysisTK=/share/apps/genomics/GenomeAnalysisTK-3.8.1.0/GenomeAnalysisTK.jar

#create input file
ls $INPUT/*rmdup.bam |awk -F "/" '{print $NF}' | awk -F "." '{print $1}' > mode.names 

# Identify targets to realign
time while read NAME3 <$3; do java -jar $GenomeAnalysisTK -T RealignerTargetCreator \
-R $REF \
-o $OUTPUT/$NAME3..intervals \
-I $INPUT/$NAME3.rmdup.bam; done 3<mode.names


# use IndelRealigner to realign the regions found in the RealignerTargetCreator step
while read NAME3 <&3; do java -jar $GenomeAnalysisTK -T IndelRealigner \
-R $REF \
-targetIntervals $INPUT/$NAME3.intervals \
-I $INPUT/$NAME3.rmdup.bam \
-o $OUTPUT/$NAME3.realn.bam; done 3<mode.names



#############################
#4. Validate bam file
#############################

#Import the newest version of java again
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar

#Define variables
TAIL4="realn.bam"

#create input file
ls $INPUT/*realn.bam |awk -F "/" '{print $NF}' | awk -F "." '{print $1}' > mode.names

#Write command to log file
echo "while read NAME4 <&4; do java -jar $PICARD ValidateSamFile \
INPUT=$INPUT/$NAME4.$TAIL4 \
OUTPUT=$OUTPUT/$NAME4.validatesam
MODE=SUMMARY; done 4<mode.names" >> 02b.ProcessBams.log

#Run validatesamfile
time while read NAME4 <&4; do java -jar $PICARD ValidateSamFile \
INPUT=$INPUT/$NAME4.$TAIL4 \
OUTPUT=$OUTPUT/$NAME4.validatesam \
MODE=SUMMARY; done 4<mode.names
