#!/bin/bash
#$ -S /bin/bash 
#$ -N E3.BBmerge  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

#run job in working directory
cd $SGE_O_WORKDIR 
pwd

#Load modules
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH

#Define variables
BBMERGE=/share/apps/genomics/bbmap-38.59/bbmerge.sh
PATH1=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/01b_musPERepaired
PATH2=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/01c_musALL_merged

NAME1=$(sed "${SGE_TASK_ID}q;d" R1.museum.names.repaired)
NAME2=$(sed "${SGE_TASK_ID}q;d" R2.museum.names.repaired)


##Run BBMerge
PREFIX=`echo ${NAME1} | awk -F "_" '{print $1}'`
#echo "BBmerge read merging started for $PREFIX" >> bbmerge.log
#echo "---------------" >> bbmerge.log

echo "time $BBMERGE in1=$PATH1/${NAME1} in2=$PATH1/${NAME2} out=$PATH2/$PREFIX.repaired.merged.fastq.gz outu1=$PATH2/$PREFIX.repaired.unmerged1.fastq.gz outu2=$PATH2/$PREFIX.repaired.unmerged2.fastq.gz" >> bbmerge.log
time $BBMERGE in1=$PATH1/${NAME1} in2=$PATH1/${NAME2} out=$PATH2/$PREFIX.repaired.merged.fastq.gz outu1=$PATH2/$PREFIX.repaired.unmerged1.fastq.gz outu2=$PATH2/$PREFIX.repaired.unmerged2.fastq.gz
