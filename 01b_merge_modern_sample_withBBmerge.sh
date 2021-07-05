#!/bin/bash
#PBS -N D2.BBMap1  ##job name
#PBS -l nodes=1:ppn=16  #nr of nodes and processors per node
#PBS -l mem=32gb #RAM
#PBS -l walltime=1:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-47

#run job in working directory
cd $PBS_O_WORKDIR 
pwd

#Load modules
module load languages/java-jdk-11.0.3

#Define variables
BBMERGE=/newhome/bzzjrb/Software/bbmap/bbmerge.sh

NAME1=$(sed "${PBS_ARRAYID}q;d" R1.names)
NAME2=$(sed "${PBS_ARRAYID}q;d" R2.names)


##Run BBMerge
PREFIX=`echo ${NAME1} | awk -F "_mod" '{print $1}'`
#echo "BBmerge read merging started for $PREFIX" >> bbmerge.log
#echo "---------------" >> bbmerge.log

echo "time $BBMERGE in1=${NAME1} in2=${NAME2} out=$PREFIX.merged.fastq outu1=$PREFIX.unmerged1.fastq outu2=$PREFIX.unmerged2.fastq" >> bbmerge.log
time $BBMERGE in1=${NAME1} in2=${NAME2} out=$PREFIX.merged.fastq outu1=$PREFIX.unmerged1.fastq outu2=$PREFIX.unmerged2.fastq
