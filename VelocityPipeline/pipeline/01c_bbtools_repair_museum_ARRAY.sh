#!/bin/bash
#PBS -N E3.BBRepair  ##job name
#PBS -l nodes=1:ppn=16  #nr of nodes and processors per node
#PBS -l mem=32gb #RAM
#PBS -l walltime=1:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-48

#run job in working directory
cd $PBS_O_WORKDIR 
pwd

#Load modules
module load languages/java-jdk-11.0.3

#Define variables
BBREPAIR=/newhome/bzzjrb/Software/bbmap/repair.sh

NAME1=$(sed "${PBS_ARRAYID}q;d" R1.museum.names.torepair)
NAME2=$(sed "${PBS_ARRAYID}q;d" R2.museum.names.torepair)


##Run BBRepair
PREFIX=`echo ${NAME2} |awk -F "_" '{print $1}'`
echo "$BBREPAIR in=01a_mus.concat_cutadapt_reads/${NAME1} in2=01a_mus.concat_cutadapt_reads/${NAME2} out=01c_musPERepaired/${PREFIX}.R1.repaired.fastq.gz out2=01c_musPERepaired/${PREFIX}.R2.repaired.fastq.gz overwrite=f tossbrokenreads " >> repair.log
time $BBREPAIR in=01a_mus.concat_cutadapt_reads/${NAME1} in2=01a_mus.concat_cutadapt_reads/${NAME2} out=01c_musPERepaired/${PREFIX}.R1.repaired.fastq.gz out2=01c_musPERepaired/${PREFIX}.R2.repaired.fastq.gz overwrite=f tossbrokenreads
