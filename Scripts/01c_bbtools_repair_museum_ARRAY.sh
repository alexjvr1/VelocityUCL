#!/bin/bash
#$ -S /bin/bash
#$ -N E3.BBRepair  ##job name
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
BBREPAIR=/share/apps/genomics/bbmap-38.59/repair.sh
SHAREDPATH=/SAN/ugi/LepGenomics
SPECIESDIR=$SHAREDPATH/E3_Aphantopus_hyperantus
PATH1=01a_museum_cutadapt_reads
PATH2=01c_musPERepaired

# create input files here by: 
#ls $SPECIESDIR/$PATH1/*R1* | awk -F "/" '{print $NF}' > R1.museum.names.torepair
#ls $SPECIESDIR/$PATH1/*R2* | awk -F "/" '{print $NF}' > R2.museum.names.torepair

NAME1=$(sed "${SGE_TASK_ID}q;d" R1.museum.names.torepair)
NAME2=$(sed "${SGE_TASK_ID}q;d" R2.museum.names.torepair)


##Run BBRepair
PREFIX=`echo ${NAME2} |awk -F "_" '{print $1}'`
echo "$BBREPAIR in=$SPECIESDIR/$PATH1/${NAME1} in2=$SPECIESDIR/$PATH1/${NAME2} out=$SPECIESDIR/$PATH2/${PREFIX}.R1.repaired.fastq.gz out2=$SPECIESDIR/$PATH2/${PREFIX}.R2.repaired.fastq.gz overwrite=f tossbrokenreads " >> repair.log
time $BBREPAIR in=$SPECIESDIR/$PATH1/${NAME1} in2=$SPECIESDIR/$PATH1/${NAME2} out=$SPECIESDIR/$PATH2/${PREFIX}.R1.repaired.fastq.gz out2=$SPECIESDIR/$PATH2/${PREFIX}.R2.repaired.fastq.gz overwrite=f tossbrokenreads
