#!/bin/bash
#PBS -N WPABWA2  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=1:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-6

#run job in working directory
cd $PBS_O_WORKDIR 
pwd
#cd WPA #uncomment when running locally
#pwd

#Load modules
module load apps/bwa-0.7.15

#Define variables

RefSeq=GCF_000002315.5_GRCg6a_genomic.fna
total_files=`find demultiplexed/ -name '*.fastq.gz' | wc -l`
NAME1=$(sed "${PBS_ARRAYID}q;d" R1.names)
NAME2=$(sed "${PBS_ARRAYID}q;d" R2.names)

#arr=("demultiplexed/"*.fastq.gz)
echo "mapping started" >> map.log
echo "---------------" >> map.log

##Check if Ref Genome is indexed by bwa
if [[ ! $RefSeq.fai ]]
then 
	echo $RefSeq" not indexed. Indexing now"
	bwa index $RefSeq
else
	echo $RefSeq" indexed"
fi


##Map using array

sample_name=`echo ${NAME1} | awk -F "." '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time bwa mem $RefSeq demultiplexed/${NAME1} demultiplexed/${NAME2} > demultiplexed/${NAME1}.sam" >> map.log
time bwa mem $RefSeq demultiplexed/${NAME1} demultiplexed/${NAME2} > demultiplexed/${NAME1}.sam
