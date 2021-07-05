#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 12/07/2019 05:49 
###########################################

## Concatenates fastq reads from a list of file names
## R1 and R2 are processed separately


#PBS -N D1.AH.concat.R2  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=10:00:00 ##wall time.
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-33  #change this to the number of samples resequenced. In most cases this is 33. 

#run job in working directory
cd $PBS_O_WORKDIR

#Define variables

SPECIESDIR="E3_Aphantopus_hyperantus_2020"

#create files with sample names listed for the 33 samples
#ls ../01a_museum2_cutadapt_reads/*R1*gz > samplenames.museum2.R1
#ls ../01a_museum2_cutadapt_reads/*R2*gz > samplenames.museum2.R2
#ls ../01a_museum_cutadapt_reads/*R1*gz >> samplenames.museum1.R1
#ls ../01a_museum_cutadapt_reads/*R2*gz >> samplenames.museum1.R2
#remove all the extra sample names from the last two files. Make sure the sample names are in the same order in all files.
#remove path before sample names

#NAMEMUS1R1=$(sed "${PBS_ARRAYID}q;d" samplenames.museum1.R1)
NAMEMUS1R2=$(sed "${PBS_ARRAYID}q;d" samplenames.museum1.R2)
#NAMEMUS2R1=$(sed "${PBS_ARRAYID}q;d" samplenames.museum2.R1)
NAMEMUS2R2=$(sed "${PBS_ARRAYID}q;d" samplenames.museum2.R2)


##Concat R2 fastq files
echo "[concatenating] ${NAMEMUS1R2} and ${NAMEMUS2R2}"
printf "\n"

echo "time cat /newhome/aj18951/$SPECIESDIR/01a_museum_cutadapt_reads/${NAMEMUS1R2} /newhome/aj18951/$SPECIESDIR/01a_museum2_cutadapt_reads/${NAMEMUS2R2} > ${NAMEMUS2R2}.concat.fastq.gz" >> concat.mus.R2.log
time cat /newhome/aj18951/$SPECIESDIR/01a_museum_cutadapt_reads/${NAMEMUS1R2} /newhome/aj18951/$SPECIESDIR/01a_museum2_cutadapt_reads/${NAMEMUS2R2} > ${NAMEMUS2R2}.concat.fastq.gz
