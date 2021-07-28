#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 28/07/2021 05:49 
###########################################

## Cmd2 to concatenate all samples that were sequenced twice

#$ -S /bin/bash ##shell to be used
#$ -N C3.mus.concat.cutadapt  ##job name
#$ -l tmem=16G #amount of memory you whish to request
#$ -l h_vmem=16G #enforced amount of shell storage
#$ -l h_rt=1:00:00 ##wall time
#$ -j y  #concatenates error and output files (with prefix job1)


#run job in working directory
cd $SGE_O_WORKDIR

#Define variables
SPECIESDIR=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus
PATH1=00_raw_reads_museum
PATH2=00_raw_reads_museum2
OUTPUT=$SPECIESDIR/00_raw_reads_museum_FINAL


##Concat fastq files

echo "[concatenating] $1 and $2" >> concat.mus.log
printf "\n"

echo "while read NAME1 <&MUS1 && read NAME2 <&MUS2; do cat $SPECIESDIR/$PATH1/$NAME1 $SPECIESDIR/$PATH2/$NAME2 > $OUTPUT/$NAME1.concat.fastq.gz; done" >> concat.mus.log
time while read NAME1 <&1 && read NAME2 <&2; do cat $SPECIESDIR/$PATH1/$NAME1 $SPECIESDIR/$PATH2/$NAME2 > $OUTPUT/$NAME1.concat.fastq.gz; done 1<mus1_toconcat 2<mus2_toconcat
