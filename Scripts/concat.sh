#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 28/07/2021 05:49 
###########################################

## Concatenate all museum samples that were sequenced twice

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
PATH1=00_raw_reads_museum/ALLSAMPLES/
PATH2=00_raw_reads_museum2/ALLSAMPLES
OUTPUT=$SPECIESDIR/00_raw_reads_museum_FINAL
TAIL1=_R1.fastq.gz
TAIL2=_R2.fastq.gz

##Concat fastq files

echo "[concatenating] $1 and $2" >> concat.mus.log
printf "\n"

echo "while read NAME1 <&MUS1 && read NAME2 <&MUS2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL1 $SPECIESDIR/$PATH2/$NAME2$TAIL1 > $OUTPUT/$NAME1_R1.concat.fastq.gz; done" >> concat.mus.log
time while read NAME1 <&1 && read NAME2 <&2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL1 $SPECIESDIR/$PATH2/$NAME2$TAIL1 > $OUTPUT/$NAME1_R1.concat.fastq.gz; done 1<museum1.names 2<museum2.names

echo "while read NAME1 <&MUS1 && read NAME2 <&MUS2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL2 $SPECIESDIR/$PATH2/$NAME2$TAIL2 > $OUTPUT/$NAME1_R2.concat.fastq.gz; done" >> concat.mus.log
time while read NAME1 <&1 && read NAME2 <&2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL2 $SPECIESDIR/$PATH2/$NAME2$TAIL2 > $OUTPUT/$NAME1_R2.concat.fastq.gz; done 1<museum1.names 2<museum2.names
