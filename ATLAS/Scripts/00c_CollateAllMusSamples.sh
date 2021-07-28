#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 12/07/2019 05:49 
###########################################

## Cmd3 for concatenating all samples sequenced twice
## Moves samples that were sequenced once to 00_raw_reads_museum_FINAL


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


##Move all the samples that were sequenced once to the output directory

ls $SPECIESDIR/$PATH1/ALLSAMPLES/*R1*gz | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}' > museum1.names
ls $OUTPUT/*R1*gz | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}' > museumconcat.names

ls $SPECIESDIR/$OUTPUT/*gz | awk -F "/" '{print $NF}' | awk '{print substr($0,14)}' | awk -F ".concat" '{print $1}' | sort |uniq > TAILS


diff museum1.names museumconcat.names | grep '^<' | sed 's/^<\ //'> museum1.tomove

wc -l museum1.tomove > count
echo "number of samples to move:" $count >> concat.mus.log

while read TAIL1 <&1; do for i in $(cat museum1.tomove); do cp $SPECIESDIR/$PATH1/ALLSAMPLES/$i$TAIL1 $OUTPUT; done; done 1<TAILS
