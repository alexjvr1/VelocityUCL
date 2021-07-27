#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 12/07/2019 05:49 
###########################################

## Concatenates fastq reads from a list of file names
## Checks that sample order is the same in both input files. 
## Writes concatenated samples to 00_raw_reads_museum_FINAL
## Moves samples that were sequenced once to 00_raw_reads_museum_FINAL
### Manually check log files that the correct samples were concatenated together. 

#$ -S /bin/bash ##shell to be used
#$ -N C3.mus.concat.cutadapt  ##job name
#$ -l tmem=16G #amount of memory you whish to request
#$ -l h_vmem=16G #enforced amount of shell storage
#$ -l h_rt=1:00:00 ##wall time
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-66  ##change this to the number of samples x 2. In most cases we resequenced 33 indivs per species x 2 (for Read 1 and Read 2). 

#run job in working directory
cd $SGE_O_WORKDIR

#Define variables
SPECIESDIR=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus
PATH1=00_raw_reads_museum
PATH2=00_raw_reads_museum2
OUTPUT=$SPECIESDIR/00_raw_reads_museum_FINAL

##Create output directory if it doesn't exist
if 
  ! ls $OUTPUT; 
then 
  mkdir -p $OUTPUT; 
fi


#Find all samples that were sequenced twice and move mus1 samples to the expected input folder

ls 00_raw_reads_museum/ALLSAMPLES/*R1*gz |awk -F "/" '{print $NF}'|awk -F "_" '{print $1}' > museum1.names
ls 00_raw_reads_museum2/ALLSAMPLES/*R1*gz |awk -F "/" '{print $NF}' |awk -F "_" '{print $1}' > museum2.names

cat museum1.names museum2.names |sort |uniq --repeated > museum1.toconcat

for i in $(cat museum1.toconcat); do cp 00_raw_reads_museum/ALLSAMPLES/$i*gz 00_raw_reads_museum/; done



#create files with sample names listed for the 33 samples
##R1
ls $SPECIESDIR/$PATH2/*R1*gz | awk -F "/" '{print $NF}' > samplenames.museum2.R1
ls $SPECIESDIR/$PATH1/*R1*gz | awk -F "/" '{print $NF}' > samplenames.museum1.R1

##R2
ls $SPECIESDIR/$PATH2/*R2*gz | awk -F "/" '{print $NF}' > samplenames.museum2.R2
ls $SPECIESDIR/$PATH1/*R2*gz | awk -F "/" '{print $NF}' > samplenames.museum1.R2


##Concatenate the two input files to create one input file
cat samplenames.museum1.R1 samplenames.museum1.R2 > mus1_toconcat
cat samplenames.museum2.R1 samplenames.museum2.R2 > mus2_toconcat

#Check that the sample numbers and sample order is the same for both input files
cat mus1_toconcat | awk -F '_' '{print $1}' > START1
cat mus2_toconcat | awk -F '_' '{print $1}' > START2

cat mus1_toconcat | awk -F '_' '{print $NF}' > TAIL1
cat mus2_toconcat | awk -F '_' '{print $NF}' > TAIL2


if 
  ! cmp START1 START2 >/dev/null 2>&1; then 
  echo "Sample order is not the same in mus1_toconcat and mus2_toconcat"
  exit 1
fi


if 
  ! cmp TAIL1 TAIL2 >/dev/null 2>&1; then 
  echo "R1 and R2 order must be the same in mus1_toconcat and mus2_toconcat"
  exit 1
fi


#Create ARRAY job
NAMEMUS1=$(sed "${SGE_TASK_ID}q;d" mus1_toconcat)
NAMEMUS2=$(sed "${SGE_TASK_ID}q;d" mus2_toconcat)



##Concat fastq files

echo "[concatenating] ${NAMEMUS1} and ${NAMEMUS2}"
printf "\n"
echo "time cat $SPECIESDIR/$PATH1/${NAMEMUS1} $SPECIESDIR/$PATH2/${NAMEMUS2} > $OUTPUT/${START1}.concat.fastq.gz" >> concat.mus.log
time cat $SPECIESDIR/$PATH1/${NAMEMUS1} $SPECIESDIR/$PATH2/${NAMEMUS2} > $OUTPUT/${START1}.concat.fastq.gz



##Move all the samples that were sequenced once to the output directory

ls $SPECIESDIR/$PATH1/*R1*gz | awk -F "/" '{print $1}' | awk -F "_" '{print $1}' > museum1.names
ls $OUTPUT/*R1*gz | awk -F "/" '{print $1}' |
awk -F "_" '{print $1}' > museumconcat.names


diff museum1.names museumconcat.names | grep '^<' | sed 's/^<\ //'> museum1.tomove

wc -l museum1.tomove > count
echo "number of samples to move:" $count >> concat.mus.log

for i in $(cat museum1.tomove); do cp $SPECIESDIR/$PATH1/ALLSAMPLES/$i*gz $OUTPUT; done


