#!/bin/bash
#$ -S /bin/bash
#$ -N C3.MODE.Downsample  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=8:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)


#Downsample indviduals to a series of depths
#Takes 3 input files of the same length:
# mode.names.100 = a list of sample names, repeated for the series of subsampling. 
# eg. if a list of samples needs to be resampled 10 times: 
# seq 1 10 | xargs -Inone cat mode.names.10 > mode.names.100
# mode.prop.downsample = the proportion of data to sample for each individual
# mode.propnames = a name for the output. I've used the expected depth after downsampling


#Run on working directory
cd $SGE_O_WORKDIR 

#Call software
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=C3_Aricia_agestis_FP
REF=$SHAREDFOLDER/$SPECIES/RefGenome2/GCA_905147365.1_ilAriAges1.1_genomic.fna 
INPUT=$SHAREDFOLDER/$SPECIES/04_Downsampled_MODE
OUTPUT=$SHAREDFOLDER/$SPECIES/04_Downsampled_MODE



paste mode.names.100 mode.prop.downsample.10 mode.propnames.10 | while IFS=$'\t' read -r var1 var2 var3; do java -jar $PICARD DownsampleSam \
I=$INPUT/"$var1" \
O=$OUTPUT/"$var1".fullLR75.Downsampled."$var3".bam \
P="$var2"; done
