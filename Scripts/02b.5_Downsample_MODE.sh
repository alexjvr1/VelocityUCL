#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODE.Downsample  ##job name
#$ -l tmem=8G #RAM
#$ -l h_vmem=8G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
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
SPECIES=E3_Aphantopus_hyperantus
FOLDER=E3_SubsetTests
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern_exp
OUTPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern_exp_Downsampled



paste mode.names.100 mode.prop.downsample mode.propnames | while IFS=$'\t' read -r var1 var2 var3; do java -jar $PICARD DownsampleSam \
I=$INPUT/"$var1" \
O=$OUTPUT/"$var1".fullLR75.Downsampled."$var3".bam \
P="$var2"; done
