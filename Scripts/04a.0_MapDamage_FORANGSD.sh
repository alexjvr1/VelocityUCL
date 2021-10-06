#!/bin/bash
# -----------------------------------------

# (c) Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified: 03/10/2018 15:39:48

# Description:
# MapDamage to assess bias DNA damage in museum vs modern data and rescale mapping quality in .bam files accordingly. 
# Each samples runs in 30min - 5hours. The script completes a species museum and modern bam files in 8 hours on bluecrystal

#$ -S /bin/bash
#$ -N E3.MapDmg.mus  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-10

#Run on working directory
cd $SGE_O_WORKDIR


# Software
##python
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH

##R
export PATH=/share/apps/R-4.0.3/bin:$PATH

##mapDamage
mapDamage="/share/apps/python-3.8.5-shared/bin/mapDamage"


# Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_SubsetTests
REF=$SHAREDFOLDER/E3_Aphantopus_hyperantus/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum_FORANGSD
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum_FORANGSD
TAIL="LR75.subset.bam"


#Create Array
NAME=$(sed "${SGE_TASK_ID}q;d" mus.names)


##Script
echo "time $mapDamage --merge-libraries -i $INPUT/${NAME}.$TAIL -d $OUTPUT -r $REF --rescale --single-stranded"
time $mapDamage --merge-libraries -i $INPUT/${NAME}.$TAIL -d $OUTPUT/${NAME} -r $REF --rescale --single-stranded  
