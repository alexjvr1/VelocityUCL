#!/bin/bash
#$ -S /bin/bash
#$ -N ESP.ATLAS.GlobDiv  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=5:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-25

#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=Spanish_Pararge_aegeria
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCF_905163445.1_ilParAegt1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped/FINALBAMS/SplitMerge
OUTPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS
RECAL=$SHAREDFOLDER/$SPECIES/02a_mapped/FINALBAMS/SplitMerge/recalFiles
PMD=$SHAREDFOLDER/$SPECIES/02a_mapped/FINALBAMS/SplitMerge
TASK=theta

NAME=$(sed "${SGE_TASK_ID}q;d" bamlist.rpts)


#Run analysis
#while IFS=  read -r line; do $ATLAS task=$TASK bam=$line regions=Regions.LR75.bam > $line.theta  -l avx2=yes; done < $INPUT
$ATLAS task=$TASK bam=$INPUT/${NAME}.bam recal=$RECAL/${NAME}_recalibrationEM.txt pmdFile=$PMD/${NAME}_PMD_input_Empiric.txt > ${NAME}.full.theta
