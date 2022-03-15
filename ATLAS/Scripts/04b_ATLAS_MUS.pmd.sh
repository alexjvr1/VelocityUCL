#!/bin/bash
#$ -S /bin/bash
#$ -N D3.MUS.pmd  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.19/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=D3_Pararge_aegeria
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCF_905163445.1_ilParAegt1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02b.2_LR_mus/bamlist
OUTPUT=$SHAREDFOLDER/$SPECIES/
TASK=PMD

##Set up ARRAY
NAME=$(sed "${SGE_TASK_ID}q;d" MUS.RG1only.bamlist)


#Run analysis
#while IFS=  read -r line; do time $ATLAS task=$TASK bam=$line region=../M_hyperantus.phast.1.766.bed minDepth=2 maxDepth=10 window=500000 verbose > $line.recal -l avx2=yes; done < $INPUT

$ATLAS task=$TASK bam=${NAME} fasta=$REF
