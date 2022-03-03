#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MUS.RG1recal  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=2:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)


#Run on working directory
cd $SGE_O_WORKDIR 



bedtools=/share/apps/genomics/bedtools-2.30.0/bin/bedtools

$bedtools intersect -abam RG1.merged.bam -b A.hyperantus_LR76only.bed > RG1.invariant.merged.bam
