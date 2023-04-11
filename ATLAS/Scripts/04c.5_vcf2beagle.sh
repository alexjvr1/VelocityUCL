#!/bin/bash
#$ -S /bin/bash
#$ -N Esp.vcf2beagle  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=5:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)


#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=Spanish_Pararge_aegeria
INPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS
VCF=ATLAS_majorMinor_majorMinor.vcf.gz
TASK=VCFToBeagle
SAMPLES=bamlist


#Run analysis

$ATLAS task=$TASK vcf=$VCF samples=$SAMPLES
