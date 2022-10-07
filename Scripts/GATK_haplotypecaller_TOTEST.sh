#!/bin/bash
#$ -S /bin/bash
#$ -N D3.C.GATK  ##job name
#$ -l tmem=32G #RAM
#$ -l h_vmem=32G #enforced limit on shell memory usage
#$ -l h_rt=50:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-39

#Run on working directory
cd $SGE_O_WORKDIR || exit


#Path to software
export PATH=/share/apps/genomics/gatk-4.2.4.1:$PATH
export PATH=/share/apps/genomics/samtools-1.14/bin:$PATH
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH

# Set paths
MASTER=/SAN/ugi/LepGenomics/D3_Pararge_aegeria
bam_in=$MASTER/04_ATLAS/MODC_splitmerge
gvcfs=$MASTER/GATK_SNPcalling/MODC
REF=$MASTER/RefGenome2/GCF_905163445.1_ilParAegt1.1_genomic.fna


#ARRAY
NAME=$(sed "${SGE_TASK_ID}q;d" MODC_bamlist.forSNPs)

## Run haplotype caller     
time gatk --java-options "-Xmx10g" HaplotypeCaller \
-R $REF \
-I $bam_in/${NAME}.bam \
-O $gvcfs/${NAME}_g.vcf.gz -ERC GVCF
