#!/bin/bash
#$ -S /bin/bash
#$ -N D3.MODC.findRG  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on memory shell usage
#$ -l h_rt=1:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)

#run job in working directory
cd $SGE_O_WORKDIR 

##Software
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH

#Variables
#CHR = the start of the chromosome names. This can be found in the reference index file. 
CHR="NC_"

#Script
for i in $(ls *mergedReads.bam); do samtools view $i | grep $CHR | awk -F ":" '{print $1, $2, $3, $4}' |sort |uniq >> ALL_RGs.txt; done
