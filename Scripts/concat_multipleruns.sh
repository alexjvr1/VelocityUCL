#$ -S /bin/bash ##shell to be used
#$ -N E2.R1.concat  ##job name
#$ -l tmem=16G #amount of memory you whish to request
#$ -l h_vmem=16G #enforced amount of shell storage
#$ -l h_rt=3:00:00 ##wall time
#$ -j y  #concatenates error and output files (with prefix job1)

#Script to concatenate sample run across multiple lanes. 
#Assumes data is split into R1 and R2 folders. This script is written for R1 - Modify and rerun in R2 folder. 

#Run on working directory
cd $SGE_O_WORKDIR 


for indiv in $(ls *.gz | cut -f1 -d_ | sort -u)
do
    cat "$indiv"*.gz > "$indiv.concat.R1.fastq.gz"
done

