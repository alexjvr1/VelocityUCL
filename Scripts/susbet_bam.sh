#$ -S /bin/bash ##shell to be used
#$ -N E3_MUS_LR75  ##job name
#$ -l tmem=16G #amount of memory you whish to request
#$ -l h_vmem=16G #enforced amount of shell storage
#$ -l h_rt=1:00:00 ##wall time
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

##InputFile:
#1. *.bam files that have been recalibrated with MapDamage

#Run on working directory
cd $SGE_O_WORKDIR 

#Set variables
#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/E3_Aphantopus_hyperantus/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02c_MapDamage_MUS
OUTPUT=$SHAREDFOLDER/$SPECIES/02c_MapDamage_MUS
TAIL=LR75.bam


#Set paths
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH

#Set up ARRAY job
NAME=$(sed "${SGE_TASK_ID}q;d" mus.names)


echo "time samtools view -b ${NAME}.bam LR761675.1 > ${NAME}.LR75.bam" &&\
time samtools view -b ${NAME}.bam LR761675.1 > ${NAME}.LR75.bam && \
echo "time samtools index ${NAME}.LR75.bam" && \
time samtools index ${NAME}.LR75.bam && \
echo "time samtools depth ${NAME}.LR75.bam > ${NAME}.LR75.depth" &&\
time samtools depth ${NAME}.LR75.bam > ${NAME}.LR75.depth
