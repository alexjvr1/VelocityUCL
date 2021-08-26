#$ -S /bin/bash ##shell to be used
#$ -N E3.ATLAS.mus  ##job name
#$ -l tmem=16G #amount of memory you whish to request
#$ -l h_vmem=16G #enforced amount of shell storage
#$ -l h_rt=1:00:00 ##wall time
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

##InputFiles:
#1. RG.txt that contains read group information: 
#  E3mus  paired
#2. *realn.bam files that have been checked with ValidateSamFile

#Run on working directory
cd $SGE_O_WORKDIR 

#Set variables
#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum
TAIL1=realn.bam
TAIL2=realn_mergedReads.bam
RGFILE=RG.txt
LENGTH=50


#Set paths
ATLAS=/share/apps/genomics/atlas-0.9/atlas
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH

#Set up ARRAY job
NAME=$(sed "${SGE_TASK_ID}q;d" mus.names)

#Run script
time $ATLAS task=splitMerge bam=$INPUT/${NAME}.$TAIL1 readGroupSettings=$RGFILE

time $ATLAS task=PMD bam=$INPUT/${NAME}.$TAIL2 fasta=$REF length=$LENGTH

