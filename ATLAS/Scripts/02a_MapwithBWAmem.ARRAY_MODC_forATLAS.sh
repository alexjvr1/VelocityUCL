#!/bin/bash
#$ -S /bin/bash
#$ -N C3.MODC_BWAmem  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on memory shell usage
#$ -l h_rt=10:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-38

#run job in working directory
cd $SGE_O_WORKDIR 


##Software
BWA=/share/apps/genomics/bwa-0.7.17/bwa
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH

#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/01a_MODC_cutadapt_reads
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_modern

#Create a list of input files
ls $INPUT/*R1.fastq.gz | awk -F "/" '{print $NF}' > modc.R1.tomap
ls $INPUT/*R2.fastq.gz | awk -F "/" '{print $NF}' > modc.R2.tomap

echo "list of samples to map = modc.R1.tomap & modc.R2.tomap" >> map_modc.log

##Define ARRAY names
NAME1=$(sed "${SGE_TASK_ID}q;d" modc.R1.tomap)
NAME2=$(sed "${SGE_TASK_ID}q;d" modc.R2.tomap)


##Run mapping
echo "mapping started" >> map_modc.log
echo "---------------" >> map_modc.log

##Check if Ref Genome is indexed by bwa
if [[ ! $REF.fai ]]
then 
	echo $REF" not indexed. Indexing now" >> map_modc.log
	$BWA index $REF
else
	echo $REF" indexed" >> map_modc.log
fi


##Map using array

sample_name=`echo ${NAME1} | awk -F "." '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time $BWA mem $REF $INPUT/${NAME1} $INPUT/${NAME2}| samtools sort -o  $OUTPUT/${NAME1}.forATLAS.bam" >> map_mus.log
time $BWA mem $REF $INPUT/${NAME1} $INPUT/${NAME2} | samtools sort -o  $OUTPUT/${NAME1}.forATLAS.bam
