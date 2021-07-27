#!/bin/bash
#$ -S /bin/bash
#$ -N C3.BWAmem_mod  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on memory shell usage
#$ -l h_rt=10:00:00 ##wall time.  
#$ -j y  #concatenates error and output files (with prefix job1)
#$ -t 1-48

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
INPUT=$SHAREDFOLDER/$SPECIES/01a_museum_cutadapt_reads
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum

#Create a list of input files
ls $INPUT/*R1.fastq.gz | awk -F "/" '{print $NF}' > museum.R1.tomap
ls $INPUT/*R2.fastq.gz | awk -F "/" '{print $NF}' > museum.R2.tomap

echo "list of samples to map = museum.R1.tomap & museum.R2.tomap" >> map_mus.log

##Define ARRAY names
NAME1=$(sed "${SGE_TASK_ID}q;d" museum.R1.tomap)
NAME2=$(sed "${SGE_TASK_ID}q;d" museum.R2.tomap)


##Run mapping
echo "mapping started" >> map.log
echo "---------------" >> map.log

##Check if Ref Genome is indexed by bwa
if [[ ! $REF.fai ]]
then 
	echo $REF" not indexed. Indexing now" >> map_mus.log
	$BWA index $REF
else
	echo $REF" indexed" >> map_mus.log
fi


##Map using array

sample_name=`echo ${NAME1} | awk -F "." '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time $BWA mem $REF $INPUT/${NAME1} $INPUT/${NAME2}| samtools sort -o  $OUTPUT/${NAME1}.forATLAS.bam" >> map_mus.log
time $BWA mem $REF $INPUT/${NAME1} $INPUT/${NAME2} | samtools sort -o  $OUTPUT/${NAME1}.forATLAS.bam
