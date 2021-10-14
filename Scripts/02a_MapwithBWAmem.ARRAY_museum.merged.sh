#!/bin/bash
#$ -S /bin/bash
#$ -N C3.BWAmem_mod  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on memory shell usage
#$ -l h_rt=15:00:00 ##wall time.  
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
INPUT=$SHAREDFOLDER/$SPECIES/01c_musALL_merged
OUTPUT=$SHAREDFOLDER/$SPECIES/02a_mapped_museum_MERGED
TAIL="R1.repaired.fastq.gz.repaired.merged.fastq.gz"

##Define ARRAY names
NAME=$(sed "${SGE_TASK_ID}q;d" mus.names)


##Run mapping
echo "mapping started" >> map.log
echo "---------------" >> map.log

##Check if Ref Genome is indexed by bwa
if [[ ! $REF.fai ]]
then 
	echo $REF" not indexed. Indexing now"
	$BWA index $REF
else
	echo $REF" indexed"
fi


##Map using array

sample_name=`echo ${NAME1} | awk -F "." '{print $1}'`
echo "[mapping running for] $sample_name"
printf "\n"
echo "time $BWA mem $REF $INPUT/${NAME}.$TAIL| samtools sort -o  $OUTPUT/${NAME}.MERGED.bam" >> $OUTPUT/map_mus_MERGED.log
time $BWA mem $REF $INPUT/${NAME}.$TAIL | samtools sort -o  $OUTPUT/${NAME}.MERGED.bam
