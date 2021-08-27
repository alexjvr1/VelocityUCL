# Raw Data Check

I identified an issue with E3 MODC with duplicate reads found in the raw fastq files (names, sequence, and quality scores). This accounts for ~23% of the reads (based on one sample tested - AH-01-2017-40). 

Picard Tools MarkDuplicates failed to run on all E3 MODC samples due to the presence of duplicate read names. 

E3 MODC was sequenced in Modern 03 (library name), along with E1_Erebia_epiphron, E2_Erebia_aethiops, H2_Miltochrista_miniata, H3_Eilema_griseola, and I1_Xanthorhoe_fluctuata. 

Modern 03 01/05/2019	HiSeq 4000 GENEWIZ



# Check for duplicates in raw fastq files from E3 MODC


Spot check some individuals 
```
for i in $(ls *fastq.gz); do gunzip $i; done
for i in $(ls *fastq); do sort $i |uniq -cd > $i.duplicates
for i in $(ls *duplicates); do grep "@" |wc -l; done
```


## SeqKit

A program designed to remove duplicate reads can be found [here](https://github.com/shenwei356/seqkit/#installation).


We can use this to count the number of duplicate reads in each sample, and create cleaned fastq files with duplicates removed. 

I'm using it to quickly check all indivs in Modern3 and Modern4: 

On BlueCrystal. 

Copy all raw mod3 and mod4 data to the home directory from where they are in long-term storage. 


```
#!/bin/bash
#PBS -N AH.mode.dups ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=24gb #RAM
#PBS -l walltime=2:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)
#PBS -t 1-80

#run job in working directory
cd $PBS_O_WORKDIR 

#Define variables
SEQKIT=/newhome/aj18951/software/seqkit

#SET UP ARRAY
NAME=$(sed "${PBS_ARRAYID}q;d" AH.mode.names)


#RUN SEQKIT

$SEQKIT rmdup ${NAME} -i -n -o ${NAME}.clean.fastqc.gz -D ${NAME}.duplicates.txt
```

## Results for each library and species


### Modern3







