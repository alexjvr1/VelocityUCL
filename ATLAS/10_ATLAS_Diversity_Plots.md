
# Theta estimates from ATLAS


See paper [here](https://www.genetics.org/content/genetics/205/1/317.full.pdf)

We can estimate a global theta, or estimate theta in windows. According to the paper, we can obtain accurate estimates in for data with <1X average depth in windows of 1Mb. 

The [ibis paper](https://www.sciencedirect.com/science/article/pii/S0960982218316099#mmc2) uses 5Mb windows. 



## 1. Estimate theta in windows of 0.5Mb 

These are smaller than the min window size suggested by the ATLAS paper. So we'll increase these window sizes for the final analysis. 

1. Create a bed file specifying the windows needed: 

```
/SAN/ugi/LepGenomics/E3_SubsetTests/02a_mapped_museum_FORANGSD

cat Regions.LR75_500kbwindow.bed 
LR761675.1	0	500000
LR761675.1	500001	1000000
LR761675.1	1000001	1500000
LR761675.1	1500001	2000000
LR761675.1	2000001	2500000
LR761675.1	2500001	3000000
LR761675.1	3000001	3500000
LR761675.1	3500001	4000000
LR761675.1	4000001	4500000
LR761675.1	4500001	5000000
LR761675.1	5000001	5500000
LR761675.1	5500001	6196582

```

ATLAS doesn't allow overlapping windows, so if we want a sliding window approach we'll need to run the analysis multiple times with shifted non-overlapping windows. 

Run ATLAS for MODE, MODC, and MUS either by submitting the following script (modified for each pop), or on the interactive server. 

Make sure to use the merged bam files (created with SplitMerge in ATLAS). You'll get errors about the read lengths exceeding the insert size if this hasn't been done. 

```
#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.ATLAS.GlobDiv  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)


#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
FOLDER=E3_SubsetTests
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern/modc.merged.bamlist
OUTPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern
TASK=theta


#Run analysis
while IFS= read -r line; do $ATLAS task=theta bam=$line window=Regions.LR75_500kbwindow.bed > $line.theta.500kbwindow; done < $INPUT

```


Download the resulting theta.txt.gz files to your home computer: 

```
/Users/alexjvr/2021postdoc/Velocity/E3_A.hyperantus/ATLAS

gunzip *gz
```

Add a column to each file specifying pop and indiv

```



