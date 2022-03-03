# Invariant sites for Recal

Goal 1: Are the invariant sites really invariant? 

Sam extracted a set of possibly invariant sites identified from a multi-species alignment. I need to check that they actually are invariant before I can use them. 

Goal 2: Do we have enough data. According to [ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Sequence%20Data%20Processing%20Tools%3A%20recal) we need at least 1Mbp of invariants sites at 0.5X


## Data

- a bed file containing all the possible invariant sites

- bam files for all the individuals. 



## Strategy: Goal 1 - Are the invariant sites invariant?

We'll intersect the invariant sites in the bed file with the modern bam files. 

```
qrsh -l tmem=8G, h_vmem=8G, h_rt=3600

cd /SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/RG1bams

```


1. Extract all the reads that cover the invariant sites [ExtractInvariantSitesFromBam.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/ExtractInvariantSitesFromBam.sh) from the MODE.unmerged and MODC.unmerged realn.bam datasets. 


2. Concatenate all the bam files together. 


Other approaches: 

This was run for the museum data RG1 dataset. 
```
bedtools=/share/apps/genomics/bedtools-2.30.0/bin/bedtools

$bedtools intersect -abam RG1.merged.bam -b A.hyperantus_LR76only.bed > RG1.invariant.merged.bam
```



OR

Extract the sites from the bed file using [Subsetbybed.sh]() 
```
java -Xmx6g -Xms6g -jar $PICARD FilterSamReads I=${NAME} O=${NAME}.RG1.bam READ_LIST_FILE=${NAME}.RG1.list FILTER=includeReadList

```




We can use bcftools mpileup to quickly find the variation at these sites across all our bam files: 

```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/02a_mapped_museum/RG1bams

bcftools=/share/apps/genomics/bcftools-1.14/bin/bcftools
$bcftools mpileup --fasta-ref ../../RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna RG1.merged.bam -a INFO/AD -R A.hyperantus_LR76only.bed > RG1.merged.mpileup 
```

Check the allelic depth (AD). We expect that these sites are invariant, but the variance that is present tells us what the error rate is within the Illumina run. 
```
grep "AD=" RG1.invariant.merged.mpileup |awk -F ";" '{print $2}'
```

We'll check what the range of AD is and choose a cut-off for the alt AD before excluding loci as invariant sites. 




## Strategy: Goal 2 - Do we have enough data for recal? 


Merge the bam files for each RG independently: 

```
samtools=/share/apps/genomics/samtools-1.14/bin/samtools

#-r adds the sample name to each RGID. This way we can keep the individual information. 
#When running ATLAS on this file, use poolReadGroups

#list all the bam files to be merged
ls *bam > RG1.bamlist  

#and merge
$samtools merge -r -b RG1.bamlist RG1.merged.bam
```


What is the coverage for the pooled data per invariant site: 
```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/02a_mapped_museum/RG1bams

$samtools bedcov A.hyperantus_LR76only.bed RG1.merged.bam > RG1merged.bedcov.out
```





