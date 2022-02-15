# ATLAS pipeline

We've been trying to piece together a pipeline which addresses all of the different sources of error in our data. But ATLAS has been specifically designed for low coverage ancient DNA samples. 

So I will rerun some samples using only the ATLAS pipeline to estimate theta. 

## Issues: 

1) The error rates differ between Illumina runs, so we need to model the error rate for each RG separately. However, I concatenated the raw data together right at the start before mapping. 

The flowcell and lane number information can be extracted from the fastq [file headers](https://help.basespace.illumina.com/articles/descriptive/fastq-files/). Genewiz follows the recommended headers exactly. 

For the species where I have already mapped and processed the bam reads (C3, E3, D3), we can extract this information from the bam file. 

First check what the read information is and what the uniqe read headers information is for all the sequences: 
```
samtools view bamfile | grep "at" |awk '{print $1, $2, $3, $4}'|sort |uniq
```

e.g. for AH-02-2019-42:
```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/02a_mapped_MODE.unmerged
samtools=/share/apps/genomics/samtools-1.14/bin/samtools

$samtools view AH-02-2019-42.realn.bam | grep "@" | awk '{print $1, $2, $3, $4}'|sort |uniq

A00410:144:HHW5LDRXX:2
A00410:136:HH7T3DRXX:1
A00410:136:HH7T3DRXX:2
```

And write these to a text file: 
```
nano E3.MODE_RG.list
ReadID	RGinbam	RGname	
A00410:144:HHW5LDRXX:2	HHW5LDRXX_lane2	RG1
A00410:136:HH7T3DRXX:1 HH7T3DRXX_lane1 RG2
A00410:136:HH7T3DRXX:2	HH7T3DRXX_lane2 RG3

```


Now split the bam file by these reads using [FilterSamReads](https://broadinstitute.github.io/picard/command-line-overview.html#FilterSamReads)
```
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar

java -jar $PICARD FilterSamReads I=AH-02-2019-43.realn.bam O=AH-02-2019-43.realn.RG1.bam READ_LIST_FILE=AH53.RG1.reads FILTER=includeReadList

Elapsed time: 3.96 minutes for one sample
```

And add the appropriate read group to each of the split files. 
```
time java -jar $PICARD AddOrReplaceReadGroups I=AH-02-2019-43.realn.RG1.bam O=AH-02-2019-43.realn.RG1.rgadded.bam RGID=RG1 RGLB=modern04 RGPL=Illumina4000 RGPU=HHW5LDRXX2 RGSM=43

real	3m36.398s
user	3m34.347s
sys	0m2.717s
```



##NOTE: Sam uses the RefSeq reference genome (annotated) to generate the bed file. We need to check that the chromosome names are the same between our reference and the RefSeq reference, and rename if necessary. 

We need to rename them for: *A. hyperantus*. See the two versions of chromosome names [here](https://www.ncbi.nlm.nih.gov/assembly/GCF_902806685.1). 






Run ATLAS: 
```
ATLAS=/share/apps/genomics/atlas-0.9/atlas
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH

#create a RG file
nano RG.txt
RG1 paired


##SplitMerge
for i in $(ls *RG1*bam); do time $ATLAS task=splitMerge bam=$i updateQuality readGroupSettings=RG.txt; done
#10mins per sample
```


Recal
```
##recal
time $ATLAS task=recal bam=AH-02-2019-42.realn.RG1.rgadded_mergedReads.bam region=../M_hyperantus.phast.1.766.bed verbose
```

Recal uses a LOT of memory! Based on the ATLAS notes, they required 55Gb RAM for an average of 2X sequencing, using the human genome. Our biggest dataset is MODE, which has a mean sequencing depth of ~3.5X. This requires a lot more memory. But we can limit the memory usage by limiting the depth (minDepth=2, and maxDepth=4), the size of the windows (default = 1Mil, windowSize=500000). If all else fails, we can limit the number of windows read in from each chromosome (eg. limitWindows=4)

Testing on 2 samples, this seems to work. I'm going to test

1) 80Gb RAM with no window limit, and maxDP = 6
```
#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODE.ATLAS.recal  ##job name
#$ -l tmem=64G #RAM
#$ -l h_vmem=64G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)


#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.19/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS/MODE/recal.bamlist
OUTPUT=$SHAREDFOLDER/$SPECIES/04_ATLAS/MODE
TASK=recal


#Run analysis
while IFS=  read -r line; do time $ATLAS task=$TASK bam=$line region=../A.hyperantus_LR76only.bam minDepth=2  maxDepth=4 limitWindows=2 verbose > $line.recal -l avx2=yes; done < $INPUT

```




##Subsample to the same depth as the museum data (0.8X)




##Estimate theta


```

###We can add lots of filters: 
readUpToDepth
- Setting window size to 1000000. (parameter 'window')
   - Will read data up to depth 1000 and ignore additional bases. (parameter 'readUpToDepth')
   - Will filter out bases with quality outside the range [1, 93] (parameters 'minQual', 'maxQual')
   - Will print qualities truncated to [33, 126] (parameters 'minOutQual', 'maxOutQual')
   - Will filter out windows with a missing data fraction > 1.000000. (parameter 'maxMissing')
   - Will filter out windows with a fraction of 'N' in reference > 1.000000. (parameter 'maxRefN')
   - Chromosomes with no further specifications are assumed to be diploid (parameters 'ploidy' or 'haploid' to change ploidy).
   - Mates that are farther than 2000 apart will be considered orphans. (parameter 'acceptedDistance')
   - Will ignore orphaned reads and not write them to BAM (use 'keepOrphans' to keep them).
   - Will keep random read for all of overlapping positions
   - Will update quality scores of prefered bases to reflect information from overlapping bases.
