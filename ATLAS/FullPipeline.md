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

##recal
time $ATLAS task=recal bam=AH-02-2019-42.realn.RG1.rgadded_mergedReads.bam region=../M_hyperantus.phast.1.766.bed verbose


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
