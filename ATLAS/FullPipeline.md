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
```


This will output the sequence names which contain the sequencing machine and lane number information. We can use this to get all the readIDs for reads associated with the different read groups. 
```

```


Now split the bam file by these reads using [FilterSamReads](https://broadinstitute.github.io/picard/command-line-overview.html#FilterSamReads)
```
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar


```
