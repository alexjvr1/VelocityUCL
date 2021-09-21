# Testing diversity and depth estimates with 10 samples from each pop, 1Mb of LR767675.1


## Subset 1Mb from each pop

MUS
```
pwd
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/02a_mapped_museum_FORANGSD

cat subset.list 
AH-01-1900-04
AH-01-1900-06
AH-01-1900-10
AH-01-1900-13
AH-01-1900-14
AH-01-1900-29
AH-01-1900-36
AH-01-1900-37
AH-01-1900-38
AH-01-1900-41
```


MODC
```
pwd
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/02a_mapped_modern

##merged by ATLAS SplitMerge. We want to test if these merge samples work in ANGSD (as merged reads) or if they're ATLAS specific
TAIL=realn_mergedReads.bam

##before merge
TAIL2=realn.bam

cat subset.list
AH-01-2016-01
AH-01-2016-04
AH-01-2016-05
AH-01-2016-06
AH-01-2016-07
AH-01-2016-08
AH-01-2016-09
AH-01-2016-10
AH-01-2016-11
AH-01-2016-12

```


MODE
```
pwd

AH-02-2019-42
AH-02-2019-43
AH-02-2019-44
AH-02-2019-45
AH-02-2019-46
AH-02-2019-47
AH-02-2019-48
AH-02-2019-50
AH-02-2019-51
AH-02-2019-52
```

```
qrsh -l tmem=16G, h_vmem=16G

/share/apps/genomics/samtools-1.9/bin/samtools
while read NAME <&4; do $samtools view $NAME.$TAIL "LR761675.1:1-1000000" -o $NAME.LR75.subset.bam; done 4<subset.list
```

## MapDamage tests

### Test MapDamage is working on bbmerge reads (ANGSD pipeline)


```
/SAN/ugi/LepGenomics/E3_SubsetTests/02a_mapped_museum_FORANGSD

# Software
##python
export PATH=/share/apps/python-3.8.5-shared/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/python-3.8.5-shared/lib:$LD_LIBRARY_PATH

##R
export PATH=/share/apps/R-4.0.3/bin:$PATH

##mapDamage
mapDamage="/share/apps/python-3.8.5-shared/bin/mapDamage"

SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna


time $mapDamage --merge-libraries -i AH-01-1900-04.LR75.subset.bam -d MAPDAMAGE -r $REF --rescale --single-stranded 


```

This works. The samples have been merged with bbMerge before mapping to the reference. 


### Test MapDamage is working ATLAS merged bams 


We'll try this with the modern samples that were merged using SplitMerge from Atlas 
```

```


## ANGSD

### Estimate SAF for each population


MODC
```
/SAN/ugi/LepGenomics/E3_SubsetTests/02a_mapped_modern

cat MODC.10.poplist 
AH-01-2016-01.LR75.subset.bam
AH-01-2016-04.LR75.subset.bam
AH-01-2016-05.LR75.subset.bam
AH-01-2016-06.LR75.subset.bam
AH-01-2016-07.LR75.subset.bam
AH-01-2016-08.LR75.subset.bam
AH-01-2016-09.LR75.subset.bam
AH-01-2016-10.LR75.subset.bam
AH-01-2016-11.LR75.subset.bam
AH-01-2016-12.LR75.subset.bam


ANGSD="/share/apps/genomics/angsd-0.935/bin/angsd"

#Set filters
#N="38"
N="10"
MININD="5"
MINMAF=""
MINQ="20"
minMAPQ="20"
minDP="2"
maxDP="100"
POP="MODC"
C="50"
POPLIST="MODC.10.poplist"
SPECIESDIR="/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus"
OUTDIR="/SAN/ugi/LepGenomics/E3_SubsetTests"
PP=1 #use all reads. Flag 1 uses only proper pairs, but	MODC has	merged reads. NB to filter for proper pair reads in the bamfiles using samtools before this point


time $ANGSD -b MODC.10.poplist -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs $PP -GL 1 -doSaf 1 -anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna -doCounts 1 -setMinDepthInd $minDP -setMaxDepth $maxDP -doMajorMinor 4 -out $OUTDIR/$POP.1Msubset.July21 -C $C -baq 1 -dumpCounts 2 -doDepth 1 -doGlf 2 -minInd $MININD


	-> Tue Sep 21 15:15:51 2021
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 927370
	-> Number of sites retained after filtering: 782145 
	[ALL done] cpu-time used =  43.22 sec
	[ALL done] walltime used =  44.00 sec

real	0m43.597s
user	0m43.112s
sys	0m0.116s

```


MODE 
```
	-> Tue Sep 21 15:26:29 2021
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 921428
	-> Number of sites retained after filtering: 807561 
	[ALL done] cpu-time used =  57.30 sec
	[ALL done] walltime used =  58.00 sec

real	0m58.554s
user	0m57.001s
sys	0m0.316s


```




MUS
```
	-> Tue Sep 21 15:36:38 2021
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 671007
	-> Number of sites retained after filtering: 112500 
	[ALL done] cpu-time used =  7.97 sec
	[ALL done] walltime used =  10.00 sec

real	0m9.517s
user	0m7.835s
sys	0m0.148s


```
