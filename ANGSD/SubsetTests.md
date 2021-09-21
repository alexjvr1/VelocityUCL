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
qrsh

/share/apps/genomics/samtools-1.9/bin/samtools
while read NAME <&4; do $samtools view $NAME.$TAIL "LR761675.1:1-1000000" -o $NAME.LR75.subset.bam; done 4<subset.list
```


## Test MapDamage is working on merged reads


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
