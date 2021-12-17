# GLs and Variants


We can use ATLAS to estimate genotype likelihoods (GLF). Subsequently we can find the two most likely alleles and write the likelihoods to a vcf file for each site and each individual (MajorMinor). 

We're running this for all three populations based on the merged bam files (MODC and MODE) and the mapdamage corrected files (MUS). 



```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/04_ATLAS/

```


Estimate GLF for each population using the [GLF.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/GLF.sh) script. This runs for 10-20 minutes per sample, run as an array (~30minutes total per population)

Next, find the MajorMinor alleles for each site for each population using the [MajorMinor.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/MajorMinor.sh) script.

This runs for ~10 hours per population, and creates really big vcf.gz files: 

```
du -sch */ATLAS*vcf.gz
17G	MODC/ATLAS_majorMinor_majorMinor.vcf.gz
17G	MODE/ATLAS_majorMinor_majorMinor.vcf.gz
8.3G	MUS/ATLAS_majorMinor_majorMinor.vcf.gz
42G	total
```

These are too big to read in, so we'll split by chromosome: 
```
/share/apps/genomics/htslib-1.14/bin/tabix 

```



Check vcf file: 

E3
```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/04_ATLAS/MODC

#MODC

$vcftools --gzvcf ATLAS_majorMinor_majorMinor.vcf.gz

VCFtools - v0.1.13
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf ATLAS_majorMinor_majorMinor.vcf.gz

Using zlib version: 1.2.7
After filtering, kept 38 out of 38 Individuals



```


### Smaller vcf

Create smaller vcf files based on 16Mb of data, minInd10 and minVariantQuality of 3. This is a likelihood ratio of Major/Major vs minor allele

```
$ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000 minVariantQual=3

```

D3
```
MODE
time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000 minVariantQual=3


time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000 minVariantQual=2


time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000


## Using minVariantQuality=3
vcftools=/share/apps/genomics/vcftools-0.1.16/bin/vcftools
$vcftools --gzvcf ATLAS_majorMinor_majorMinor.vcf.gz --exclude-positions remove.singletons --recode --recode-INFO-all --out ATLAS.nosingletons

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf ATLAS_majorMinor_majorMinor.vcf.gz
	--exclude-positions remove.singletons
	--recode-INFO-all
	--out ATLAS.nosingletons
	--recode

Using zlib version: 1.2.7
After filtering, kept 36 out of 36 Individuals
Outputting VCF file...
After filtering, kept 358447 out of a possible 1037042 Sites




MODC
##Test LRT

time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000 minVariantQual=3


time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000 minVariantQual=2


time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000





## Using minVariantQuality=3
- Atlas terminated successfully in 7.63333 min!

$vcftools --gzvcf ATLAS_majorMinor_majorMinor.vcf.gz --singletons

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf ATLAS_majorMinor_majorMinor.vcf.gz
	--singletons

Using zlib version: 1.2.7
After filtering, kept 33 out of 33 Individuals
Outputting Singleton Locations
After filtering, kept 864331 out of a possible 864331 Sites
Run Time = 9.00 seconds

#write first two columns and use to remove singletons and private doubletons

awk '{print $1, $2}' out.singletons > toremove.singletons

$vcftools --gzvcf ATLAS_majorMinor_majorMinor.vcf.gz --exclude-positions toremove.singletons --recode --recode-INFO-all --out MODC.nosingletons

After filtering, kept 33 out of 33 Individuals
Outputting VCF file...
After filtering, kept 302824 out of a possible 864331 Sites


```


MUS
```
time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000 minVariantQual=3
- Atlas terminated successfully in 4.41667 min!


time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000 minVariantQual=2


time $ATLAS task=$TASK glf=$INPUT/$GLFLIST minSamplesWithData=10 limitSites=16000000


$vcftools --gzvcf ATLAS_majorMinor_majorMinor.vcf.gz --freq2
After filtering, kept 1010512 out of a possible 1010512 Sites

```
