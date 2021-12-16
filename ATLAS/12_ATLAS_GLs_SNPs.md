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

