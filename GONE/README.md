# Ne history

Based on our estimates of genetic variance over time, it looks like all species are declining. We need to independently validate this using by estimating the demographic history. 

[GONE](https://github.com/esrud/GONE/blob/master/USER´S%20GUIDE.pdf) was developed specifically to reconstruct Ne over time over the most recent 100-200 generations. 

They use LD at different distances between SNPs to estimate drift over time as a proxy for Ne. 

We will use the SNPs called using ATLAS

```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/04_ATLAS/GLs/*vcf.gz
```

ATLAS outputs information for every locus, thus we have very "gappy" data with missingness across loci and individuals. 

GONE can work with a lot of missingness and requires that no MAF filter is applied. 

However, a quality filter needs to be applied. The following reduces the dataset from ~400Mil SNPs to 
```
bcftools filter -O z -o MODC_filtered_Qual20.vcf.gz -i '%QUAL>20 ' MODC_ATLAS_majorMinor_majorMinor.vcf.gz

vcftools --gzvcf MODC_filtered_Qual20.vcf.gz 

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf MODC_filtered_Qual20.vcf.gz

Using zlib version: 1.2.7
After filtering, kept 38 out of 38 Individuals
After filtering, kept 12454911 out of a possible 12454911 Sites
Run Time = 94.00 seconds


##Or if we want to filter any genotypes where individuals have a GQ<20 and DP<8: 
bcftools filter -O z -o MODC_filtered_Qual20_DP8.vcf.gz -i 'FMT/GQ>20 & FMT/DP>7' MODC_filtered_Qual20_Autosomes.vcf.gz

vcftools --gzvcf MODC_filtered_Qual20_DP8.vcf.gz

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf MODC_filtered_Qual20_DP8.vcf.gz

Using zlib version: 1.2.7
After filtering, kept 38 out of 38 Individuals
After filtering, kept 10671495 out of a possible 10671495 Sites
Run Time = 61.00 seconds

## Delete the Z chromosome. For E3 this is LR761650.1. We're also excluding all the non-chromosome contigs/scaffolds. 

vcftools --gzvcf MODC_filtered_Qual20_DP8.vcf.gz --chr LR761647.1 --chr LR761648.1 --chr LR761649.1 --chr LR761651.1 --chr LR761652.1 --chr LR761653.1 --chr LR761654.1 --chr LR761655.1 --chr LR761656.1 --chr LR761657.1 --chr LR761658.1 --chr LR761659.1 --chr LR761660.1 --chr LR761661.1 --chr LR761662.1 --chr LR761663.1 --chr LR761664.1 --chr LR761665.1 --chr LR761666.1 --chr LR761667.1 --chr LR761668.1 --chr LR761669.1 --chr LR761670.1 --chr LR761671.1 --chr LR761672.1 --chr LR761673.1 --chr LR761674.1 --chr LR761675.1 --recode --recode-INFO-all --stdout | gzip -c > MODC_filtered_Qual20_DP8_Autosomes.vcf.gz



## The final dataset needs to be in plink format
## Create a chromosome map file
grep LR ../../RefGenome/*fai | awk '{print $1, $1, sep="\t"}' > E3.chrom-map.txt


vcftools --gzvcf MODC_filtered_Qual20_Autosomes.vcf.gz --recode --plink --chrom-map E3.chrom-map.txt --out MODC_filtered_forGONE

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--gzvcf MODC_filtered_Qual20_Autosomes.vcf.gz
	--chrom-map E3.chrom-map.txt
	--out MODC_filtered_forGONE
	--plink

Using zlib version: 1.2.7
After filtering, kept 38 out of 38 Individuals
Writing PLINK PED and MAP files ... 
	Read 29 chromosome mapping file entries.
Done.
After filtering, kept 12060531 out of a possible 12060531 Sites
Run Time = 409.00 seconds

```

The maximum number of SNPs allowed is 10Mil, and 1Mil per chromosome. We have about 2Mil too many. We can randomly select 10Mil SNPs in plink: 






Too see how many we have per chr: 
```
for i in $(cat Chr.map); do echo $i && grep $i MODC_filtered_forGONE.map |wc -l; done

LR761647.1
532249
LR761648.1
529196
LR761649.1
489314
LR761651.1
478192
LR761652.1
494363
LR761653.1
506271
LR761654.1
499994
LR761655.1
481173
LR761656.1
471907
LR761657.1
496805
LR761658.1
478793
LR761659.1
445319
LR761660.1
445054
LR761661.1
427766
LR761662.1
442035
LR761663.1
456870
LR761664.1
445389
LR761665.1
423448
LR761666.1
438023
LR761667.1
428606
LR761668.1
418423
LR761669.1
404461
LR761670.1
357258
LR761671.1
348420
LR761672.1
322363
LR761673.1
277181
LR761674.1
281148
LR761675.1
240510

```



## Input files: map and ped

We need plink ped and map files

The map file needs chromosomes named 1-xx, so we need to rename our LR chromosomes. 
```
#create a file mapping chr names to numbers
head Chr.map
LR761647.1 1 	
LR761648.1 2 	
LR761649.1 3 	
LR761651.1 4 	
LR761652.1 5 	
LR761653.1 6 	
LR761654.1 7 	
LR761655.1 8 	
LR761656.1 9 	
LR761657.1 10 
...


##Rename
awk '

NR==FNR {

a[$1]=$2

next

}

{

for(i=1;i<=NF;i++) 

if($i in a)

$i=a[$i]

}1

' FS=" " Chr.map FS="\t" MODC_filtered_forGONE.map > MODC_filtered_forGONE_renamed1.map
```

The snps need to be renamed to remove the colon
```
sed -i 's/:/_/' MODC_filtered_forGONE_renamed1.map

```

The ped file needs "-9" in column 6 (the column just before the genotypes). All preceding columns are ignored. 
```
awk '$6=="-9"' MODC_filtered_forGONE.ped
```


## Input files: [INPUT_PARAMETERS_FILE](https://github.com/alexjvr1/VelocityUCL/edit/main/GONE/INPUT_PARAMETERS_FILE)

Choosing a recombination rate: 

We don't have species specific estimates of recombination for any of the species, but we have estimates for other Lepidoptera. 

1. 2.97 - 4.0 cM / Mb in the silkmoth (Yamamoto et al., 2008; Yasukochi, 1998) 

2. 5.5 - 6.0 cM / Mb in different Heliconius species (Jiggins et al., 2005; Tobler et al., 2005).

3. 3.81-4.05 cM / Mb in [Vanessa cardui](https://www.biorxiv.org/content/10.1101/2022.04.14.488360v1.full.pdf)

We'll use the average: ~4.5 cM/Mb


