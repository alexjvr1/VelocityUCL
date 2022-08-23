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


vcftools --gzvcf MODC_filtered_Qual20.vcf.gz --chr LR761647.1 --chr LR761648.1 --chr LR761649.1 --chr LR761651.1 --chr LR761652.1 --chr LR761653.1 --chr LR761654.1 --chr LR761655.1 --chr LR761656.1 --chr LR761657.1 --chr LR761658.1 --chr LR761659.1 --chr LR761660.1 --chr LR761661.1 --chr LR761662.1 --chr LR761663.1 --chr LR761664.1 --chr LR761665.1 --chr LR761666.1 --chr LR761667.1 --chr LR761668.1 --chr LR761669.1 --chr LR761670.1 --chr LR761671.1 --chr LR761672.1 --chr LR761673.1 --chr LR761674.1 --chr LR761675.1 --recode --recode-INFO-all --out MODC_filtered_Qual20_Autosomes 
 

```



## Input files: map and ped

We need plink ped and map files



## Input files: [INPUT_PARAMETERS_FILE](https://github.com/alexjvr1/VelocityUCL/edit/main/GONE/INPUT_PARAMETERS_FILE)

Choosing a recombination rate: 

We don't have species specific estimates of recombination for any of the species, but we have estimates for other Lepidoptera. 

1. 

2. 

3. [Vanessa cardui](https://www.biorxiv.org/content/10.1101/2022.04.14.488360v1.full.pdf): 


