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
```



## Input files: map and ped

We need plink ped and map files



## Input files: [INPUT_PARAMETERS_FILE](https://github.com/alexjvr1/VelocityUCL/edit/main/GONE/INPUT_PARAMETERS_FILE)




