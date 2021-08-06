# Variant calling

Using bcftools call on bluecrystal3

```
/newhome/aj18951/C3_Aricia_agestis_2020/03.2_Variant_calling_modern

```

Variants are called per chromosome
```
job0001.20210806-132508.report.txt  job0011.var_calling.log             job0021.variants.raw.bcf
job0001.var_calling.log             job0011.variants.raw.bcf            job0021.variants.raw.bcf.csi
job0001.variants.raw.bcf            job0011.variants.raw.bcf.csi        job0022.20210806-133223.report.txt
job0001.variants.raw.bcf.csi        job0012.20210806-132813.report.txt  job0022.var_calling.log
job0002.20210806-132508.report.txt  job0012.var_calling.log             job0022.variants.raw.bcf
job0002.var_calling.log             job0012.variants.raw.bcf            job0022.variants.raw.bcf.csi
job0002.variants.raw.bcf            job0012.variants.raw.bcf.csi        job0023.20210806-133224.report.txt
```

## Chr 

We're working with the smallest chromosome first to test the diploS/HIC analysis

1. The raw bcf files need to be processed to "see" missing data
```
module load apps/bcftools-1.8
bcftools filter -S . -O u -e 'FMT/DP=0' job0022.variants.raw.bcf |bcftools view -O b -o job0022.withmissing.bcf
```

2. Filters

We're applying some hard filters to keep only the best loci in the final dataset. 

a) mapping quality and sequencing quality of >20; b) depth of > 5x; c) max missingness of 0.5
```
vcftools --bcf 
```


