# Invaraint sites for Recal

Goal: Are the invariant sites really invariant? 

Sam extracted a set of possibly invariant sites identified from a multi-species alignment. I need to check that they actually are invariant before I can use them. 

## Data

a bed file containing all the possible sites

bam files for all the individuals. 



## Strategy

We'll intersect these with the modern bam files. 

```
qrsh -l tmem=8G, h_vmem=8G, h_rt=3600

cd /SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/

```

We can use bcftools mpileup to quickly find the variation at these sites across all our bam files: 

```
bcftools=
```
