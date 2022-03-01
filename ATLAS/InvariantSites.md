# Invaraint sites for Recal

Goal 1: Are the invariant sites really invariant? 

Sam extracted a set of possibly invariant sites identified from a multi-species alignment. I need to check that they actually are invariant before I can use them. 

Goal 2: Do we have enough data. According to [ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Sequence%20Data%20Processing%20Tools%3A%20recal) we need at least 1Mbp of invariants sites at 0.5X


## Data

- a bed file containing all the possible invariant sites

- bam files for all the individuals. 



## Strategy: Goal 1 - Are the invariant sites invariant?

We'll intersect these with the modern bam files. 

```
qrsh -l tmem=8G, h_vmem=8G, h_rt=3600

cd /SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/

```

Extract the sites from the bed file: 
```

```


We can use bcftools mpileup to quickly find the variation at these sites across all our bam files: 

```
bcftools=
```



## Strategy: Goal 2 - Do we have enough data for recal? 






