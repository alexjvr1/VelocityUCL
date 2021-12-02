# Call GLs for final data


## 1. Subset bam and estimate depth

We're initially working with only one chromosome. We'll subset the bam file to extract only the smallest chromsome. And we'll calculate the depth per sample for these. 

Use the [subset_bam.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/susbet_bam.sh) script. 


Initial results allowing a lot of missingness: 

```
==> E3.SAF.C.o4267657 <==
	-> Wed Dec  1 16:33:49 2021
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 5775262
	-> Number of sites retained after filtering: 5389598 
	[ALL done] cpu-time used =  1289.71 sec
	[ALL done] walltime used =  1302.00 sec

real	21m42.038s
user	21m22.069s
sys	0m7.692s

==> E3.SAF.E.o4267655 <==
	-> Wed Dec  1 16:34:27 2021
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 5679903
	-> Number of sites retained after filtering: 5233482 
	[ALL done] cpu-time used =  1318.59 sec
	[ALL done] walltime used =  1334.00 sec

real	22m14.027s
user	21m40.970s
sys	0m17.656s

==> E3.SAF.mus.o4267654 <==
	-> Wed Dec  1 16:18:05 2021
	-> Arguments and parameters for all analysis are located in .arg file
	-> Total number of sites analyzed: 4820076
	-> Number of sites retained after filtering: 2268882 
	[ALL done] cpu-time used =  529.55 sec
	[ALL done] walltime used =  533.00 sec

real	8m53.602s
user	8m46.574s
sys	0m2.988s

```


But we need to filter for missingness

#### 1. 80% genotyping rate: 

Results
```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_Dec2021/03.1_SAF_0.8GenotypingRate

```


```
MUS

-> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 4820076
        -> Number of sites retained after filtering: 3432 
        
real    2m36.221s
user    2m27.710s
sys     0m2.244s


MODC

-> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 5775262
        -> Number of sites retained after filtering: 4469532 
        [ALL done] cpu-time used =  980.42 sec
        [ALL done] walltime used =  1051.00 sec


MODE

        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 5679903
        -> Number of sites retained after filtering: 4381798 
        [ALL done] cpu-time used =  1093.45 sec
        [ALL done] walltime used =  1120.00 sec

real    18m40.268s
user    17m56.608s
sys     0m16.881s


```

#### 2. 100% genotyping rate

Results
```
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_Dec2021/03.1_SAF_

```


```
MUS
        -> Arguments and parameters for all analysis are located in .arg file
        -> Total number of sites analyzed: 4820076
        -> Number of sites retained after filtering: 87 
        [ALL done] cpu-time used =  147.46 sec
        [ALL done] walltime used =  155.00 sec

real    2m34.426s
user    2m21.534s
sys     0m5.938s


```



## minDP

We also need to check if the depth makes a difference. Ideally we want minDP 5x, but this would mean that many museum individuals drop out. 

