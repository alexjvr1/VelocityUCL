# Diversity estimates

We're estimating genetic diversity using the new MODC dataset. I've also checked that MapDamage is working properly on the MUS samples.

I'll test

1: The effect of GL 1 vs GL 2. GL 2 (GATK model) performs better at very low coverage. 

See Section 4 [here](https://onlinelibrary.wiley.com/doi/10.1111/mec.16077)

```
"The two GL models differ in performance because both the gatk model and our simulation model assume that each base quality score reflects an independent and unbiased measurement of the probability of sequencing error (Huang et al., 2012; McKenna et al., 2010), whereas the samtools model assumes that if one sequencing error occurs at a certain locus, subsequent errors are more likely (Li, 2011; Li et al., 2009). As a result, with the samtools model, lower frequency mutations are less likely to be identified as polymorphic sites and more likely to be interpreted as sequencing errors when the coverage is low. This leads to an underestimation of the number of singleton mutations when using the samtools model, and therefore Watterson's θ tends to be underestimated, at least for our simulated data. We note, however, that these low-frequency SNPs have minimal impact on many other types of population genomic analyses and, in fact, are often filtered out. Consistent with this, we did not observe any strong discrepancies between the two GL models in other types of analysis that we performed in this study (Figures S4-S7). We also stress that the sequencing errors modelled in our simulations may not accurately represent the sequencing error profile in real life, so our result should not be interpreted as a recommendation of one GL model over the other."
```

2. Test downsampled vs full 



## Data Pipeline

### Museum

1. Cutadapt

2. bbRepair

3. bbMerge

4. Map to genome with bwa mem

5. MapDamage (on merged reads only) - estimate PMD and recalibrate base quality scores

6. markDuplicates, addRG and Realign


### Modern

1. Cutadapt

2. Map to genome with bwa mem

3. Add RG, mark Duplicates, and local realignment 

4. SplitMerge with ATLAS to merge overlapping reads. 



### Data quality

Check the bam files to estimate global depth, individual depth distribution, and distribution of base qualities. 


MODERN & MUS
```
#Run in interactive node
qrsh -l tmem=32G, h_vmem=32G

#Define angsd
angsd=/share/apps/genomics/angsd-0.935/bin/angsd

#Create two input files that list the bam file names
ls 02a_mapped_modern/*realn_mergedReads.bam >> ANGDS_mod.names
ls 02a_mapped_modern_exp/*realn_mergedReads.bam >> ANGDS_mode.names
ls 02a_mapped_museum_FORANGSD/*realn.bam >> ANGSD_mus.names

$angsd -b ANGDS_mod.names -ref RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -out 03.1_ANGSD_2021/MODC.qc -r LR761675.1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500 &> /dev/null

$angsd -b ANGDS_mode.names -ref RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -out 03.1_ANGSD_2021/MODE.qc -r LR761675.1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 -minMapQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500 &> /dev/null

$angsd -b ANGSD_mus.names -ref RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna -out 03.1_ANGSD_2021/MUS.qc -r LR761675.1 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 0 -trim 0 -C 50 -baq 1 -minMapQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 500 &> /dev/null

##Use NGStools script (plotQC.R) to calculate distributions
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_2021

export PATH=/share/apps/R-4.0.3/bin:$PATH

Rscript plotQC.R MODC.qc 2> /dev/null
Rscript plotQC.R MODE.qc 2> /dev/null
Rscript plotQC.R MUS.qc 2> /dev/null
```

MODE


![alt_txt][MODE.1]

[MODE.1]:https://user-images.githubusercontent.com/12142475/134003399-b9f6b15e-8c35-47fe-9062-6f05c73c1468.png


MODC

![alt_txt][MODC.1]

[MODC.1]:https://user-images.githubusercontent.com/12142475/134003414-733d2649-67c1-4aea-9288-7a0b0dae29c6.png



MUS

![alt_txt][MUS.1]

[MUS.1]:https://user-images.githubusercontent.com/12142475/134025020-5240366a-5245-48f7-93c6-41774325d27d.png

![alt_txt][MUS.2]

[MUS.2]:https://user-images.githubusercontent.com/12142475/134025047-d5700ae4-6464-4863-a792-5dbc72bdec63.png



### GL 1 - Full Data


#### Estimate SAF

```
SPECIESDIR=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus
POP=mod
#POP=mode
#POP=mus
REGION=LR761675.1
PP=1
##PP=0 for museum samples because sequences are merged into single reads. 
#PP=0
##Depth filters obtained from figures above
minDP=17
maxDP=350
OUTDIR=$SPECIESDIR/03.1_ANGSD/03.1_SAF
C=50
MININD=10
GL=1

time $angsd -b ANGDS_$POP.names -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs $PP -r $REGION
-GL $GL -doSaf 1 -anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna -doCounts 1 -setMinDepth $minDP -setMaxDepth $maxDP -doMajorMinor 4 -out /newhome/aj18951/E3_Aphantopus_hyperantus_2020/03.1_ANGSD_2021/03.1_SAF/$POP.$REGION.SEPT20 -C $C -baq 1 -dumpCounts 2 -doDepth 1 -doGlf 2 -minInd $MININD
```


MODC
```
-> Total number of sites analyzed: 5780437  (previous 5715589)
	-> Number of sites retained after filtering: 5377972  (4964197)
	[ALL done] cpu-time used =  1420.85 sec
	[ALL done] walltime used =  1425.00 sec
real	23m45.665s
user	23m32.493s
sys	0m8.385s
```


MODE
```
-> Total number of sites analyzed: 5688313 (5688111)
	-> Number of sites retained after filtering: 5227896  (5085827)
	[ALL done] cpu-time used =  1172.12 sec
	[ALL done] walltime used =  1176.00 sec

real	19m35.369s
user	19m20.978s
sys	0m11.173s
```

MUS   - we have an order of magnitude more SNPs! This is most likely due to the change in mapDamage (i.e. explicitly using merged reads only), and perhaps to do with the depth filters that were quite strict in the first run. 
```
	-> Total number of sites analyzed: 4784988  (4786808)
	-> Number of sites retained after filtering: 2978827 (276245)
	[ALL done] cpu-time used =  476.66 sec
	[ALL done] walltime used =  479.00 sec

real	7m59.237s
user	7m54.824s
sys	0m1.867s
```


#### Estimate SFS

```
realSFS=/share/apps/genomics/angsd-0.935/bin/realSFS 

cd /SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/03.1_ANGSD_2021/03.1_SAF/

$realSFS mode.LR761675.1.1.SEPT20.saf.idx mod.LR761675.1.SEPT20.saf.idx -fold 1 > MODE.MODC.fold.sfs

MODE-MODC:
Will run optimization on nSites: 5170774

MODE-MUS
Will run optimization on nSites: 2945538
```




### GL 2 - Full Data


```
SPECIESDIR=/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus
POP=mod
#POP=mode
#POP=mus
REGION=LR761675.1
PP=1
##PP=0 for museum samples because sequences are merged into single reads. 
#PP=0
##Depth filters obtained from figures above
minDP=17
maxDP=350
OUTDIR=$SPECIESDIR/03.1_ANGSD/03.1_SAF.FULL.GL2
C=50
MININD=10
GL=2

time $angsd -b ANGDS_$POP.names -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs $PP -r $REGION
-GL $GL -doSaf 1 -anc $SPECIESDIR/RefGenome/*fna -ref $SPECIESDIR/RefGenome/*fna -doCounts 1 -setMinDepth $minDP -setMaxDepth $maxDP -doMajorMinor 4 -out /newhome/aj18951/E3_Aphantopus_hyperantus_2020/03.1_ANGSD_2021/03.1_SAF/$POP.$REGION.$GL.SEPT20 -C $C -baq 1 -dumpCounts 2 -doDepth 1 -doGlf 2 -minInd $MININD

```


### GL 1 - Downsampled

```


```


### GL 2 - Downsampled

```

```


