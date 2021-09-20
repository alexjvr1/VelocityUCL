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





### GL 1 - Full Data


### GL 2 - Full Data


### GL 1 - Downsampled


### GL 2 - Downsampled




