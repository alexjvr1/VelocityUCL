# Plots of MapDamage for Modern and museum data

1. What level of correction is applied to the modern and museum datasets? 

2. Does C -> T and G -> A have the same proportion of correction in each dataset? 

3. How does MapDamage decide to correct for C->T over G->A? 


## Run MapDamage

Mapdamge is run for both modern and museum samples, with input files here: 

```
/SAN/ugi/LepGenomics/E3_SubsetTests/02a_mapped_museum_FORANGSD

/SAN/ugi/LepGenomics/E3_SubsetTests/02a_mapped_modern_Downsampled

/SAN/ugi/LepGenomics/E3_SubsetTests/02a_mapped_modern_exp_Downsampled
```

Using this script: [04a.0_MapDamage_FORANGSD.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/04a.0_MapDamage_FORANGSD.sh), modified for MODC and MODE as necessary. 



## Extract output information for plots


#### 1. How much correction is applied across the sequence?

Plot the original and corrected mean frequencies per bp. 


