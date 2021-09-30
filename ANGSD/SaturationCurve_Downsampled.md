# Saturation of nuc div

Plot the saturation of estimates of nuc div from downsampled data


/SAN/ugi/LepGenomics/E3_SubsetTests

```
#Call software
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
FOLDER=E3_SubsetTests
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern_exp
OUTPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern_exp_Downsampled
TAIL="fullLR75.subset.bam"

paste mode.names.100 mode.prop.downsample mode.propnames | while IFS=$'\t' read -r var1 var2 var3; do java -jar $PICARD DownsampleSam \
I=$INPUT/"$var1".$TAIL \
O=$OUTPUT/"$var1".fullLR75.Downsampled."$var3".bam \
P="$var2"; done


##Where the following input files are used: 
# mode.names.100 = 100 times the mode.names 
#seq 1 10 | xarg -Inone cat mode.names >> mode.names.100
# mode.prop.downsample and mode.propnames are the proportions of sequences to subsample from each bam file, and what coverage that represents
#this info is in the excel sheet: VelocitySampleStats_2021_AJvR.xlsx

```

We now have 10 sets of downsampled bam files for each of MODE and MODC. We need to generate the folded SFS files for each of these populations: 

```


```


## minDP 1 vs minDP 2

Compare nucDiv estimated when allowing 1x and 2x minDP for each population
```
realSFS=/share/apps/genomics/angsd-0.935/bin/realSFS
thetaStat=/share/apps/genomics/angsd-0.935/bin/thetaStat


$realSFS MODC.fullLR75.minDP1.saf.idx -fold 1 > MODC.0.7X.minDP1.fold.SFS    ##1128156 sites
$realSFS MODC.fullLR75.minDP2.saf.idx -fold 1 > MODC.0.7X.minDP2.fold.SFS    ##38075 sites

$realSFS MODE.fullLR75.minDP1.saf.idx -fold 1 > MODE.0.7X.minDP1.fold.SFS    ##1164387 sites
$realSFS MODE.fullLR75.minDP2.saf.idx -fold 1 > MODE.0.7X.minDP2.fold.SFS    ##42161 sites

$realSFS MUS.fullLR75.minDP1.saf.idx -fold 1 > MUS.0.7X.minDP1.fold.SFS      ##726361 sites
$realSFS MUS.fullLR75.minDP2.saf.idx -fold 1 > MUS.0.7X.minDP2.fold.SFS      ##13947 sites
```


Estimate diversity for each of the SFS files
```
##MODC
$realSFS saf2theta MODC.fullLR75.minDP1.saf.idx -sfs MODC.0.7X.minDP1.fold.SFS -outname MODC.0.7X.minDP1
$thetaStat do_stat MODC.0.7X.minDP1.thetas.idx -win 50000 -step 10000 -outnames MODC.0.7X.minDP1.theta.window.gz

$realSFS saf2theta MODC.fullLR75.minDP2.saf.idx -sfs MODC.0.7X.minDP2.fold.SFS -outname MODC.0.7X.minDP2
$thetaStat do_stat MODC.0.7X.minDP2.thetas.idx -win 50000 -step 10000 -outnames MODC.0.7X.minDP2.theta.window.gz

##MODE
$realSFS saf2theta MODE.fullLR75.minDP1.saf.idx -sfs MODE.0.7X.minDP1.fold.SFS -outname MODE.0.7X.minDP1
$thetaStat do_stat MODE.0.7X.minDP1.thetas.idx -win 50000 -step 10000 -outnames MODE.0.7X.minDP1.theta.window.gz

$realSFS saf2theta MODE.fullLR75.minDP2.saf.idx -sfs MODE.0.7X.minDP2.fold.SFS -outname MODE.0.7X.minDP2
$thetaStat do_stat MODE.0.7X.minDP2.thetas.idx -win 50000 -step 10000 -outnames MODE.0.7X.minDP2.theta.window.gz


##MUS
$realSFS saf2theta MUS.fullLR75.minDP1.saf.idx -sfs MUS.0.7X.minDP1.fold.SFS -outname MUS.0.7X.minDP1
$thetaStat do_stat MUS.0.7X.minDP1.thetas.idx -win 50000 -step 10000 -outnames MUS.0.7X.minDP1.theta.window.gz

$realSFS saf2theta MUS.fullLR75.minDP2.saf.idx -sfs MUS.0.7X.minDP2.fold.SFS -outname MUS.0.7X.minDP2
$thetaStat do_stat MUS.0.7X.minDP2.thetas.idx -win 50000 -step 10000 -outnames MUS.0.7X.minDP2.theta.window.gz
```


copy to mac to plot: 
```
#in first window run 
ssh -l ajansen -L 3000:morecambe.cs.ucl.ac.uk:22 ajansen@tails.cs.ucl.ac.uk


#open a second window and navigate to the Velocity folder on your comp

/Users/alexjvr/2021postdoc/Velocity/E3_A.hyperantus/03.ANGSD
rsync -auve "ssh -p 3000" $i ajansen@localhost:/SAN/ugi/LepGenomics/E3_SubsetTests/03a_ANGSD/*PG .


```


It looks like diversity is underestimated in museum data at 2X already. This could be due to missing data, as the denominator in Watterson's theta is the harmonic mean of the number of samples for which there is data. However, in ANGSD this is the total samples in the dataset without taking into account missing data. Thus the denominator is too large for most sites in the museum dataset, which will result in a lower tW. 

We could use pixy to correct for this. 


![alt_txt][Fig1_hist]

[Fig1_hist]:https://user-images.githubusercontent.com/12142475/135264361-d2142190-6674-4e4e-9803-dd844c53741b.png



![alt_txt][Fig2_linegraph]

[Fig2_linegraph]:https://user-images.githubusercontent.com/12142475/135264347-20f3547a-b52b-4447-be51-cc3f4311cdba.png






# ATLAS


I can estimate the global theta in ATLAS from the bam files. 

## Modern data

First the modern bam files need to be processed with SplitMerge to merge any overlapping reads. 

Modify this scripts, where the input file is a list of bam file names: 

[04b_ATLAS_SplitMerge.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04b_ATLAS_SplitMerge.sh)



## Estimate theta

The global theta (here for LR75.1) can be estimated for each chromosome or region as specified in the bam file. 

Run the following command in the interactive node: 
```
#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
FOLDER=E3_SubsetTests
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_museum_FORANGSD/mus.bamlist
OUTPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_museum_FORANGSD/
TASK=theta

#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas

while IFS=  read -r line; do $ATLAS task=$TASK thetaGenomeWide bam=$line regions=Regions.LR75.bed > $line.theta; done < $INPUT


## Or run it as a script, but add in -l avx2=yes
while IFS=  read -r line; do $ATLAS task=$TASK thetaGenomeWide bam=$line regions=Regions.LR75.bed > $line.theta -l avx2=yes; done < $INPUT
```



Where the bed file is just a list of regions for which to estimate the globalTheta. Here chromosome LR76..75
```
cat Regions.LR75.bed 

LR761675.1	0	6196582

```

The output writes all the chromosomes (0 lines read, etc), but we can extract the data specifically for our chromsomes: 
```
for i in $(ls *theta); do grep -A 7 "LR761675.1" > $i.window.LR75.theta; done
```



These outputs include information on each 100kb window for 1) mean depth, the % of sites with 2x data, and % sites with no data. 

We can use this to plot the depth and data across the chromosome per individual. 

```


```
