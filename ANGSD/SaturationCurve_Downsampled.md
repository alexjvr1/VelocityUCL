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
