
# Theta estimates from ATLAS


See paper [here](https://www.genetics.org/content/genetics/205/1/317.full.pdf)

We can estimate a global theta, or estimate theta in windows. According to the paper, we can obtain accurate estimates in for data with <1X average depth in windows of 1Mb. 

The [ibis paper](https://www.sciencedirect.com/science/article/pii/S0960982218316099#mmc2) uses 5Mb windows. 


## 1. Downsample modern data

Modern data was sequenced to a much higher coverage than the museum data. Downsample these data to a) determine how sequencing depth affects estimates of theta, b) choose a deth comparable to the MUS dataset. 

a) Create folders for downsampled data

```
mkdir 04_Downsampled_MODC
mkdir 04_Downsampled_MODE
```

b) Generate the three input files needed for the Downsample.sh script: 

```
#1) mode.names.4sets = a list of sample names, repeated for the series of subsampling (here to 4 different depths). 
seq 1 4 | xargs -Inone cat mode.names > mode.names.4sets

#2) mode.prop.downsample = the proportion of data to sample for each individual

0.2, 0.4, 0.6, 0.8

#3) mode.propnames = a name for the output. I've used the expected depth after downsampling

0.2X, 0.4X, 0.6X, 0.8X


```


c) Run the [02b.5_Downsample.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/02b.5_Downsample_MODE.sh) script for each population. This will write results to the Downsampled folders. 


## 2. Run SpitMerge for MODC and MODE

Make sure to use the merged bam files (created with SplitMerge in ATLAS). You'll get errors about the read lengths exceeding the insert size if this hasn't been done. 


a) Run SplitMerge from ATLAS for 
  
        1) The full Modern bam files (.realn.bam) 

        2) all downsampled populations using the [04b_ATLAS_SplitMerge.sh (https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/04b_ATLAS_SplitMerge.sh) script. 

Run this in the 02a_mapped_MODC.unmerged and 02a_mapped_MODE.unmerged folders. 

You'll need to 
        
        1) index all the bam files before running ALTAS. 
        
        ATLAS expects indexed files to have the extension bam.bai, so rename files accordingly:
        
```
/share/apps/perl-5.30.0/bin/perl5.30.0 /home/ajansen/prename.pl 's:realn.bai:realn.bam.bai:g' *bai
```

        2) Create an RG.txt file for the SplitMerge run. This file contains the RG and if the reads are PE or SE: 
        
```
        #e.g.
        cat RG.txt
        E3modc paired
```


b) Check bam files using [ValidateSamFiles.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/02b.3_ValidateSamFile_MODC.sh) and fix if necessary. 



## 3. Estimate theta in windows of 1Mb for MODE, MODC, and MUS

1Mb windows are suggested for low coverage data. 

The script runs for 3-5 hours per population. 

Run ATLAS for MODE, MODC, and MUS using the [Theta.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/Scripts/04b_ATLAS_Theta_Window.sh) script. 


If the analysis should be run on a specific regions or chromosomes, create a bed file to specify the windows. Skip this note when analysing the whole genome:
```
/SAN/ugi/LepGenomics/E3_SubsetTests/02a_mapped_museum_FORANGSD

cat Regions.LR75_500kbwindow.bed 
LR761675.1	0	500000
LR761675.1	500001	1000000
LR761675.1	1000001	1500000
LR761675.1	1500001	2000000
LR761675.1	2000001	2500000
LR761675.1	2500001	3000000
LR761675.1	3000001	3500000
LR761675.1	3500001	4000000
LR761675.1	4000001	4500000
LR761675.1	4500001	5000000
LR761675.1	5000001	5500000
LR761675.1	5500001	6196582

```

ATLAS doesn't allow overlapping windows, so if we want a sliding window approach we'll need to run the analysis multiple times with shifted non-overlapping windows. 



### 4. Estimate global theta for each chromosome  (Skip this step initially) 

Regions for which to estimate global theta can be specified using a bed file. 

We'll estimate a global (mean) theta for each chromosome and each individual. 

```
#!/bin/bash
#$ -S /bin/bash
#$ -N E3.MODC.ATLAS.GlobDiv  ##job name
#$ -l tmem=16G #RAM
#$ -l h_vmem=16G #enforced limit on shell memory usage
#$ -l h_rt=1:00:00 ##wall time.
#$ -j y  #concatenates error and output files (with prefix job1)


#Run on working directory
cd $SGE_O_WORKDIR 


#Path to software
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
ATLAS=/share/apps/genomics/atlas-0.9/atlas


#Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
FOLDER=E3_SubsetTests
REF=$SHAREDFOLDER/$SPECIES/RefGenome/GCA_902806685.1_iAphHyp1.1_genomic.fna
INPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern/modc.merged.bamlist
OUTPUT=$SHAREDFOLDER/$FOLDER/02a_mapped_modern
TASK=theta


#Run analysis
while IFS= read -r line; do $ATLAS task=theta bam=$line window=Regions.LR75_500kbwindow.bed > $line.theta.500kbwindow; done < $INPUT

```








## Read into R and plot


Download the resulting theta.txt.gz files to your home computer: 

```
/Users/alexjvr/2021postdoc/Velocity/E3_A.hyperantus/ATLAS

```

And plot. 

You'll need to know how your chromosomes and scaffolds are named. Have a look (using cat or less) at the indexed reference genome (xx.fna.fai) in the RefGenome folder on the server. You'll see several longer chromosomes (typically ~28), and some scaffolds that are shorter and start with a different name. 

e.g. for Ringlet (Aphantopus hyperantus). 

Chromosomes = LRxxx

Scaffolds = CADCXxxx


Brown Argus (Aricia agestis)

Chromosomes = LR

Scaffolds = CAJH


Speckled wood (Pararge aegeria)

Chromosomes = NC

Scaffolds = NW

###MUS
```
#Read population into a list
MUS.files <- list.files(pattern="AH-01-1900") 
MUS.myfiles <- lapply(MUS.files, read.table, header=T)  ##read all files into a list
colnames.new <- c("Chr", "start", "end", "depth", "fracMissing", "fracTwoOrMore", "pi.A", "pi.C", "pi.G", "pi.T", "theta_MLE", "theta_C95_l", "theta_C95_u", "LL")
MUS.myfiles2 <- lapply(MUS.myfiles, setNames, nm=colnames.new)   #rename columns

#Add pop and sample columns 
MUS.files <- gsub("_theta_estimates.txt.gz", "", MUS.files)
library(purrr)
MUS.myfiles <- Map(cbind, MUS.myfiles, Sample=MUS.files)
MUS.myfiles <- Map(cbind, MUS.myfiles, Pop="MUS")


#Convert list to dataframe
library(reshape2)
allvars <- colnames(MUS.myfiles[[1]])
MUS.ll <- melt(MUS.myfiles, id.vars=allvars)

#Add a column with the midposition of the estimate (calculated from the start and end of each window)
MUS.ll$midpos=((MUS.ll$end-MUS.ll$start)/2)+MUS.ll$start

#Keep only the chromosomes. Replace "LR" with your chromosome name. See above
MUS.ll.Chrsonly <- (MUS.ll %>% filter(grepl("LR", Chr)))


##Check the proportion of missing data: 
pdf("MUS.missingdata.pdf")
ggplot(MUS.ll.Chrsonly, aes(x=fracMissing, y=theta_MLE, colour=Chr))+geom_point()     #Plot of theta vs fraction missing data in each window coloured by chromosome
ggplot(MUS.ll.Chrsonly, aes(x=fracTwoOrMore, y=theta_MLE, colour=Sample))+geom_point()  #Plot of theta vs fraction of loci with 2x+ coverage
ggplot(MUS.ll.Chrsonly, aes(x=fracMissing, y=theta_MLE, colour=Sample))+geom_point()  #Plot of theta vs fraction missing data in each window, coloured by sample. 
dev.off()

##Select a threshold based on the previous plots that maximise data, but remove any poor quality individuals: 
##
Thresh <- 0.6  #replace this threshold based on the above figures. 
summary(MUS.ll.Chrsonly[which(ll.Chrsonly$fracMissing>Thresh),])

#List of samples to keep based on the above threshold
MUS.tokeep <- c("AH-01-1900-04", "AH-01-1900-05", "AH-01-1900-06", "AH-01-1900-08", "AH-01-1900-09", "AH-01-1900-10", "AH-01-1900-11", "AH-01-1900-13", "AH-01-1900-14", "AH-01-1900-15", "AH-01-1900-16", "AH-01-1900-20", "AH-01-1900-21", "AH-01-1900-22", "AH-01-1900-22", "AH-01-1900-23", "AH-01-1900-24", "AH-01-1900-25", "AH-01-1900-27", "AH-01-1900-28", "AH-01-1900-29", "AH-01-1900-32", "AH-01-1900-33", "AH-01-1900-34", "AH-01-1900-35", "AH-01-1900-37", "AH-01-1900-38", "AH-01-1900-39", "AH-01-1900-40", "AH-01-1900-41", "AH-01-1900-42", "AH-01-1900-43", "AH-01-1900-45", "AH-01-1900-46", "AH-01-1900-47")

MUS.ll.Chronly.0.6miss <- filter(MUS.ll.Chrsonly, Sample %in% MUS.tokeep)
summary(MUS.ll.Chronly.0.6miss$Sample)

```


###MODC
```
#Read population into a list
MODC.files <- list.files(pattern="AH-01-20") 
MODC.myfiles <- lapply(MODC.files, read.table, header=T)  ##read all files into a list
colnames.new <- c("Chr", "start", "end", "depth", "fracMissing", "fracTwoOrMore", "pi.A", "pi.C", "pi.G", "pi.T", "theta_MLE", "theta_C95_l", "theta_C95_u", "LL")
MODC.myfiles2 <- lapply(MODC.myfiles, setNames, nm=colnames.new)   #rename columns

#Add pop and sample columns 
MODC.files <- gsub("_theta_estimates.txt.gz", "", MODC.files)
library(purrr)
MODC.myfiles <- Map(cbind, MODC.myfiles, Sample=MODC.files)
MODC.myfiles <- Map(cbind, MODC.myfiles, Pop="MODC")


#Convert list to dataframe
library(reshape2)
allvars <- colnames(MODC.myfiles[[1]])
MODC.ll <- melt(MODC.myfiles, id.vars=allvars)
MODC.ll$midpos=((MODC.ll$end-MODC.ll$start)/2)+MODC.ll$start

#Keep only the chromosomes
MODC.ll.Chrsonly <- (MODC.ll %>% filter(grepl("LR", Chr)))

##Check the proportion of missing data: 
pdf("MODC.missingdata.pdf")
ggplot(MODC.ll.Chrsonly, aes(x=fracMissing, y=theta_MLE, colour=Chr))+geom_point()   
ggplot(MODC.ll.Chrsonly, aes(x=fracTwoOrMore, y=theta_MLE, colour=Sample))+geom_point()
ggplot(MODC.ll.Chrsonly, aes(x=fracMissing, y=theta_MLE, colour=Sample))+geom_point()
dev.off()

##Select a threshold and identify samples to remove (if any) 
summary(MODC.ll.Chrsonly[which(MODC.ll.Chrsonly$fracMissing>0.6),])

#MODC.tokeep <- c("")

#MODC.ll.Chronly.0.6miss <- filter(MODC.ll.Chrsonly, Sample %in% MODC.tokeep)
#summary(MODC.ll.Chronly.0.6miss$Sample)

```




###MODE
```
#Read population into a list
MODE.files <- list.files(pattern="AH-02") 
MODE.myfiles <- lapply(MODE.files, read.table, header=T)  ##read all files into a list
colnames.new <- c("Chr", "start", "end", "depth", "fracMissing", "fracTwoOrMore", "pi.A", "pi.C", "pi.G", "pi.T", "theta_MLE", "theta_C95_l", "theta_C95_u", "LL")
MODE.myfiles2 <- lapply(MODE.myfiles, setNames, nm=colnames.new)   #rename columns

#Add pop and sample columns 
MODE.files <- gsub("_theta_estimates.txt.gz", "", MODE.files)
library(purrr)
MODE.myfiles <- Map(cbind, MODE.myfiles, Sample=MODE.files)
MODE.myfiles <- Map(cbind, MODE.myfiles, Pop="MODE")


#Convert list to dataframe
library(reshape2)
allvars <- colnames(MODE.myfiles[[1]])
MODE.ll <- melt(MODE.myfiles, id.vars=allvars)
MODE.ll$midpos=((MODE.ll$end-MODE.ll$start)/2)+MODE.ll$start

#Keep only chromosomes
MODE.ll.Chrsonly <- (MODE.ll %>% filter(grepl("LR", Chr)))

##Check the proportion of missing data: 

pdf("MODE.missingdata.pdf")
ggplot(MODE.ll.Chrsonly, aes(x=fracMissing, y=theta_MLE, colour=Chr))+geom_point()   
ggplot(MODE.ll.Chrsonly, aes(x=fracTwoOrMore, y=theta_MLE, colour=Sample))+geom_point()
ggplot(MODE.ll.Chrsonly, aes(x=fracMissing, y=theta_MLE, colour=Sample))+geom_point()
dev.off()

##Select a threshold and identify samples to remove (if any) 
#summary(MODE.ll.Chrsonly[which(MODE.ll.Chrsonly$fracMissing>0.6),])

#MODE.tokeep <- c("")

#MODE.ll.Chronly.0.6miss <- filter(MODE.ll.Chrsonly, Sample %in% MODE.tokeep)
#summary(MODE.ll.Chronly.0.6miss$Sample)

```


Concat and Plot
```
# Concatenate the three data.frames
library(dplyr)

E3.data <- bind_rows(MUS.ll.Chrsonly, MODC.ll.Chrsonly, MODE.ll.Chrsonly)
dim(data)
[1] 47347    18



# Make chromosome a factor 
E3.data$Chr <- as.factor(E3.data$Chr)


# Plot
pdf("E3.theta.pdf")
ggplot(E3.data, aes(x=midpos, y=theta_MLE, colour=Pop))+geom_point(size=0.2)+facet_wrap(~Chr)
dev.off()



```





