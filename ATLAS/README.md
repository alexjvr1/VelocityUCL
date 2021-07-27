# ATLAS

[ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Home) is a pipeline developed in Daniel Wegmann's lab to process raw bam files from ancient DNA to obtain accurate variant calls and estimates of genetic diversity. 
The pipeline takes into account post-mortem damage (PMD) and low sequence coverage, and recovers variatns more accurately than the state-of-the-art method: MapDamage + GATK. 

ATLAS can be used to estimate genetic diversity. We can also generate input for ANGSD using ATLAS. 


## Test-species: 

E3 Aphantopus hyperantus

I'm testing the raw data pre-processing and ATLAS pipeline on the UCL server (CS) 


## Pre-processing: 



2. Remove adapter sequence using Cutadapt

3. Map to genome using BWA mem   #Note that unmerged reads are mapped as ATLAS has an in-built function to annotate and merge these reads. We need to keep this information because ATLAS also uses these data for the PMD recalibration

4. Add RG & Remove duplicates  

5. Local realignment using GATK3.8

6. Validate bam file using PicardTools ValidateSamFile

7. ATLAS: splitMerge

8. ATLAS: PMD

9. ATLAS: recalibrate

10. ATLAS: output ANGSD inputs

11. ATLAS: estimate genetic diversity in windows


What is the state of the art with poolseq data? Should we use these approaches in windows for our data? 



## Running Pipeline: 
### 1. Concat raw museum samples 

33 individuals from each museum species has been sequenced twice to increase sequence depth. We're concatenating these raw data together, then moving all samples to a folder called 01a_raw_museum_FINAL

Use the [01a_ConcatMusRpts.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/01a_ConcatMusRpts.sh)



### 4. Add RG, remove duplicates

### 5. Local realignment using GATK3.8

First create a dictionary for the reference genome if this is not available yet.

I'm running this in the interactive screen: 
```
qrsh -l tmem=34G,h_vmem=34G

cd /SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/RefGenome

export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar

java -jar $PICARD CreateSequenceDictionary R=GCA_902806685.1_iAphHyp1.1_genomic.fna O=GCA_902806685.1_iAphHyp1.1_genomic.fna.dict

```


Use GATK3.8 for local realignment using [02b.2_LocalRealignment.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.2_LocalRealignment.sh)



### 6. Validate with PicardTools ValidateSamFile






