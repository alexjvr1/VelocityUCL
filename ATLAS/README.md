# ATLAS

[ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Home) is a pipeline developed in Daniel Wegmann's lab to process raw bam files from ancient DNA to obtain accurate variant calls and estimates of genetic diversity. 
The pipeline takes into account post-mortem damage (PMD) and low sequence coverage, and recovers variatns more accurately than the state-of-the-art method: MapDamage + GATK. 

ATLAS can be used to estimate genetic diversity. We can also generate input for ANGSD using ATLAS. 


## Test-species: 

E3 Aphantopus hyperantus

I'm testing the raw data pre-processing and ATLAS pipeline on the UCL server (CS) 





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

### 2. Remove adapter sequence using Cutadapt

We're using a script for each population. Run these in the working directory to create the submission script in each case: 

```
./script.sh
```

[01a_MUS_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/01a_MUS_cutadapt_filtering_trimming.sh)

[01a_MODC_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/01a_MODC_cutadapt_filtering_trimming.sh)

[01a_MODE_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/01a_MODE_cutadapt_filtering_trimming.sh)



### 3. Map to Reference genome

Note that to use ATLAS correctly we will map unmerged reads.

ATLAS has an in-built function to annotate and merge these reads after mapping, and uses this sequence information for the PMD recalibration step. 


[02a_MapwithBWAmem.ARRAY_MUS_forATLAS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02a_MapwithBWAmem.ARRAY_MUS_forATLAS.sh)

[02a_MapwithBWAmem.ARRAY_MODC_forATLAS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02a_MapwithBWAmem.ARRAY_MODC_forATLAS.sh)

[02a_MapwithBWAmem.ARRAY_MODE_forATLAS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02a_MapwithBWAmem.ARRAY_MODE_forATLAS.sh)



### 4. Add Read Groups and MarkDuplicates

We'll use PicardTools and GATK to add read group information to the bam files, and then to mark any PCR duplicates. 


A script for each of the populations. The input files and RG information is changed for each: 

[02b.0_AddRG_MUS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.0_AddRG_MUS.sh)

[02b.0_AddRG_MODC.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.0_AddRG_MODC.sh)

[02b.0_AddRG_MODE.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.0_AddRG_MODE.sh)



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






