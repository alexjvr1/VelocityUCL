# ATLAS

[ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Home) is a pipeline developed in Daniel Wegmann's lab to process raw bam files from ancient DNA to obtain accurate variant calls and estimates of genetic diversity. 
The pipeline takes into account post-mortem damage (PMD) and low sequence coverage, and recovers variatns more accurately than the state-of-the-art method: MapDamage + GATK. 

ATLAS can be used to estimate genetic diversity. We can also generate input for ANGSD using ATLAS. 


## Test-species: 

E3 Aphantopus hyperantus

I'm testing the raw data pre-processing and ATLAS pipeline on the UCL server (CS) 




What is the state of the art with poolseq data? Should we use these approaches in windows for our data? 



## Running Pipeline: 
### 0. Concat raw museum samples 

33 individuals from each museum species has been sequenced twice to increase sequence depth. We're concatenating these raw data together, then moving all samples to a folder called 01a_raw_museum_FINAL

Since the concat step runs fairly quickly we'll submit these commands sequentially rather than as an array. (Too many small jobs can slow down or crash jobs on the server, so arrays should only be used for bigger jobs >1 hour each).


Use these three scripts in this order: 

[00a_CreateInputs.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/00a_CreateInputs.sh)

[00b_ConcatFastq.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/00b_ConcatFastq.sh)

[00c_CollateAllMusSamples.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/00c_CollateAllMusSamples.sh)

### 1. Remove adapter sequence using Cutadapt

We're using a script for each population. Run these in the working directory to create the submission script in each case: 

```
./script.sh
```

[01a_MUS_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/01a_MUS_cutadapt_filtering_trimming.sh)

[01a_MODC_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/01a_MODC_cutadapt_filtering_trimming.sh)

[01a_MODE_cutadapt_filtering_trimming.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/01a_MODE_cutadapt_filtering_trimming.sh)



### 2a. Map to Reference genome

Note that to use ATLAS correctly we will map unmerged reads.

ATLAS has an in-built function to annotate and merge these reads after mapping, and uses this sequence information for the PMD recalibration step. 


[02a_MapwithBWAmem.ARRAY_MUS_forATLAS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02a_MapwithBWAmem.ARRAY_MUS_forATLAS.sh)

[02a_MapwithBWAmem.ARRAY_MODC_forATLAS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02a_MapwithBWAmem.ARRAY_MODC_forATLAS.sh)

[02a_MapwithBWAmem.ARRAY_MODE_forATLAS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02a_MapwithBWAmem.ARRAY_MODE_forATLAS.sh)



### 2b.0 Add Read Groups & 02b.1 MarkDuplicates

We'll use PicardTools and GATK to add read group information to the bam files, and then to mark any PCR duplicates. 


A script for each of the populations. The input files and RG information is changed for each: 

[02b.0_AddRG_MUS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.0_AddRG_MUS.sh)

[02b.0_AddRG_MODC.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.0_AddRG_MODC.sh)

[02b.0_AddRG_MODE.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.0_AddRG_MODE.sh)


Then mark duplicates

[02b.1_MarkDup_MUS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.1_MarkDup_MUS.sh)

[02b.1_MarkDup_MODC.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.1_MarkDup_MODC.sh)

[02b.1_MarkDup_MODE.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.1_MarkDup_MODE.sh)




### 02b.2 Local realignment using GATK3.8

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


Use GATK3.8 for local realignment using:

[02b.2_LocalRealignment_MUS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.2_LocalRealignment_MUS.sh)

[02b.2_LocalRealignment_MODC.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.2_LocalRealignment_MODC.sh)

[02b.2_LocalRealignment_MODE.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.2_LocalRealignment_MODE.sh)





### 02b.3 Validate with PicardTools ValidateSamFile

This can be run interactively if there are only a few samples. Or use these scripts: 

[02b.3_ValidateSamFile_MUS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.3_ValidateSamFile_MUS.sh)

[02b.3_ValidateSamFile_MODC.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.3_ValidateSamFile_MODC.sh)

[02b.3_ValidateSamFile_MODE.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b.3_ValidateSamFile_MODE.sh)


Outputs are written to the 02a_mapped_* folders, and are names ${NAME}.validatesam. If there are no errors they'll simply say: "No errors found"


```
cd 02a_mapped_museum
cat *validatesam

No errors found
No errors found
No errors found

```



### 03a.1 ATLAS: SplitMerge





### 03a.2 ATLAS: PMD



### 03a.3 ATLAS: recal





### 10. ATLAS: global diversity



### 11. ATLAS: Individual heterozygosity




#### 12. ATLAS: Output ANGSD input







