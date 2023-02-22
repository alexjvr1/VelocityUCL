# ATLAS

[ATLAS](https://bitbucket.org/wegmannlab/atlas/wiki/Home) is a pipeline developed in Daniel Wegmann's lab to process raw bam files from ancient DNA to obtain accurate variant calls and estimates of genetic diversity. 
The pipeline takes into account post-mortem damage (PMD) and low sequence coverage, and recovers variatns more accurately than the state-of-the-art method: MapDamage + GATK. 

ATLAS can be used to estimate genetic diversity. We can also generate input for ANGSD using ATLAS. 


## Test-species: 

E3 Aphantopus hyperantus

I'm testing the raw data pre-processing and ATLAS pipeline on the UCL server (CS) 

Since many of the steps in this pipeline run fairly quickly we'll submit these commands sequentially rather than as an array. Too many small jobs can slow down or crash jobs on the server, so arrays should only be used for bigger jobs >1 hour each. Previous scripts written as arrays are kept here: [smallARRAYscripts](https://github.com/alexjvr1/VelocityUCL/tree/main/ATLAS/Scripts/smallARRAYscripts)



## Running Pipeline: 
### 0. Concat raw museum samples 

33 individuals from each museum species has been sequenced twice to increase sequence depth. We're concatenating these raw data together, then moving all samples to a folder called 01a_raw_museum_FINAL



Use these three scripts in this order: 

[00a_CreateInputs.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/00a_CreateInputs.sh)

[00b_ConcatFastq.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/00b_ConcatFastq.sh)

[00c_CollateAllMusSamples.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/00c_CollateAllMusSamples.sh)


These take ~40minutes to run for E3



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




### 2b Process bam files before using ATLAS

Step1. Add Read groups

Step2. Mark Duplicates

Step3. Local Realignment

Step4. Validate sam file



We're using GATK3.8 for the local realignment. First we need to create a dictionary for the reference genome if this is not available yet.

I'm running this in the interactive screen: 
```
qrsh -l tmem=34G,h_vmem=34G

cd /SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/RefGenome

export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar

java -jar $PICARD CreateSequenceDictionary R=GCA_902806685.1_iAphHyp1.1_genomic.fna O=GCA_902806685.1_iAphHyp1.1_genomic.fna.dict

```

Then run the processing scripts. There is a script for each of the populations. The input files and RG information is changed for each: 

[02b_processbams_MUS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b_processbams_MUS.sh)

[02b_processbams_MODC.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b_processbams_MODC.sh)

[02b_processbams_MODE.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/02b_processbams_MODE.sh)



Outputs are written to the 02a_mapped_* folders. The final validation files are named ${NAME}.validatesam. If there are no errors they'll simply say: "No errors found"

e.g.
```
cd 02a_mapped_museum
cat *validatesam

No errors found
No errors found
No errors found

```

If any errors are found at this point they need to be corrected. Run ValidateSam on previous versions of the bam file to establish where the error started. 

eg. 

MATE_NOT_FOUND error means sequences were removed at one step in the processing. Sequences should be marked but not removed if they are problematic (e.g. mark duplicates rather than remove duplicates). 

Mate name or information errors can be corrected with [FixMateInformation](https://gatk.broadinstitute.org/hc/en-us/articles/360036713471-FixMateInformation-Picard-)

ERROR:MISMATCH_MATE_ALIGNMENT_START	1

ERROR:MISMATCH_MATE_CIGAR_STRING	1

e.g for E3 MODE
```
pwd
/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/02a_mapped_modern_exp

cat *validatesam

#Most of the files have "No errors found", but for four samples we have: 
## HISTOGRAM	java.lang.String
Error Type	Count
ERROR:MISMATCH_MATE_CIGAR_STRING	1

#To find the files: 
grep "MISMATCH_MATE_CIGAR_STRING" *validatesam

AH-02-2019-47.validatesam:ERROR:MISMATCH_MATE_CIGAR_STRING	1
AH-02-2019-48.validatesam:ERROR:MISMATCH_MATE_CIGAR_STRING	2
AH-02-2019-58.validatesam:ERROR:MISMATCH_MATE_CIGAR_STRING	1
AH-02-2019-69.validatesam:ERROR:MISMATCH_MATE_CIGAR_STRING	1

#Then we can run FixMateInformation for these four samples

#Set path
export PATH=/share/apps/java/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/java/lib:$LD_LIBRARY_PATH
PICARD=/share/apps/genomics/picard-2.20.3/bin/picard.jar

#Run this for each of the four samples. 
java -jar $PICARD FixMateInformation \
       I=AH-02-2019-47.realn.bam \
       O=AH-02-2019-47.realn.fixed_mate.bam \
       ADD_MATE_CIGAR=true

#Run ValidateSam on the fixed files to check that this has worked. 

```


Rerun ValidateSam on the final bams to make sure there are no errors. 





## ATLAS

The version numbers decrease, so v.0.9 is the latest version (not v.1.0). 

To use ATLAS: 
```
ATLAS=/share/apps/genomics/atlas-0.9/atlas
export LD_LIBRARY_PATH=/share/apps/openblas-0.3.6/lib:/share/apps/armadillo-9.100.5/lib64:$LD_LIBRARY_PATH
```

### 04a.0a ATLAS: Find all read groups

We need to split all the bam files by read group. 

First find all the read groups in each population. Use the [04.0a_FindRGs.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04.0_FindRG.sh)

```
for i in $(ls *mergedReads.bam); do samtools view $i | grep "NC_" | awk -F ":" '{print $1, $2, $3, $4}' |sort |uniq; done

K00124 207 HKWH3BBXX 3
```

Once we have this information, we can make two RG.txt files: 

1. RG.txt which contains the RG for the population
```
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH
export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH

for i in $(ls *realn.bam); do samtools view -H $i |grep "@RG"; done 

nano RG.txt
D3mus  paired
```

2. ALL.RG.txt which contains the RG information for each of the sequencing lanes
```
nano ALL.RG.txt
K00124:207:HKWH3BBXX:3      HKWH3BBXX_lane3      RG1
```


### 04a.0 ATLAS: SplitMerge

**The latest version of ATLAS is 0.9, not 1.0. Ed upgraded to v.0.9 on the UCL shared apps folder so that I can run splitMerge and PMD on the museum samples. 

Create a text file with the ReadGroup names, cycle length (if single end), and paired/single. 



```
ATLAS=/share/apps/genomics/atlas-0.9/atlas
for i in $(ls *realn.bam); do $ATLAS task=splitMerge bam=$i; done
```



### 04a.1b Split bams into RGs

Split all the individual bam files into different RGs

1. Write all of the reads associated with a particular RG to file 

Modify the [04a.1b_WriteRG1_list.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04a.0b_WriteRG1_list.sh) for each RG for each population. 


2. Split each individual bam file by RG

Modify the [04a.1c_SplitBAMS_RG1.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04a.0c_SplitBAMS_RG1.sh) script for each RG within each population


3. Calculate the total size of all the bam files for each RG

```
du -sch *RG1*merged*bam
```

If the total bams is <10Gb for a read group, merge the bam files. Use the [MergeBAMS.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/MergeBAMS.sh) script


### 03a.2 ATLAS: PMD - Museum individuals


Estimate the PMD per sample for each of the museum samples. This needs to be run on the full museum samples (with all RGs combined). Don't filter the reads. 

Modify the [04b_ATLAS_MUS.pmd.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04b_ATLAS_MUS.pmd.sh) script

This will output a pmd correction for each individual.


### 03a.3 ATLAS: recal - All RGs

This should be run on each of the read groups within each of the populations. 

Modify the [04_ATLAS_recal.sh](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/Scripts/04b_ATLAS_recal.sh)

Files needed to run this: 

pmdFile *Empiric* for the museum bam file

File containing a set of invaraint sites

mergeTheseRGs.txt file with all the RGs within the 

e.g. for E3 museum: 
```
cat mergeTheseRGs.txt 
AH-01-1900-01 AH-01-1900-04 AH-01-1900-05 AH-01-1900-06	AH-01-1900-08	AH-01-1900-09	AH-01-1900-10	AH-01-1900-11	AH-01-1900-13	AH-01-1900-14	AH-01-1900-15	AH-01-1900-16	AH-01-1900-20	AH-01-1900-21	AH-01-1900-22	AH-01-1900-23	AH-01-1900-24	AH-01-1900-25	AH-01-1900-27	AH-01-1900-28	AH-01-1900-29	AH-01-1900-32	AH-01-1900-33	AH-01-1900-34	AH-01-1900-35	AH-01-1900-37	AH-01-1900-38	AH-01-1900-39	AH-01-1900-40	AH-01-1900-41	AH-01-1900-45	AH-01-1900-46	AH-01-1900-47
```


Run recal independently for each sample. In some populations we find that not all the samples run to completion or the run fails to find an appropriate recalibration score. See [here](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/recal_tests.md) for my exploration of these errors. 

Based on the above results (and based on communication with the developers) we will use the medial recal score for the samples that did work as the recal score for the samples that failed, given that at least 5 samples produced valid recal scores. 

##### Calculate median recal score

First find all the samples that have run to completion and that do not contain the unconverged/un-initiated values ("1.00000, 0.00000..." see [here](https://github.com/alexjvr1/VelocityUCL/blob/main/ATLAS/recal_tests.md))
```
grep -L "1.00000" *EM.txt >> worked 
for i in $(cat worked); do cat $i >> worked.EM.txt; done

#Remove all header lines
sed -i '/^readGroup/d' worked.EM.txt

#replace all commas with white space
sed -i 's:,:\t:g' worked.EM.txt

#split the file into fwd and reverse reads and cut the text columns
grep first worked.EM.txt |cut -f4-27 > worked.EM_first.txt
grep second worked.EM.txt |cut -f4-27 > worked.EM_second.txt
```

Calculate the median for all columns and write a median file
```
#open python and import pandas
python
import pandas as pd

##read in files
df1 = pd.read_csv("worked.EM_first.txt", sep="\t")
df2 = pd.read_csv("worked.EM_second.txt", sep="\t")

##Find median (50% in the table)
df1_med = df1.describe()
df2_med = df2.describe()

##write tables
df1_med.to_csv("worked.first.describe", sep=",")
df2_med.to_csv("worked.second.describe", sep=",")
```

Extract the median values and create a median_recalibrationEM.txt file
```
#Extract the median for fwd and reverse
grep "50%" worked.*.describe |cut -d, -f2-27 > worked.median

#Extract text for first three rows
cut -f1-3 worked.EM.txt |head -n 2 > worked.rownames

#paste files together
paste worked.rownames worked.median > median_recalibrationEM_1.txt

#create a colnames file
nano colnames
readGroup	mate	model	quality	position	context

#Add headers to EM file
cat colnames median_recalibrationEM_1.txt > median_recalibrationEM.txt

#replace the comma field separators where expected
sed -i 's/,/ /2;s/,/ /3' median_recalibrationEM.txt

```

Find all the samples that this file will apply to
```
awk -F "_" '{print $1}' worked > bamlist.worked
diff bamlist bamlist.worked | grep '^<' | sed 's/^<\ //' > samples_with_med_recal
```



### 10. ATLAS: global diversity



### 11. ATLAS: Individual heterozygosity




#### 12. ATLAS: Output ANGSD input







