# PE read merging

After plotting MapDamage profiles, I've found a strange distribution of G>A substitutions in the Ringlet data after merging PE reads. This pattern holds in the Brown Argus and Pararge aegeria museum data as well (run by Frances and Zhirou). 

![alt_txt][Fig1]

[Fig1]:https://user-images.githubusercontent.com/12142475/141117504-0e913f6f-fe53-4ccb-960f-d9c2220edf5b.png


Possible issues: 

1) Adapters haven't been filtered out properly

2) As are kept preferentially over Gs based on base quality during merging. Is there a bias in base quality? Are we more certain of As in general? 

3) Reads are merged before any error correction (e.g. filtering for low quality reads, or removing the first 10bp of the reads. 


## IGV

We can compare the sequences that worked (unmerged PE reads) and the elevated G>A using IGV to view the sequences. 

This can't be launched on the UCL servers, so I'll download the relevant AH-01-1900-02 bam files to my computer:


#Copy data
```
#Open a tunnel to copy data


#Create a folder to work in: 
/Users/alexjvr/2021postdoc/Velocity/E3_A.hyperantus/IGV

## Sample that worked. Reads merged after mapping. 
rsync -auve "ssh -p 3000" $i ajansen@localhost:/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/02a_mapped_museum/AH-01-1900-02.realn_mergedReads.bam .
##See log file here: less E3.MapDmg.mus.o3966910.2
##Here the reads were merged after mapping, but this isn't seen properly by MapDamage, so MapDamage assess only the inward facing PP reads: 
##WARNING Processed 5090528 paired reads, assumed to be non-overlapping, facing inwards and correctly paired; 4385412 of these were excluded as improperly paired.
#However, this results in the expected substitution frequencies. 
00:05:04 mapdamage.rescale INFO Expected substition frequencies before and after rescaling:
00:05:04 mapdamage.rescale INFO     C>T    0.0039    0.0027
00:05:04 mapdamage.rescale INFO     T>C    0.0015    0.0015
00:05:04 mapdamage.rescale INFO     G>A    0.0031    0.0031
00:05:04 mapdamage.rescale INFO     A>G    0.0014    0.0014

## Sample with elevated G>A
rsync -auve "ssh -p 3000" $i ajansen@localhost:/SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/TrimmomaticTest/AH-01-1900-02.realn.bam* .
#This was generated by merging data before mapping using AdapterRemoval. Using default --collapse function
14:31:58 mapdamage.rescale INFO Expected substition frequencies before and after rescaling:
14:31:58 mapdamage.rescale INFO     C>T    0.0033    0.0028
14:31:58 mapdamage.rescale INFO     T>C    0.0021    0.0021
14:31:58 mapdamage.rescale INFO     G>A    0.0071    0.0071
14:31:58 mapdamage.rescale INFO     A>G    0.0020    0.0020

##Find the regions to look at by focusing on PE unmerged regions processed for the working file first. Find read names 
```


So far the only sites that generate believable G>A are the non-overlapping PE reads filtered from the full dataset. 

Test G>A for the following: 

Prepare data: Trimmomatic followed by AdapterRemoval
a) Remove adapters
b) 

1) R1 only. And distribution of nucleotide substitutions

2) R2 only. And distribution of nucleotide substitutions

3) Should we remove ends?

4) 


#### Adapters

```
#Default (Nextera)
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
#Nextera library prep
AGATCGGAAGAGCACACGTCTGAACTCCAGTC	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
GTGTAGATCT	GTGTAGATCT	
#TruSeq (Used in the Velocity project)
AATGATACGGCGACCACCGAGATCTACAC	CAAGCAGAAGACGGCATACGAGAT
TACACTCTTTCCCTACACGACGCTCTTCCGATCT	GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
#Nextera and Illumina adapters
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
AATGATACGGCGACCACCGAGATCTACAC	CAAGCAGAAGACGGCATACGAGAT
```

### Variables and wd

```
qrsh -l tmem=5G h_vmem=5G
cd /SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/TrimmomaticTest

AdapterRemoval=/SAN/ugi/LepGenomics/Software/adapterremoval-2.3.2/build/AdapterRemoval
AdapterList=/SAN/ugi/LepGenomics/Software/adapterremoval-2.3.2/Barcodes/VelocityAllAdapters
```



## Test1

AH02 PE adapter trimming with list of default, Nextera, Illumina, and TruSeq adapters.

```
 $AdapterRemoval --file1 AH-01-1900-02.trimtest_1P.fastq.gz --file2 AH-01-1900-02.trimtest_2P.fastq.gz --adapter-list $AdapterList
```


OUT
```
[Adapter trimming]
RNG seed: 2646895129
Alignment shift value: 2
Global mismatch threshold: 0.333333
Quality format (input): Phred+33
Quality score max (input): 41
Quality format (output): Phred+33
Quality score max (output): 41
Mate-number separator (input): '/'
Trimming 5p: 0
Trimming 3p: 0
Trimming Ns: No
Trimming Phred scores <= 2: No
Trimming using sliding windows: No
Minimum genomic length: 15
Maximum genomic length: 4294967295
Collapse overlapping reads: No
Deterministic collapse: No
Conservative collapse: No
Minimum overlap (in case of collapse): 11


[Trimming statistics]
Total number of read pairs: 3163204
Number of unaligned read pairs: 22895
Number of well aligned read pairs: 3140309
Number of discarded mate 1 reads: 142
Number of singleton mate 1 reads: 0
Number of discarded mate 2 reads: 142
Number of singleton mate 2 reads: 0
Number of reads with adapters[1]: 18840
Number of reads with adapters[2]: 12
Number of reads with adapters[3]: 9257
Number of reads with adapters[4]: 5203
Number of reads with adapters[5]: 7809
Number of reads with adapters[6]: 205
Number of reads with adapters[7]: 0
Number of retained reads: 6326124
Number of retained nucleotides: 315878771
Average length of retained reads: 49.9324


```



## Test2 

AH02 PE adapter trimming with list of default, Nextera, Illumina, and TruSeq adapters, PLUS collapse-conservatively

```
$AdapterRemoval --file1 AH-01-1900-02.trimtest_1P.fastq.gz --file2 AH-01-1900-02.trimtest_2P.fastq.gz --adapter-list $AdapterList --trimqualities --collapse-conservatively
```


OUTPUT
```
[Adapter trimming]
RNG seed: 3011570292
Alignment shift value: 2
Global mismatch threshold: 0.333333
Quality format (input): Phred+33
Quality score max (input): 41
Quality format (output): Phred+33
Quality score max (output): 41
Mate-number separator (input): '/'
Trimming 5p: 0
Trimming 3p: 0
Trimming Ns: No
Trimming Phred scores <= 2: Yes
Trimming using sliding windows: No
Minimum genomic length: 15
Maximum genomic length: 4294967295
Collapse overlapping reads: Yes
Deterministic collapse: No
Conservative collapse: Yes
Minimum overlap (in case of collapse): 11


[Trimming statistics]
Total number of read pairs: 3163204
Number of unaligned read pairs: 22895
Number of well aligned read pairs: 3140309
Number of discarded mate 1 reads: 142
Number of singleton mate 1 reads: 0
Number of discarded mate 2 reads: 142
Number of singleton mate 2 reads: 0
Number of reads with adapters[1]: 18840
Number of reads with adapters[2]: 12
Number of reads with adapters[3]: 9257
Number of reads with adapters[4]: 5203
Number of reads with adapters[5]: 7809
Number of reads with adapters[6]: 205
Number of reads with adapters[7]: 0
Number of full-length collapsed pairs: 3128610
Number of truncated collapsed pairs: 86
Number of retained reads: 3197428
Number of retained nucleotides: 167111238
Average length of retained reads: 52.2643

```






