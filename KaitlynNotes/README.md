# e-labbook 🦋

## Tuesday 06/07/2021 🦋
- Editing scripts to suit UCL server
     - Things to change in the current scripts:
          - path to wrapper or tools script
          - change PBS to #$
          - change array jobs to SGE_TASK_ID
          - remove 'module load ..'
          - VMEM added
- Edits made to **01a_modern_cutadapt_filtering_trimming.sh** and **01a_parallel_cutadapt_UCL.sh**:
     - made relevant changes with no syntax errors (using./) 
     - unsure how to submit script using qsub to check it's correct
- Edits made to **01b_concat_fastq_R1.sh** and **01b_concat_fastq_R2.sh**:
     - probably need to change *cd $PBS_O_WORKDIR* but not sure what to
     - is /newhome etc correct?
- Edits made to **01b_modern_trimmomatic_filtering_trimming.sh**: 
     - cannot find file **parallel_trimmomatic_bluecp3.sh** which this file should run to
- Edits made to **01c_bbtools_repair_museum_ARRAY.sh** and **01d_bbtools_merge_museum_ARRAY.sh**:
     - same run job in work directory code as before (unsure what to change this to)
     - commented out running java for now *'#/share/apps/java/bin/java'*, not sure if it's the right thing to replace *'module load...'* with
- Edits made to **02a_MapwithBWAmem.ARRAY.sh**:
     - commented out running bwa for now *'#/share/apps/genomics/bwa-0.7.17/bwa'*, not sure if it's the right thing to replace *'module load...'* with
- Edits made to **03a_variant_calling_bluecp.sh and 03a_call_SNVs_bluecp3.sh** in wrapper:
     - edited paths of module versions to be loaded
     - changed the path of *'CALLER'* to UCL server tools
     - commented out path to *samtools* as not sure if it replaces *SAMTOOLS='samtools'*
     - edited UCL server options 
     - trying to run in home directory:
          - ERROR: You must specify an input directory
          - ERROR: You must specify an output directory
          - ERROR: You must specify a reference file
- Edits made to **03b_summary_variant_calling.sh** and **03b_call_SNVs_summarize_bluecp3.sh**:
     - edited definitions of modules to be loaded
     - did not add VMEM as no other memory variable was present (maybe it's not needed)
- Have now edited all pathways from pipeline that directly lead to tools & wrapper but there are still files in both tools & wrapper because wrapper can lead to tools:
- Edits made to **popgenstats.pl**:
     - edited paths to external programs (don't know what estpEM software is or whether its current path is correct)
     - permission denied to run

## Wednesday 07/07/2021 🦋
- Edited **01a_modern_cutadapt_filtering_trimming.sh** by comparing to museum script
     - no errors using ./ in the home directory but cannot submit yet as we need to upload data 
- Edits made to **01b_concat_fastq_R1.sh** include:
     - changed from running to the work directory of UoB to running to the home directory of UCL server pchuckle
     - edited the scheduler directives to *$* instead of *PBS*
     - removed nodes scheduler directive and added *-l vmem* and *-S /bin/bash* 
     - changed from *-j oe* to *-j y* in scheduler directives
     - changed the definition of *SPECIESDIR* to the current pathway to the E3_Aphantopus_hyperantus directory
     - changed the array job submissions to *{SGE_TASK_ID}*, current form is *sed "${SGE_TASK_ID}q;d"* but may need to change to sed *-n ${SGE_TASK_ID}'{q:d}' input.data* form
     - **NOTE:** when running using ./ I get the error *cat: //SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/01a_museum_cutadapt_reads/: No such file or directory*, so maybe it's because the data hasn't been uploaded yet or is it because 01a_museum_cutadapt_reads is in a different folder?
- Edits made to **01b_concat_fastq_R2.sh** include:
     - same as above for **01b_concat_fastq_R1.sh**
     - **NOTE:** syntax error looking for missing " , but it is written in the same way as R1 which doesn't have this syntax error ??

## Thursday 08/07/2021 🦋
- Edits made to **01b_concat_fastq_R1.sh** include:
     - change *vmem* to *h_vmem* for compatibility to UCL server
     - changed the pathway for 01a_museum_catadapt_reads to *local/...* but error of *No such file or directory* still for 01a_museum2_catadapt_reads
     - changed walltime to *h_rt* in the scheduler descriptive
     - qsub the job after deleting path for 01a_museum2_catadapt_reads which does create 33 sample files but no log files 
     - changed work directory to *SGE_O_WORKDIR*
- Edits made to **01b_concat_fastq_R2.sh** include:
     - changed walltime to *h_rt* in scheduler descriptive
     - changed pathway for 01a_museum_catadapt_reads and 01a_museum2_catadapt_read to *local/...*
     - **NOTE:** still the same syntax error but cannot identify where as it is the exact same as 01b_concat_fastq_R1.sh but with R1 changed to R2
     - changed work directory to *SGE_O_WORKDIR*
- Edits made to **01c_bbtools_repair_museum_ARRAY.sh** include:
     - changed scheduler descriptives from *#PBS* to *#$*
     - deleting nodes descriptive and added *h_vmem* descriptive
     - change *mem* to *tmem*, *oe* to *y* and *walltime* to *h_rt*
     - add */bin/bash* descriptive
     - changed work directory to *SGE_O_WORKDIR*
     - define path to bbrepair through bbmap
     - changed path to load java
     - **NOTE:** errors of invalid maximum heap size and could not create the Java Virtual Machine
     - changed maximum heap size using *export _JAVA_OPTIONS="-Xmx2g* but still the same error occurs of invalid max heap size
- Edits made to **01d_bbtools_merge_museum_ARRAY.sh** include:
     - same as above for **01c_bbtools_repair_museum_ARRAY.sh** for all the scheduler descriptives
     - edited array jobs to *SGE_TASK_ID*
     - edited path to bbmerge
     - **NOTE:** can't read *R1.museum.names.repaired: No such file or directory*
- Edits made to **02a_MapwithBWAmem.ARRAY.sh** include:
     - change from *#PBS* to *#$*
     - delete nodes descriptive and add *h_vmem* descriptive
     - change *mem, walltime* and *oe* to *tmem, h_rt* and *y*
     - add *-S /bin/bash* descriptive
     - change working directory
     - load bwa module 
     - change array jobs to *SGE_TASK_ID*
     - **NOTE:** error *'demultiplexed/': No such file or directory*
- [Created a pipeline script for museum2 samples to use cutadapt](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/01a_museum2_cutadapt_filtering_trimming.sh) 
- [Created a pipeline script for modern_exp to use cutadapt](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/01a_modern_exp_cutadapt_filtering_trimming.sh)
- [Created a pipeline script for museum 2 sample to use fastqc](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/00_fastqc_raw_museum2.sh)
- Successfully ran concat script (I think....!)

## Friday 09/07/2021 🦋
- Edits made to **03a_variant_calling_bluecp.sh** include:
     - renamed file to **03a_variant_calling_UCL.sh**
     - edited path to samtools
     - changed definition of working directory
     - **NOTE:** will need to change input and output files 
- Changed all the permissions in VelocityPipeline folders so that everyone can execute all files
- Edits made to **03a_call_SNVs_bluecp.sh** include:
     - renamed file to **03a_call_SNVs_UCL.sh**
     - edited UCL server options such that
          - echo '#$ -S /bin/bash' >> $SMSJOB\
          echo '#$ -N '$JOBNAME'' >> $SMSJOB\
echo '#$ -l tmem='$MEM'G' >> $SMSJOB\
echo '#$ -l h_vmem='$VMEM'G' >> $SMSJOB\
echo '#$ -l h_rt='$HRS'1:00:00' >> $SMSJOB\
echo '#$ -j y' >> $SMSJOB #concatenates error and output files (with prefix job1)\
echo '#$ -t 1-'$N >> $SMSJOB\
echo '#$ -o '$LOG >> $SMSJOB
- Edits made to **01c_bbtools_repair_museum_ARRAY.sh** include:
     - added a *PATH* variable for a newly created folder containing the concatenated files 
     - changed the infiles to include *$PATH*
     - created a new folder in my local directory called 01c_musPERepaired which will be where the outfiles of **01c_bbtools_repair_museum_ARRAY.sh** will go
     - **NOTE:** outfile appears empty after 'qsub'bing the script 

## Monday 12/07/2021 🦋
- Edits made to **01c_bbtools_repair_museum_ARRAY.sh** include:
     - moved R2 concat files to 01a.mus.concat.cutadapt folder
     - coded for the creation of sample files
     - **NOTE:** now successfully running
- Edits made to **01d_bbtools_merge_museum_ARRAY.sh** include:
     - added variables for the infile and outfile paths
     - created a new folder in local for the outfile
     - coded for the creation of sample names files
     - edited the *PREFIX* for the outfiles
- Edits made to **02a_MapwithBWAmem.ARRAY_museum.sh**:
     - changed job name and walltime in scheduler descriptive
     - edited paths to the reference, and input and output directories
     - coded for the creation of sample file names to map
     - hashed out NAME2 and deleted R2 paths 
     - **NOTE:** takes a long time to run. 
          - UPDATE: finished running with a total size of 7.3G
- Running scripts through from the start and transferring files across to shared folder:
     - Successfully ran **00_fastqc_raw_museum.sh** and **00_fastqc_raw_museum2.sh**
     - successfully ran **01a_museum_cutadapt_filtering_trimming.sh** and **01a_museum2_cutadapt_filtering_trimming.sh**
     - successfully ran **01b_concat_fastq_R1.sh and 01b_concat_fastq_R2.sh** but manually move both sets of outfiles into a folder called 01a.mus.concat.cutadapt
     - **01c_bbtools_repair_museum_ARRAY.sh** does not seem to be running after submitting the job
          - UPDATE: it simply takes quite a while to run
     - succesfully ran **01c_bbtools_repair_museum_ARRAY.sh**
     - succesfully ran **01d_bbtools_merge_museum_ARRAY.sh**
     - succesfully ran **02a_MapwithBWAmem.ARRAY_museum.sh** 

## Tuesday 13/07/2021 🦋
- Running all scripts on the full museum dataset:
     - Successfully ran **00_fastqc_raw_museum.sh** and **00_fastqc_raw_museum2.sh**
     - created a folder in the shared directory called **00_raw_reads_museum_FINAL**
     - successfully ran **01b_concat_fast1_R1.sh** and **01b_concat_fastq.R2.sh** on the files which were sequences twice
     - moved all the files that were sequenced once to 00_raw_reads_museum_FINAL
     - successfully ran **01a_museum_cutadapt_filtering_trimming.sh**
     - ran **01c_bbtools_repair_museum_ARRAY.sh** but a few files are missing
     - successfully ran **01d_bbtools_merge_museum_ARRAY.sh**
     - successfully ran **02a_MapwithBWAmem.ARRAY_museum.sh**
- Reading on selective sweeps and made handwritten notes on models for selective sweeps and tests for sweep detection

## Wednesday 14/07/2021 🦋:
- Re-ran **01c_bbtools_repair_museum_ARRAY.sh** but with the lines in the script to create the sample file names hashed out which now seems to run successfully with no files missing
- Successfully re-ran **01d_bbtools_merge_museum_ARRAY.sh** and **02a_MapwithBWAmem.ARRAY_museum.sh** without any missing files
- [Created a markdown file of useful bash commands](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/Bash.md) which I will continuously add to
- Edited README.md so that it is more readable
- Create an array to obtain mapping stats
     - *ls /SAN/ugi/LepGenomics/C3_Aricia_agestis/02a_mapped_museum/\*.gz\* | awk -F "/" '{print $NF}' > museum.toflagstat*
     - *NAME=$(sed "${SGE_TASK_ID}q;d" museum.toflagstat)*
     - */share/apps/genomics/samtools-1.9/bin/samtools flagstat ${NAME}.sorted.bam >> flagstat.log*
- Check size of sam file size in **02a_MapwithBWAmem.ARRAY_museum.sh**
     - export *PATH=/share/apps/genomics/bcftools-1.9/bin:$PATH*
     - export *LD_LIBRARY_PATH=/share/apps/genomics/bcftools-1.9/lib:$LD_LIBRARY_PATH*
     - **ERROR**: unknown file type
     - **NOTE**: in main README it checks the size of bam files but we have sam files outputted here
- Load samtools module to gather statistics from **02a_MapwithBWAmem.ARRAY_museum.sh** data:
     - *export PATH=/share/apps/genomics/samtools-1.9/bin:$PATH*
     - *export LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH*
     - *samtools flagstat file.sam*
          - **NOTE**: returning 0 + 0 for all variables
     - converting sam file to bam file and then running flagstat:
          - *samtools view -S -b file.sam > file.bam*
          - *samtools flagstat file.bam*
          - **NOTE**: again returning all variables as 0 + 0  
     - made a flagstat log file for all of the samples
          - for i in $(ls *bam); do ls $i >>flagstat.log && samtools flagstat $i >> flagstat.log; done
- Running modern data through the scripts:
     - successfully ran **01a_modern_cutadapt_filtering_trimming.sh**
- Running modern exp data through the scripts:
     -  successfully ran **01a_modern_exp_cutadapt_filtering_trimming.sh**

## Thursday 15/07/2021 🦋
- Edits made to **02a_MapwithBWAmem.ARRAY_modern.sh.save** include:
     - edited the input and output paths
     - defined array names
     - coded for the ls of the R1 and 42 cutadapt files 
     - added *${NAME2}* to line 52 starting *echo*
     - write output to bam instead of sam
          - bwa mem $ref input1 input2 | samtools sort -o output.bam
     - export samtools module into script
- Created a new directory in C3_Aricia_agestis called 02a_mapped_modern
- **NOTE**: flagstat showing a really low mapping rate of the museum data to the reference genome
- Successfully ran **02a_MapwithBWAmem.ARRAY_modern.sh.save** and flagstat showed that the mapping rates were very good >95%
- Edits made to **02a_MapwithBWAmem.ARRAY_museum.sh** include:
     - write output to bam instead of sam
     - export samtools module into script
     - new reference genome path
          - do indexing step separately /share/apps/genomics/bwa-0.7.17/bwa index /SAN/ugi/LepGenomics/C3_Aricia_agestis/RefGenome/GCA_905147365.1_ilAriAges1.1_genomics.fna 
- Mapping the test files straight from the raw data:
     - with the old reference genome mapping rates were ~75%
     - with the new reference genome mapping rates were also ~75%
- (ALWAYS USING NEW REF GENOME FROM NOW)
- Mapping the test files from raw data &#8594; cutadapt data:
     -  mapping rates were between 60-66%
- Mapping the test files from raw data &#8594; cutadapt data &#8594; concat data: (JOB ID 3564892)
     - mapping rates were ranging from 74-88%
- Mapping the test files from raw data &#8594; concat data: (JOB ID 3564895)
     - mapping rates were ~75%
- Mapping the test files from raw data &#8594; concat data &#8594; cutadapt data: (JOB ID 3564897)
     - mapping rates were ranging from 74-88% (exactly the same as above)
- Running overnight **02a_MapwithBWAmem.ARRAY_modern.sh.save** with the new reference genome
- Running overnight **02a_MapwithBWAmem.ARRAY_modern.sh.save** for the modern_exp data 
- Run sample files of E3 species to check for any errors along the pipeline in terms of mapping rates:
     - Added 3 samples that needed concatenating from mus and mus2 to separate directories
     - Successfully ran **01b_concat_fastq_R1.sh** and **01b_concat_fastq_R2.sh**
     - Successfully ran **01a_museum_cutadapt_filtering_trimming.sh**

## Friday 16/07/2021 🦋
- Continuing to run the sample files for the E3 species through the scripts:
     - Successfully ran **01c_bbtools_repair_museum_ARRAY.sh**
     - Successfully ran **01d_bbtools_merge_museum_ARRAY.sh**
     - Successfully ran **02a_MapwithBWAmem.ARRAY_museum.sh**
          - exported samtools 
          - used flagstat to check the mapping rates which were 90-95%, and so far better than those from the C3 museum data suggesting that the pipeline scripts are working as expected and so the quality of the C3 data is comparatively low
- Checked flagstat mapping rates for modern and modern_exp bam files obtained from the new ref genome
     - made a flagstat log file for each of modern and modern_exp
          -  for i in $(ls *bam); do ls $i >>flagstat.log && samtools flagstat $i >> flagstat.log; done
     - mapping rates were good ~99%
     - some files had 0 reads so re-running those 
     - some samples had multiple temp files meaning bwa was not completed so re-running those 
     - collected data in shared excel sheet of number of reads and read quality at this step
- Updated my bash commands list 

## Monday 19/07/2021 🦋
- Request installation of ANGSD 
- Missing 13, 15 modern mapped files, so re-running those
- Missing 33, 36, & 40 modern exp mapped files, so re-running those
- Filling in the missing statistics from the sample files that were showing 0 reads last time after having re-run them 
- Edits made to **03a_variant_calling_UCL.sh** include:
     - updating the input and output paths
     - updating the reference 
     - added -m and -l (MEM and VMEM) variable conditions 
     - specified all new variables
          - SHAREDFOLDER="/SAN/ugi/LepGenomics"\
SPECIES="C3_Aricia_agestis"\
INPUT="02a_mapped_modern"\
OUTPUT="03.2_Variant_calling_museum"\
REF="RefGenome/GCA_905147365.1_ilAriAges1.1_genomic.fna"
SAMTOOLS="/share/apps/genomics/samtools-1.9/bin/samtools"
JOBNAME="C3_mod_mpileup"
CALLER="/SAN/ugi/LepGenomics/VelocityPipeline/wrapper/03a_call_SNVs_UCL.sh"
- ran **03a_variant_calling_UCL.sh**:
     - **NOTE:** indexing step takes a while
     - **ERROR:** script was calling a UCL.pl file which we hadn't re-named yet in the UCL pipeline, so I renamed that .pl file
- added ANGSD path to main README markdown
- ran *samtools flagstat* on the modern exp files that needed re-running and updated excel sheet with statistics
- Set-up the ANGSD analysis:
     - Find all the regions (i.e. all chromosomes and contigs) from the reference index file (.fasta.fai):
          - awk '{print $1}' ../RefGenome/*.fna.fai >> regions
          - cat regions |wc -l 
          - there are 28 regions

## Tuesday 20/07/2021 🦋
- Running new version of **03a_variant_calling_UCL.sh**:
     - several syntax error to be corrected 
     - compare and contrast with old version of the script
     - replicated the new script under a new name which didn't result in any syntax errors even though it's exactly the same ??
- Edits made to **03a_variant_calling_UCL.sh**:
     - removed *module1 $SAMTOOLS* from the variables list
- **ERROR:** cannot locate Parallel/ForkManager or Term/ProgressBar (maybe need to be installed on the server)
- Added export paths to perl package and library to **03a_call_SNVs_UCL.sh** but still getting the same error of cannot locate the perl modules
- Added the last sample's (40 of modern exp) statistics to the excel sheet using *samtools flagstat*
- Added *gcc* path to main README markdown doc
- Requsted installation of perl modules *Parallel:ForkManager* and *Term:ProgressBar*
- Estimate SAF for each population (unfolded):
     - mod core population:
          - */share/apps/genomics/angsd-0.935/bin/angsd -b /SAN/ugi/LepGenomics/C3_Aricia_agestis/02a_mapped_modern -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r LR990278.1: -GL 1 -doSaf 1 -anc /SAN/ugi/LepGenomics/C3_Aricia_agestis/RefGenome/GCA_905147365.1_ilAriAges1.1_genomic.fna -ref /SAN/ugi/LepGenomics/C3_Aricia_agestis/RefGenome/GCA_905147365.1_ilAriAges1.1_genomic.fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth 144 -doMajorMinor 4 -out MODC -C 50 -baq 1*
     - mod exp population:
          - */share/apps/genomics/angsd-0.935/bin/angsd -b /SAN/ugi/LepGenomics/C3_Aricia_agestis/02a_mapped_modern -checkBamHeaders 1 -minQ 20 -minMapQ 20 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -r LR990278.1: -GL 1 -doSaf 1 -anc /SAN/ugi/LepGenomics/C3_Aricia_agestis/RefGenome/GCA_905147365.1_ilAriAges1.1_genomic.fna -ref /SAN/ugi/LepGenomics/C3_Aricia_agestis/RefGenome/GCA_905147365.1_ilAriAges1.1_genomic.fna -doCounts 1 -setMinDepthInd 2 -setMaxDepth 144 -doMajorMinor 4 -out MODE -C 50 -baq 1*
     - commands taking very long time to run and cat on the arg files are resulting as empty

## Wednesday 21/07/2021 🦋
- Created new directory in 02a_mapped_modern containing four test bam files to run ANSGD on
- Created a script called [**03a_estimateSAF**](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/03a_estimateSAF) to estimate SAF for each population
- Edits made to **03a_estimateSAF** include:
     - change *-only_proper_pairs* to 1 for modern populations
     - use the smallest chromosome
- Export ANSGD into the server:
     - *export PATH=/share/apps/genomics/angsd-0.935/bin:$PATH*
     - *export LD_LIBRARY_PATH=/share/apps/genomics/angsd-0.935/lib:$LD_LIBRARY_PATH*
- Created new directory in 02a_mapped_modern_exp containing four test bam files to run ANSGD on
- created files called MODCTEST.poplist and MODETEST.poplist containing a list of the test bam files to run through ANGSD
     - *ls /SAN/ugi/LepGenomics/C3_Aricia_agestis/02a_mapped_modern_exp/TEST/*bam* >> MODETEST.poplist*
- used *samtools index [filename]* on the samples I moved into new file
- successfully ran **03a_estimateSAF** for both MODC and MODE populations
- Next step is for unfolded SAF to be used to produce folded 2D SFS:
     - generates the folded SFS for each population pair per chromosome
     - [INSERT CODE HERE]
- Request installation of ANGSD program RealSFS

## Thursday 22/07/2021 🦋
- 13 individuals in the 01d_musALL_merged didn't run properly with file sizes of ~32K 
     - individuals to re-run are 07, 09, 11, 13, 14, 15, 19, 20, 22, 23, 25, 27, 31 
     - re-ran successfully with file sizes as expected
- Read articles on population genomic analyses for modern and museum samples
- Re-run **02a_MapwithBWAmem.ARRAY_museum.sh** to include the newly merged samples
- Example code to un-zip and view red .gz files created after running **03a_estimateSf**:
     - *gunzip MODCTEST.beagle.gz |head*
     - *head MOCTEST.beagle*
- Compute estimates of mean depths of the BAM files in each population:
     - export *samtools* first
     - example of one file:
          - samtools depth AAg-19-2016-01_L002_cutadapt_filtered_R1.fastq.gz.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
     - for loop:
          - for i in $(ls *bam); do ls $i >>depth.log && samtools depth $i |awk '{sum+=$3} END { print "Average = ",sum/NR}' >> depth.log; done
- Perform PCA using ANGSD on the BAM files for each population:
     - Export R 
          - export PATH=/share/apps/R-4.0.3/bin:$PATH
          - export LD_LIBRARY_PATH=/share/apps/R-4.0.3/lib64:$LD_LIBRARY_PATH
     - Find the number of loci and number of chromosomes in each population:
          ```
          MODCTEST <- read.table(gzfile("MODCTEST.beagle.gz"), header=T)
          dim(MODCTEST)
          [1] 7347910      15
          MODETEST <- read.table(gzfile("MODETEST.beagle.gz"), header=T)
          dim(MODETEST)
          [1] 7123055      15
          ```
     - Identify loci within MODE population that have <80% genotyping rate and remove them:
          ```
          library(reshape2)
          imiss_MODETEST <- MODETEST
          imiss_MODETEST[imiss_MODETEST=="0.333333"] <- NA
          mean(is.na(imiss_MODETEST)) ## calculates the overall proportion of missingness
          [1] 0.2663371
          
          MODE.loci.propNA <- rowMeans(is.na(imiss_MODETEST))*100 #calculates the proportion of NA in each row, i.e. for each locus
          melt.MODE.loci.propNA <- melt(MODE.loci.propNA)
          melt.MODE.loci.propNA$marker <- MODETEST$marker
          dim(melt.MODE.loci.propNA[which(melt.MODE.loci.propNA$value>20),]) #loci with >20% missingness
          [1] 2956687       2
          
          MODETEST.clean <- MODETEST[which(!MODETEST$marker %in% MODE.markerstoremove$marker),] #remove problematic loci
          dim(MODETEST.clean)
          [1] 4166368      15
          
          imiss_MODETEST <- MODETEST.clean
          imiss_MODETEST[imiss_MODETEST=="0.333333"] <- NA
          MODE.indiv.propNA <- colMeans(is.na(imiss_MODETEST))*100 #calculates the proportion of NA in each column, i.e. for each individual
          melt.MODE.indiv.propNA <- melt(MODE.indiv.propNA)  #plotting is easier with long data
          melt.MODE.indiv.propNA$Indiv <- rownames(melt.MODE.indiv.propNA)
          melt.MODE.indiv.propNA.new <- melt.MODE.indiv.propNA[seq(4,nrow(melt.MODE.indiv.propNA),3)] #get rid of the first 3 lines, and keep only one row per indiv. 
          melt.MODE.indiv.propNA.new[which(melt.MODE.indiv.propNA.new$value>20),]
          [1] value Indiv
          <0 rows> (or 0-length row.names)
          
          melt.MODE.indiv.propNA.new
          value Indiv
          Ind0 16.475645  Ind0
          Ind1  9.640435  Ind1
          Ind2 11.484079  Ind2
          Ind3 18.728974  Ind3
          ```
          - Only 4 individuals with more than 20% missing data so keep these in for now
     - Identify loci within MODC population that have <80% genotyping rate and remove them:
          - Same code as for the MODE population, but replace MODE with MODC 
          - loci with 20% missingness: 2685671
          - loci remaining about problematic loci removed: 4662239
          - four indivs with over 20% missing data 

## Friday 23/07/2021
- Updated the excel sheet:
     - with average depth statistics for all populations
     - with reads and mapping rates for museum BAM files
- Re-ran same 13 missing samples as before but outputted to 02a_mapped_museum.new
- Continue to perform PCA using ANGSD on the BAM files for modern populations:
     - **ERROR:** 'Cannot allocate memory'
     - Continue code on a development cluster *qrsh -l tmem=62G,h_vmem=62G*
     - Join the cleaned data by site and allele
          ```
          library(dplyr)
          MODE.marker <- as.data.frame(MODETEST.clean$marker)
          colnames(MODE.marker) <- "marker"
          MODC.marker <- as.data.frame(MODCTEST.clean$marker)
          colnames(MODC.marker) <- "marker"
          POP2clean.markers <- intersect(MODC.marker, MODE.marker)
          dim(POP2clean.markers)
          [1] 3241513       1
          
          ##Subset each dataset to keep only the overlapping markers
          MODC.clean.sub  <- MODCTEST.clean[which(MODCTEST.clean$marker %in% POP2clean.markers$marker),] 
          MODE.clean.sub  <- MODETEST.clean[which(MODETEST.clean$marker %in% POP2clean.markers$marker),]
           
          ##Join datasets together
          pops2.clean <- left_join(MODC.clean.sub, MODE.clean.sub, by="marker", suffix=c(".c", ".e"))
          cols.toremove <- as.data.frame(grep("allele", colnames(pops2.clean)))[-(1:2),]
          cols.toremove
          [1] 16 17
          pops2.clean2 <- pops2.clean[-cols.toremove]
          ncol(pops2.clean2)
          [1] 27
          ```
     - Need to join with the museum data too (skip this for now)
  
## Monday 02/08/2021 🦋
- clone discoal repository into shared folder 
    - *git clone https://github.com/kr-colab/discoal.git*
- install diploS/HIC 
    - *wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh*
    - *bash Anaconda3-5.0.1-Linux-x86_64.sh*
    - *export PATH=/share/apps/anaconda3-5/bin:$PATH*
    - *export LD_LIBRARY_PATH=/share/apps/anaconda3-5/lib:$LD_LIBRARY_PATH*
    - *pip install tensorflow*
         - **ERROR**: Command "/share/apps/anaconda3-5/bin/python -u -c "import setuptools, tokenize;__file__='/tmp/pip-build-ntth0we5/grpcio/setup.py';f=getattr(tokenize, 'open', open)(__file__);code=f.read().replace('\r\n', '\n');f.close();exec(compile(code, __file__, 'exec'))" install --record /tmp/pip-6jcu2ki2-record/install-record.txt --single-version-externally-managed --compile" failed with error code 1 in /tmp/pip-build-ntth0we5/grpcio/
    - *pip install keras*
         - **ERROR**: PermissionError: [Errno 13] Permission denied: '/share/apps/anaconda3-5/lib/python3.6/site-packages/Keras-2.4.3.dist-info'
         - Cache entry deserialization failed, entry ignored
    - *git clone https://github.com/kern-lab/diploSHIC.git*
    - *cd diploSHIC*
    - *python setup.py install*
         - **ERROR**: error: could not create '/share/apps/anaconda3-5/lib/python3.6/site-packages/shicstats.cpython-36m-x86_64-linux-gnu.so': Permission denied
- run code in schmep node:
     - *qrsh -l tmem=14G,h_vmem=14G*
     - re-try installing keras and tensorflow
          - keras installed successfully
          - tensorflow outputting some error messages still
     - re-try running *python setup.py install*
          - permission granted but disk quota exceeded  
 

              
