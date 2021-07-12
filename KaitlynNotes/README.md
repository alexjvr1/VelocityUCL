# e-labbook

## Tuesday 06/07/2021
### editing scripts to suit UCL server
Things to change:

- path to wrapper or tools script
- change PBS to #$
- change array jobs to SGE_TASK_ID
- remove 'module load ..'
- VMEM added

#### 01a_modern_cutadapt_filtering_trimming.sh and 01a_parallel_cutadapt_UCL.sh
- made relevant changes with no syntax errors (using./) 
- unsure how to submit script using qsub to check it's correct

#### 01b_concat_fastq_R1.sh and 01b_concat_fastq_R2.sh
- probably need to change cd $PBS_O_WORKDIR but not sure what to
- is /newhome etc correct?

#### 01b_modern_trimmomatic_filtering_trimming.sh 
- cannot find file parallel_trimmomatic_bluecp3.sh which this file should run to

#### 01c_bbtools_repair_museum_ARRAY.sh and 01d_bbtools_merge_museum_ARRAY.sh
- same run job in work directory code as before (unsure what to change this to)
- commented out running java for now '#/share/apps/java/bin/java', not sure if it's the right thing to replace 'module load...' with

#### 02a_MapwithBWAmem.ARRAY.sh
- commented out running bwa for now '#/share/apps/genomics/bwa-0.7.17/bwa', not sure if it's the right thing to replace 'module load...' with

#### 03a_variant_calling_bluecp.sh and 03a_call_SNVs_bluecp3.sh in wrapper
- edited paths of module versions to be loaded
- changed the path of 'CALLER' to UCL server tools
- commented out path to samtools as not sure if it replaces SAMTOOLS='samtools'
- edited UCL server options 
- trying to run in home directory:
ERROR: You must specify an input directory
ERROR: You must specify an output directory
ERROR: You must specify a reference file

#### 03b_summary_variant_calling.sh and 03b_call_SNVs_summarize_bluecp3.sh
- edited definitions of modules to be loaded
- did not add VMEM as no other memory variable was present (maybe it's not needed)

Have now edited all pathways from pipeline that directly lead to tools & wrapper but there are still files in both tools & wrapper because wrapper can lead to tools:

#### 04_fst_scans_etc_SNP_based.sh in wrapper 
the path to tools script has already been changed

#### popgenstats.pl
- edited paths to external programs (don't know what estpEM software is or whether its current path is correct)
- permission denied to run

## Wednesday 07/07/2021
- Edited 01a_modern_cutadapt_filtering_trimming.sh by comparing to museum script, no errors using ./ in the home directory but cannot submit yet as we need to upload data 
- Edits made to 01b_concat_fastq_R1.sh include:
   
   - changed from running to the work directory of UoB to running to the home directory of UCL server pchuckle
    - edited the scheduler directives to $ instead of PBS
    - removed nodes scheduler directive and added -l vmem and -S /bin/bash 
    - changed from -j oe to -j y in scheduler directives
    - changed the definition of SPECIESDIR to the current pathway to the E3_Aphantopus_hyperantus directory
    - changed the array job submissions to {SGE_TASK_ID}, current form is sed "${SGE_TASK_ID}q;d" but may need to change to sed -n ${SGE_TASK_ID}'{p;q}' input.data form
    - NOTE: when running using ./ I get the error cat: //SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/01a_museum_cutadapt_reads/: No such file or directory , so maybe it's because the data hasn't been uploaded yet or is it because 01a_museum_cutadapt_reads is in a different folder

- Edits made to 01b_concat_fastq_R2.sh include:
  
  - same as above for 01b_concat_fastq_R1.sh
  - NOTE: syntax error looking for missing " , but it is written in the same way as R1 which doesn't have this syntax error ??

### Thursday 08/07/2021
- Edits made to 01b_concat_fastq_R1.sh include:
   
   - change vmem to h_vmem for compatibility to UCL server
   - changed the pathway for 01a_museum_catadapt_reads to local/... but error of No such file or directory still for 01a_museum2_catadapt_reads
   - changed walltime to h_rt in the scheduler descriptive
   - qsub the job after deleting path for 01a_museum2_catadapt_reads which does create 33 sample files but no log files 
   - changed work directory to SGE_O_WORKDIR

- Edits made to 01b_concat_fastq_R2.sh include:

   - changed walltime to h_rt in scheduler descriptive
   - changed pathway for 01a_museum_catadapt_reads and 01a_museum2_catadapt_read to local/...
   - NOTE: still the same syntax error but cannot identify where as it is the exact same as 01b_concat_fastq_R1.sh but with 1 changed to 2
   - changed work directory to SGE_O_WORKDIR

- Edits made to 01c_bbtools_repair_museum_ARRAY.sh include:

   - changed scheduler descriptives from #PBS to #$
   - deleting nodes descriptive and added h_vmem descriptive
   - change mem to tmem, oe to y and walltime to h_rt
   - add /bin/bash descriptive
   - changed work directory to SGE_O_WORKDIR
   - define path to bbrepair through bbmap
   - changed path to load java
   - NOTE: errors of invalid maximum heap size and could not create the Java Virtual Machine
   - changed maximum heap size using export _JAVA_OPTIONS="-Xmx2g but still the same error occurs of invalid max heap size

- Edits made to 01d_bbtools_merge_museum_ARRAY.sh include:

   - same as above for 01c_bbtools_repair_museum_ARRAY.sh for all the scheduler descriptives
   - edited array jobs to SGE_TASK_ID
   - edited path to bbmerge
   - NOTE: can't read R1.museum.names.repaired: No such file or directory

- Edits made to 02a_MapwithBWAmem.ARRAY.sh include:

   - change from #PBS to #$
   - delete nodes descriptive and add h_vmem descriptive
   - change mem, walltime and oe to tmem, h_rt and y
   - add -S /bin/bash descriptive
   - change working directory
   - load bwa module 
   - change array jobs to SGE_TASK_ID
   - NOTE: error 'demultiplexed/': No such file or directory

- [Created a pipeline script for museum2 samples to use cutadapt](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/01a_museum2_cutadapt_filtering_trimming.sh) 
- [Created a pipeline script for modern_exp to use cutadapt](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/01a_modern_exp_cutadapt_filtering_trimming.sh)
- [Created a pipeline script for museum 2 sample to use fastqc](https://github.com/alexjvr1/VelocityUCL/blob/main/KaitlynNotes/00_fastqc_raw_museum2.sh)
- Successfully ran concat script (I think....!)

### Friday 09/07/2021
- Edits made to 03a_variant_calling_bluecp.sh include:

   - renamed file to 03a_variant_calling_UCL.sh
   - edited path to samtools
   - changed definition of working directory
   - NOTE: will need to change input and output files 

- Changed all thw permissions in VelocityPipeline folders so that everyone can execute all files
- Edits made to 03a_call_SNVs_bluecp.sh include:

   - renamed file to 03a_call_SNVs_UCL.sh
   - edited UCL server options such that

echo '#$ -S /bin/bash' >> $SMSJOB\
echo '#$ -N '$JOBNAME'' >> $SMSJOB\
echo '#$ -l tmem='$MEM'G' >> $SMSJOB\
echo '#$ -l h_vmem='$VMEM'G' >> $SMSJOB\
echo '#$ -l h_rt='$HRS'1:00:00' >> $SMSJOB\
echo '#$ -j y' >> $SMSJOB #concatenates error and output files (with prefix job1)\
echo '#$ -t 1-'$N >> $SMSJOB\
echo '#$ -o '$LOG >> $SMSJOB

- Edits made to 01c_bbtools_repair_museum_ARRAY.sh include:

   - added a PATH variable for a newly created folder containing the concatenated files 
   - changed the infiles to include $PATH
   - created a new folder in my local directory called 01c_musPERepaired which will be where the outfiles of 01c_bbtools_repair_museum_ARRAY.sh will go
   - NOTE: outfile appears empty after 'qsub'bing the script 

### Monday 12/07/2021

- Edits made to 01c_bbtools_repair_museum_ARRAY.sh include:

   - moved R2 concat files to 01a,mus.concat.cutadapt folder
   - coded for the creation of sample files
   - NOTE: now successfully running

- Edits made to 01d_bbtools_merge_museum_ARRAY.sh include:

   - added variables for the infile and outfile paths
   - created a new folder in local for the outfile
   - coded for the creation of sample names files
   - edited the PREFIX for the outfiles
  
- Edits made to 02a_MapwithBWAmem.ARRAY_museum.sh

   - changed job name and walltime in scheduler descriptive
   - edited paths to the reference, and input and output directories
   - coded for the creation of sample file names to map
   - hashed out NAME2 and deleted R2 paths 
   - NOTE: takes a long time to run. UPDATE: finished running with a total size of 7.3G

- Running scripts through from the start and transferring files across to shared folder:

   - Succesfully ran 00_fastqc_raw_museum.sh and 00_fastqc_raw_museum2.sh
   - succesfully ran 01a_museum_cutadapt_filtering_trimming.sh and 01a_museum2_cutadapt_filtering_trimming.sh
   - successfully ran 01b_concat_fastq_R1.sh and 01b_concat_fastq_R2.sh but manually move both sets of outfiles into a folder called 01a.mus.concat.cutadapt
   - 01c_bbtools_repair_museum_ARRAY.sh does not seem to be running after submitting the job; UPDATE: it simply takes quite a while to run
   - succesfully ran 01c_bbtools_repair_museum_ARRAY.sh
