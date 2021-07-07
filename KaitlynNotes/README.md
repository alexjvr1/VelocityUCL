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
    -NOTE: when running using ./ I get the error cat: //SAN/ugi/LepGenomics/E3_Aphantopus_hyperantus/01a_museum_cutadapt_reads/: No such file or directory , so maybe it's because the data hasn't been uploaded yet or is it because 01a_museum_cutadapt_reads is in a different folder

- Edits made to 01b_concat_fastq_R2.sh include:
- 
