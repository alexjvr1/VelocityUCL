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
- changed the path of CALLER to UCL server tools
- commented out path to samtools as not sure if it replaces SAMTOOLS='samtools'
- edited UCL server options 
- trying to run in home directory:
ERROR: You must specify an input directory
ERROR: You must specify an output directory
ERROR: You must specify a reference file
