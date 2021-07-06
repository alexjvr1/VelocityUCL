# e-labbook

## editing scripts to suit UCL server
Things to change:

- path to wrapper or tools script
- change PBS to #$
- change array jobs to SGE_TASK_ID
- remove 'module load ..'
- VMEM added

#### 01a_modern_cutadapt_filtering_trimming.sh and 01a_parallel_cutadapt_UCL.sh
- made relevant changes with no syntax errors (using./) 
- unsure how to submit script using qsub to check it's correct
