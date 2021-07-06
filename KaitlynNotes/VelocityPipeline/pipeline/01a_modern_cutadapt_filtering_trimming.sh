
#!/bin/bash
##################################
# Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified 06/07/2021 11:22
##################################

# Creates submission script to use cutadapt to remove adapters from demultiplexed Illumina libraries.
# Poly-A and -T filters have been removed as we decided that these would be dropped due to mapping qual$
# Additional filters included in Lymantria monacha dataset
# -fwad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -fwad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
# -rvad2 CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -rvad3 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA



/SAN/ugi/LepGenomics/VelocityPipeline/tools/01a_parallel_cutadapt_UCL.sh \
-i /panfs/panasas01/bisc/aj18951/Pararge_aegeria/PA-01-2016 \
-o /panfs/panasas01/bisc/aj18951/Pararge_aegeria/PA-01-2016/01a_filtered_reads_modern -n 1 -t 8 -m 16 -$
-fwad1 AGATCGGAAGAGCACACGTCTGAACTCCAGTC -rvad1 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-minl 20 -phredq 20;

