#!/bin/bash
#################################
# Alexandra Jansen van Rensburg
# Last modified: 24/10/2018
#################################

# Creates submission script for mapping reads to a reference genome using bwa-mem



/usr/local/extras/Genomics/scripts/rfseq_pipeline/05_map_reads.sh \
-i 01c_modern_mapping_input_folder/R1 \
-p 01c_modern_mapping_input_folder/R2 \
-u 01c_modern_mapping_input_folder/unpaired/ \
-o 02_modern_mapped_reads \
-r ../../reference_genomes/Lymantria_monacha_v1.0_minsize3kb.fasta \
-s 1 -n 1 -xi 0 -mi 1000 -t 50 -m 16;
