#!/bin/bash
# v1 - modified path to tools script 
# Last modified 05/07/2021 16:23

#PBS -N 03b_Summary  ##job name
#PBS -l nodes=1:ppn=1  #nr of nodes and processors per node
#PBS -l mem=16gb #RAM
#PBS -l walltime=8:00:00 ##wall time.  
#PBS -j oe  #concatenates error and output files (with prefix job1)


#run job in working directory
cd $PBS_O_WORKDIR


/SAN/ugi/LepGenomics/VelocityPipeline/tools/popgenstats.pl -i 04a_Fst_scans_and_other_scatistics_indir/G1.subset.museum.modern.bcf -p 04a_Fst_scans_and_other_scatistics_indir/Thymelicus_acteon_population_file.dsv -o 04b_Fst_scans_and_other_scatistics_outdir/Fst_outfile.dsv -n 1 -ms 0.5 -msp 0.4 -mf 0 -pr 1 -am 1 -fm 0;
