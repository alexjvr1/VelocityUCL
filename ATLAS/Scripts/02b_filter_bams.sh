#Filter mapped bams for quality, proper pairs, and primary mapped position (no secondary mapped pairs). 
###




#Export path to samtools and samtools libraries
export PATH=/share/apps/genomics/samtools-1.9/bin:$PATHexport 
LD_LIBRARY_PATH=/share/apps/genomics/samtools-1.9/lib:$LD_LIBRARY_PATH


##Define variables
SHAREDFOLDER=/SAN/ugi/LepGenomics
SPECIES=E3_Aphantopus_hyperantus
INFOLDER=02a_mapped_modern
OUTFOLDER=02a_mapped_modern
INPUT=$SHAREDFOLDER/$SPECIES/$INFOLDER
OUTPUT=$SHAREDFOLDER/$SPECIES/$OUTFOLDER
TAIL1=
TAIL2="flt.bam"

##Samtools flags can be found here: https://broadinstitute.github.io/picard/explain-flags.html
##2028: read unmapped, mate unmapped, not primary alignment, read fails platform/vendor quality checks; supplementary reads
##This should leave only properly paired mapped reads with a primary alignment only. 

while read NAME in <&1; do samtools view -b -f 2 -F 2028 -mapq>20 $INPUT/$NAME$TAIL1 > $OUTPUT/$NAME.$TAIL2; done 1<modc.names

for i in $(ls $INPUT/*flt.bam); do samtools flagstat $i >> $INPUT/flagstat.log; done
