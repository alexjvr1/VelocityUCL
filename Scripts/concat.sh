#!/bin/bash
###########################################
# (c) Alexandra Jansen van Rensburg
# last modified 22/11/2021 
###########################################

## Description: Given two input directories with raw fastq files (gzip allowed)
## concatenate all museum samples that appear in both folders. 
## The script is written to concatenate 33 museum samples that were sequenced twice for each species. 

VERSION='1-2021.11.22'

# Default values for optional variables
EXTRA=''
EMAIL=''
NCORES=1
HRS=1
MEM=16
QUEUE=CS
JOBNAME=Mus.concat
SHARED=/SAN/ugi/LepGenomics/
INDIR1=$SPECIESDIR/00_raw_reads_museum/ALLSAMPLES/
INDIR2=$SPECIESDIR/00_raw_reads_museum2/ALLSAMPLES
OUTDIR=$SPECIESDIR/00_raw_reads_museum_FINAL


function author {
	echo
	echo "#########################################"
	echo "  $(basename $0)"
	echo "  version $VERSION"
	echo "  (c) Alexandra Jansen van Rensburg"
	echo "  alexjvr@gmail.com"
	echo "#########################################"
	echo
}

function usage {
	echo
	echo "Usage:"
	echo "  $(basename $0)"
	echo "      -I <First input directory> => (optional, default=$INDIR1). Folder with first batch of sequencing (museum1) with fastq.gz files. Paired files shall have the same name with following extension: *R1.fastq* - *R2.fastq*"
	echo "      -i <Second input directory> => (optional, default=$INDIR2). Folder with second batch of sequencing (museum2) with fastq.gz files. Paired files shall have the same name with following extension: *R1.fastq* - *R2.fastq*"
	echo "      -o <output directory> => (optional, default=$OUTDIR). Output folder to save concatenated files"
	echo "      -S <species directory> => Name of species directory"
	echo "      -n <number of processors> => processors per trimming (optional, default=$NCORES)"
	echo "      -t <allocated time> => Allocated time (in hours) for each analysis (optional: default=$HRS)"
	echo "      -m <allocated memory> => Allocated memory (in gigabytes) for each analysis (optional: default=$MEM)"
	echo "      -q <queue> => SGE queue: uclcs | iceberg | popgenom | molecol (optional: default=$QUEUE)"
	echo "      -h => show this help"
	echo "      -N <JOBNAME> => Name given to job in queue. (optional: default=$JOBNAME)"
	echo ""
	echo "  Example:"
	echo "      $(basename $0) -S E3_Aphantopus_hyperantus -N E3.Mus.concat"
	echo ""
	echo ""
	exit 0
}

author

if [ "$#" -ge "1" ]; # min 1 args: 1 for -S <species directory>
then 
	while [ $# -gt 0 ]; do
		case "$1" in
			-h|-help) usage
			  	;;
      			-I) shift
			        INDIR1=$(readlink -f $1)
          			;;
		        -i)	shift
				INDIR2=$(readlink -f $1)
				;;
			-o)	shift
				OUTDIR=$(readlink -f $1)
				;;
			-S)	shift
				SPECIESDIR=$(readlink -f $1)
				;;
			-n)	shift
				NCORES=$1
				;;
			-t)	shift
				HRS=$1
				;;
			-m)	shift
				MEM=$1
				;;
			-q)	shift
				QUEUE=$1
				;;
			-N)	shift
				JOBNAME=$1
				;;
			*)	echo 
				echo "ERROR - Invalid option: $1"
				echo
				usage
				;;
		esac
		shift
	done
else
	usage
fi


N1=$(find $INDIR1 -maxdepth 1 -name "*R1.fastq*" | wc -l | cut -f1 -d" ")
N2=$(find $INDIR2 -maxdepth 1 -name "*R1.fastq*" | wc -l | cut -f1 -d" ")

TIMESTAMP=$(date +%Y%m%d-%H%M%S)

SMSJOB="$OUTDIR/Concat.$JOBNAME.$TIMESTAMP.smsjob.sh"
LOG="$OUTDIR/$JOBNAME.$TIMESTAMP.smsjob.log"



# Initialize submission script
mkdir -p $OUTDIR
echo '#!/bin/bash' > $SMSJOB

# UCL Server options
# -----------------------------------------
echo '#$ -S /bin/bash' >> $SMSJOB
echo '#$ -N '$JOBNAME >> $SMSJOB
echo '#$ -l h_rt='$HRS':00:00' >> $SMSJOB
echo '#$ -l tmem='$MEM'G' >> $SMSJOB
echo '#$ -l h_vmem='$MEM'G' >> $SMSJOB
if [ ! -z "$EMAIL" ];
then
	echo '#$ -m bea' >> $SMSJOB
	echo '#$ -M '$EMAIL >> $SMSJOB
fi

if (($NCORES > 1 ));
then
	echo '#$ -pe openmp '$NCORES >> $SMSJOB
fi

echo '#$ -t 1-'$N2 >> $SMSJOB
echo '#$ -j y' >> $SMSJOB
echo '#$ -o '$LOG >> $SMSJOB
# -----------------------------------------

cat >> $SMSJOB <<EOF


INDIR1=$INDIR1
INDIR2=$INDIR2
OUTDIR=$OUTDIR
FQFILES1=(\$INDIR2/*R1*.fastq*)
FQFILES2=(\$INDIR2/*R2*.fastq*)
INDEX=\$((SGE_TASK_ID-1))
FQ1=\${FQFILES1[\$INDEX]}
FQ2=\${FQFILES2[\$INDEX]}


FQ1=\$(basename \$FQ1)
FQ2=\$(basename \$FQ2)

LOG="\$OUTDIR/"\${FQ1%%_R1.*}".log"


echo "Concatenating $INDIR1\$FQ1 AND $INDIR2\$FQ1 files..." > \$LOG
echo >> \$LOG 
echo "CMD: " >> \$LOG
echo "cat \\\">> \$LOG
echo "$INDIR1\$FQ1 $INDIR2\$FQ1 \\\">> \$LOG
echo ">$OUTDIR\${FQ1%%.fastq*}_paired.fastq.gz \${FQ1%%.fastq*}_unpaired.fastq.gz \\\">> \$LOG
echo "\${FQ2%%.fastq*}_paired.fastq.gz \${FQ2%%.fastq*}_unpaired.fastq.gz \\\">> \$LOG
echo >> \$LOG
echo "---------------------------------------------------" >> \$LOG
echo >> \$LOG

cd \$SPECIESDIR


##Concat fastq files

echo "[concatenating] $1 and $2" >> concat.mus.log
printf "\n"

echo "while read NAME1 <&MUS1 && read NAME2 <&MUS2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL1 $SPECIESDIR/$PATH2/$NAME2$TAIL1 > $OUTPUT/$NAME1_R1.concat.fastq.gz; done" >> concat.mus.log
time while read NAME1 <&1 && read NAME2 <&2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL1 $SPECIESDIR/$PATH2/$NAME2$TAIL1 > $OUTPUT/$NAME1_R1.concat.fastq.gz; done 1<museum1.names 2<museum2.names

echo "while read NAME1 <&MUS1 && read NAME2 <&MUS2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL2 $SPECIESDIR/$PATH2/$NAME2$TAIL2 > $OUTPUT/$NAME1_R2.concat.fastq.gz; done" >> concat.mus.log
time while read NAME1 <&1 && read NAME2 <&2; do cat $SPECIESDIR/$PATH1/$NAME1$TAIL2 $SPECIESDIR/$PATH2/$NAME2$TAIL2 > $OUTPUT/$NAME1_R2.concat.fastq.gz; done 1<museum1.names 2<museum2.names


EOF


chmod +x $SMSJOB


echo "Command to submit the job to UCL ($QUEUE queue):"
echo
echo "qsub $SMSJOB"
echo
