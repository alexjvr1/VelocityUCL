#!/bin/bash

# (c) Alexandra Jansen van Rensburg
# alexjvr@gmail.com
# Last modified: 17/11/2021 
#############################################

# Description:
# Given an input directory with adapter trimmed fastq files (gzip compression allowed)
# it will run AdapterRemoval in parallel

# Changelog
# !!!! The script is only designed for paired-end reads at the moment !!! try something else at your own risk !!! 

VERSION='1-2021.11.17'
CMD="/SAN/ugi/LepGenomics/Software/adapterremoval-2.3.1/build/AdapterRemoval"


# Default values for optional variables
EXTRA=''
EMAIL=''
NCORES=1
HRS=1
MEM=8
QUEUE=CS
TRIMN=TRUE
TRIMQUAL=TRUE
COLLAPSE=''
OVERLAPLENGTH=11
MINLENGTH=15
MAXLENGTH=4294967295


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
	echo "      -i <input directory> => Folder with fastq* files !! only work for paired-end reads at the moment !! Paired files shall have the same name with following extension: *_R1_paired.fastq* - *_R2_paired.fastq*"
	echo "      -o <output directory> => Output folder to save trimmed files"
	echo "      -n <number of processors> => processors per trimming (optional, default=$NCORES)"
	echo "      -t <allocated time> => Allocated time (in hours) for each analysis (optional: default=$HRS)"
	echo "      -m <allocated memory> => Allocated memory (in gigabytes) for each analysis (optional: default=$MEM)"
	echo "      -q <queue> => SGE queue: uclcs | iceberg | popgenom | molecol (optional: default=$QUEUE)"
	echo "      -e <email> => Notification email address (default=$EMAIL)"
	echo "      -h => show this help"
	echo "      -N <JOBNAME> => Name given to job in queue"
	echo ""
	echo "  Example:"
	echo "      $(basename $0) -i raw_reads -o trimmed_reads -q CS"
	echo ""
	echo ""
	exit 0
}

author

if [ "$#" -ge "2" ]; # min 2 args: 1 for -i <input directory>, 1 for -o <output directory>
then 
	while [ $# -gt 0 ]; do
		case "$1" in
			-h|-help) usage
					  ;;
			-i)	shift
				INDIR=$(readlink -f $1)
				;;
			-o)	shift
				OUTDIR=$(readlink -f $1)
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
			-e)	shift
				EMAIL=$1
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

N=$(find $INDIR -maxdepth 1 -name "*_R1_paired.fastq*" | wc -l | cut -f1 -d" ")

TIMESTAMP=$(date +%Y%m%d-%H%M%S)
TMPDIR="parallel_adapterremoval_"$TIMESTAMP"_"$RANDOM

INPUT_LOCALDIR="$TMPDIR""_in"

SMSJOB="$OUTDIR/parallel_adapterremoval.$TIMESTAMP.smsjob.sh"
LOG="$OUTDIR/parallel_adapterremoval.$TIMESTAMP.smsjob.log"

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

echo '#$ -t 1-'$N >> $SMSJOB
echo '#$ -j y' >> $SMSJOB
echo '#$ -o '$LOG >> $SMSJOB
# -----------------------------------------

cat >> $SMSJOB <<EOF

INDIR=$INDIR
OUTDIR=$OUTDIR
FQFILES1=(\$INDIR/*_R1_paired.fastq*)
FQFILES2=(\$INDIR/*_R2_paired.fastq*)
INDEX=\$((SGE_TASK_ID-1))
FQ1=\${FQFILES1[\$INDEX]}
FQ2=\${FQFILES2[\$INDEX]}

INPUT_TMPDIR=$OUTDIR/$INPUT_LOCALDIR


# Create temporary local directories
if [ ! -e \$INPUT_TMPDIR ];
then
	mkdir -p \$INPUT_TMPDIR
fi

# copy to local temporary directory
cp \$FQ1 \$INPUT_TMPDIR/
cp \$FQ2 \$INPUT_TMPDIR/

FQ1=\$(basename \$FQ1)
FQ2=\$(basename \$FQ2)

LOG="\$OUTDIR/"\${FQ1%%_R1.*}".log"


echo "AdapterRemoval \$FQ1 \$FQ2 files..." > \$LOG
echo >> \$LOG 
echo "CMD: " >> \$LOG
echo "$CMD --basename \${FQ1%%_R1_paired.fastq*} \\\">> \$LOG
echo "--file1 \$FQ1 --file2 \$FQ2 \\\">> \$LOG
echo "--trimqualities\\\">> \$LOG
echo "--trimns\\\">> \$LOG
echo "--collapse\\\">> \$LOG
echo >> \$LOG
echo "---------------------------------------------------" >> \$LOG
echo >> \$LOG

cd \$INPUT_TMPDIR

time $CMD --basename \${FQ1%%_R1_paired.fastq*} --file1 \$FQ1 --file2 \$FQ2 --trimqualities --trimns --collapse >> \$LOG 2>&1

# Copy results back to output directory
# and check all files are good
COLLAPSED=\$FQ1.collapsed;
COLTRUNC=\$FQ1.collapsed.truncated
DISCARDED=\$FQ1.discarded
PAIR1TRUNC=\$FQ1.pair1.truncated
PAIR2TRUNC=\$FQ1.pair2.truncated
SINGLETON=\$FQ1.singleton.truncated
SETTINGS=\$FQ1.settings

#settings  
while [[ ! -e \$INPUT_TMPDIR/\$COLLAPSED || ! -e \$OUTDIR/\$COLLAPSED || \
     "\$(md5sum \$INPUT_TMPDIR/\$COLLAPSED | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$COLLAPSED | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$COLLAPSED \$OUTDIR/ >& /dev/null
done

#collapsed truncated  
while [[ ! -e \$INPUT_TMPDIR/\$COLTRUNC || ! -e \$OUTDIR/\$COLTRUNC || \
     "\$(md5sum \$INPUT_TMPDIR/\$COLTRUNC | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$COLTRUNC | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$COLTRUNC \$OUTDIR/ >& /dev/null
done

#paired 1 
while [[ ! -e \$INPUT_TMPDIR/\$PAIR1TRUNC || ! -e \$OUTDIR/\$PAIR1TRUNC || \
     "\$(md5sum \$INPUT_TMPDIR/\$PAIR1TRUNC | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$PAIR1TRUNC | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$PAIR1TRUNC \$OUTDIR/ >& /dev/null
done

#paired 2 
 while [[ ! -e \$INPUT_TMPDIR/\$PAIR2TRUNC || ! -e \$OUTDIR/\$PAIR2TRUNC || \
      "\$(md5sum \$INPUT_TMPDIR/\$PAIR2TRUNC | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$PAIR2TRUNC | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$PAIR2TRUNC \$OUTDIR/ >& /dev/null
done

#singleton 
while [[ ! -e \$INPUT_TMPDIR/\$SINGLETON || ! -e \$OUTDIR/\$SINGLETON || \
      "\$(md5sum \$INPUT_TMPDIR/\$SINGLETON | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$SINGLETON | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$SINGLETON \$OUTDIR/ >& /dev/null
done

#settings 
while [[ ! -e \$INPUT_TMPDIR/\$SETTINGS || ! -e \$OUTDIR/\$SETTINGS || \
     "\$(md5sum \$INPUT_TMPDIR/\$SETTINGS | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$SETTINGS | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$SETTINGS \$OUTDIR/ >& /dev/null
done

#remove temporary files
rm -f \$INPUT_TMPDIR/\$PAIRED1 \$INPUT_TMPDIR/\$PAIRED2 \$INPUT_TMPDIR/\$UNPAIRED1 \$INPUT_TMPDIR/\$UNPAIRED2
EOF


chmod +x $SMSJOB


echo "Command to submit the job to UCL ($QUEUE queue):"
echo
echo "qsub $SMSJOB"
echo
