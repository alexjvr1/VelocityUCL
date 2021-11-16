#!/bin/bash

# (c) Romain Villoutreix
# romain.villoutreix@gmail.com
# Last modified: 10/10/2016 
#############################################
# Modified by Alexandra Jansen van Rensburg
#
# Last modified 15/11/2021
#############################################

# Description:
# Given an input directory with raw fastq files (gzip compression allowed)
# it will run Trimmomatic in parallel

# Changelog
# !!!! The script is only designed for paired-end reads at the moment !!! try something else at your own risk !!! 

VERSION='2-2021.11.15'
CMD="/SAN/ugi/LepGenomics/Software/Trimmomatic-0.39/trimmomatic-0.39.jar"

# Default values for optional variables
EXTRA=''
EMAIL=''
NCORES=1
HRS=8
MEM=16
ADAPTPATH='/SAN/ugi/LepGenomics/Software/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa'
ILLUMINACLIP='2:30:8:1:true'
LEADING=20
TRAILING=20
SLIDINGWINDOW='4:20'
MINLEN=20
PHREDSCORE=33
CROP=150
HEADCROP=0
AVGQUAL=20
QUEUE=uclcs

function author {
	echo
	echo "#########################################"
	echo "  $(basename $0)"
	echo "  version $VERSION"
	echo "  (c) Romain Villoutreix"
	echo "  romain.villoutreix@gmail.com"
	echo "#########################################"
	echo
}

function usage {
	echo
	echo "Usage:"
	echo "  $(basename $0)"
	echo "		!!! trimmomatic is processing task as order of submission. The help order is submission task order here !!!"
	echo "      -i <input directory> => Folder with fastq.gz files !! only work for paired-end reads at the moment !! Paired files shall have the same name with following extension: *R1.fastq* - *R2.fastq*"
	echo "      -o <output directory> => Output folder to save trimmed files"
	echo "      -n <number of processors> => processors per trimming (optional, default=$NCORES)"
	echo "      -t <allocated time> => Allocated time (in hours) for each analysis (optional: default=$HRS)"
	echo "      -m <allocated memory> => Allocated memory (in gigabytes) for each analysis (optional: default=$MEM)"
	echo "      -ph <phred score> => Phred score encoding (either 33 or 64 ; optional: default=phred$PHREDSCORE)"
	echo "      -c <CROP option> => see trimmomatic manual (optional: default=$CROP, but put your reads length here)"
	echo "      -hc <HEADCROP option> => see trimmomatic manual (optional: default=$HEADCROP)"
	echo "      -ad <adaptor sequences> => fasta files with adaptors sequences (optional: default=homemade Illumina TruSeq DNA methylation adaptors)"
	echo "      -illclip <ILLUMINACLIP options> => see trimmomatic manual (optional: $ILLUMINACLIP)"
	echo "      -le <LEADING options> => see trimmomatic manual (optional: default=$LEADING)"
	echo "      -tr <TRAILING options> => see trimmomatic manual (optional: default=$TRAILING)"
	echo "      -sw <SLIDINGWINDOW options> => see trimmomatic manual (optional: default=$SLIDINGWINDOW)"
	echo "      -minl <MINLEN options> => see trimmomatic manual (optional: default=$MINLEN)"
	echo "      -avg <AVGQUAL options> => see trimmomatic manual (optional: default=$AVGQUAL)"
	echo "      -q <queue> => SGE queue: uclcs | iceberg | popgenom | molecol (optional: default=$QUEUE)"
	echo "      -e <email> => Notification email address (default=$EMAIL)"
	echo "      -h => show this help"
	echo ""
	echo "  Example:"
	echo "      $(basename $0) -i raw_reads -o trimmed_reads -q uclcs"
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
		      -avg)	shift
				AVGQUAL=$1
				;;
			-t)	shift
				HRS=$1
				;;
			-m)	shift
				MEM=$1
				;;
			-ph) 	shift
				PHREDSCORE=$1
				;;
			-ad) 	shift
				ADAPTPATH=$(readlink -f $1)
				;;
			-illclip) shift
				ILLUMINACLIP=$1
				;;
			-le) 	shift
				LEADING=$1
				;;
			-tr) 	shift
				TRAILING=$1
				;;
			-sw) 	shift
				SLIDINGWINDOW=$1
				;;
			-minl) shift
				MINLEN=$1
				;;
			-avg) shift
				AVGQUAL=$1
				;;
			-c) shift
				CROP=$1
				;;
			-hc) shift
				HEADCROP=$1
				;;
			-q)	shift
				QUEUE=$1
				;;
			-e)	shift
				EMAIL=$1
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

N=$(find $INDIR -maxdepth 1 -name "*R1.fastq*" | wc -l | cut -f1 -d" ")

TIMESTAMP=$(date +%Y%m%d-%H%M%S)
TMPDIR="parallel_trimmomatic_"$TIMESTAMP"_"$RANDOM

INPUT_LOCALDIR="$TMPDIR""_in"

SMSJOB="$OUTDIR/parallel_trimmomatic.$TIMESTAMP.smsjob.sh"
LOG="$OUTDIR/parallel_trimmomatic.$TIMESTAMP.smsjob.log"

# Initialize submission script
mkdir -p $OUTDIR
echo '#!/bin/bash' > $SMSJOB

# UCL Server options
# -----------------------------------------
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
FQFILES1=(\$INDIR/*R1.fastq*)
FQFILES2=(\$INDIR/*R2.fastq*)
INDEX=\$((SGE_TASK_ID-1))
FQ1=\${FQFILES1[\$INDEX]}
FQ2=\${FQFILES2[\$INDEX]}

INPUT_TMPDIR=$INPUT_LOCALDIR

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

echo "trimmomatic \$FQ1 \$FQ2 files..." > \$LOG
echo >> \$LOG 
echo "CMD: " >> \$LOG
echo "$CMD PE -phred$PHREDSCORE \\\">> \$LOG
echo "\$FQ1 \$FQ2 \\\">> \$LOG
echo "\${FQ1%%.fastq*}_paired.fastq.gz \${FQ1%%.fastq*}_unpaired.fastq.gz \\\">> \$LOG
echo "\${FQ2%%.fastq*}_paired.fastq.gz \${FQ2%%.fastq*}_unpaired.fastq.gz \\\">> \$LOG
echo "CROP:$CROP\\\">> \$LOG
echo "HEADCROP:$HEADCROP\\\">> \$LOG
echo "ILLUMINACLIP:$ADAPTPATH:$ILLUMINACLIP \\\">> \$LOG
echo "LEADING:$LEADING \\\">> \$LOG
echo "TRAILING:$TRAILING\\\">> \$LOG
echo "SLIDINGWINDOW:$SLIDINGWINDOW\\\">> \$LOG
echo "MINLEN:$MINLEN\\\">> \$LOG
echo "AVGQUAL:$AVGQUAL\\\">> \$LOG
echo >> \$LOG
echo "---------------------------------------------------" >> \$LOG
echo >> \$LOG

cd \$INPUT_TMPDIR

$CMD PE -phred$PHREDSCORE \$FQ1 \$FQ2 \${FQ1%%.fastq*}_paired.fastq.gz \${FQ1%%.fastq*}_unpaired.fastq.gz \${FQ2%%.fastq*}_paired.fastq.gz \${FQ2%%.fastq*}_unpaired.fastq.gz CROP:$CROP HEADCROP:$HEADCROP ILLUMINACLIP:$ADAPTPATH:$ILLUMINACLIP LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$SLIDINGWINDOW MINLEN:$MINLEN AVGQUAL:$AVGQUAL>> \$LOG 2>&1

# Copy results back to output directory
# and check all files are good
PAIRED1=\${FQ1%%.fastq*}_paired.fastq.gz;
PAIRED2=\${FQ2%%.fastq*}_paired.fastq.gz;
UNPAIRED1=\${FQ1%%.fastq*}_unpaired.fastq.gz
UNPAIRED2=\${FQ2%%.fastq*}_unpaired.fastq.gz

#paired 1 
while [[ ! -e \$INPUT_TMPDIR/\$PAIRED1 || ! -e \$OUTDIR/\$PAIRED1 || \
     "\$(md5sum \$INPUT_TMPDIR/\$PAIRED1 | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$PAIRED1 | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$PAIRED1 \$OUTDIR/ >& /dev/null
done

#paired 2 
 while [[ ! -e \$INPUT_TMPDIR/\$PAIRED2 || ! -e \$OUTDIR/\$PAIRED2 || \
      "\$(md5sum \$INPUT_TMPDIR/\$PAIRED2 | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$PAIRED2 | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$PAIRED2 \$OUTDIR/ >& /dev/null
done

#unpaired 1 
while [[ ! -e \$INPUT_TMPDIR/\$UNPAIRED1 || ! -e \$OUTDIR/\$UNPAIRED1 || \
      "\$(md5sum \$INPUT_TMPDIR/\$UNPAIRED1 | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$UNPAIRED1 | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$UNPAIRED1 \$OUTDIR/ >& /dev/null
done

#unpaired 2 
while [[ ! -e \$INPUT_TMPDIR/\$UNPAIRED2 || ! -e \$OUTDIR/\$UNPAIRED2 || \
     "\$(md5sum \$INPUT_TMPDIR/\$UNPAIRED2 | awk '{print \$1}')" != "\$(md5sum \$OUTDIR/\$UNPAIRED2 | awk '{print \$1}')" ]];
do
     sleep \$((30+\$RANDOM%90)) # wait 30-120 secs to try again
     cp \$INPUT_TMPDIR/\$UNPAIRED2 \$OUTDIR/ >& /dev/null
done

#remove temporary files
rm -f \$INPUT_TMPDIR/\$PAIRED1 \$INPUT_TMPDIR/\$PAIRED2 \$INPUT_TMPDIR/\$UNPAIRED1 \$INPUT_TMPDIR/\$UNPAIRED2
EOF


chmod +x $SMSJOB


echo "Command to submit the job to UCL ($QUEUE queue):"
echo
echo "qsub $SMSJOB"
echo
