#!/bin/bash

# (c) Romain Villoutreix
# romain.villoutreix@gmail.com
# Last modified: 10/10/2016 
# Modified by Alexandra Jansen van Rensburg
# Date: 24/10/2018 

# Description:
# Given an input directory with raw fastq files (gzip compression allowed)
# it will run cutadapt in parallel

# Changelog
# !!!! The script is only designed for paired-end reads at the moment !!! try something else at your own risk !!! 
# 1.0.1 - AJvR 24/10/2018
# 	- Change script to make compatible with blue crystal p3 server
# 1.0.2 - AJvR 5/2/2019
#	- Remove references to scratch directory. Tmp directory now made in local working directory. 
#	- write outdir to local dir


VERSION='1.0.1-2018.10.24'
CMD="cutadapt"

# Default values for optional variables
EXTRA=''
EMAIL=''
NCORES=1
HRS=8
MEM=16
#FWADAPT1='AGATCGGAAGAGCACACGTCTGAACTCCAGTC';
#FWADAPT2='CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT';
#FWADAPT3='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';
#RVADAPT1='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT';
#RVADAPT2='CTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT';
#RVADAPT3='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA';
MINLEN=35;
PHREDSCORE=33;
PHREDQUAL=20;
JOBNAME='Cutadapt';


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
	echo "		!!! cutadapt is processing task in a given order. The cutadapt help order is task order here !!! Not all options have been implemented here, but modifying this script should not be too painful !!!"
	echo "      -i <input directory> => Folder with fastq.gz files !! only work for paired-end reads at the moment !! Paired files shall have the same name with following extension: *R1.fastq* - *R2.fastq*"
	echo "      -o <output directory> => Output folder to save trimmed files"
	echo "      -n <number of processors> => processors per trimming (optional, default=$NCORES)"
	echo "      -t <allocated time> => Allocated time (in hours) for each analysis (optional: default=$HRS)"
	echo "      -m <allocated memory> => Allocated memory (in gigabytes) for each analysis (optional: default=$MEM)"
	echo "      -ph <phred score> => Phred score encoding (either 33 or 64 ; optional: default=phred$PHREDSCORE)"
	echo "      -fwad1 <1st adaptor sequence to look in forward read/R1> (optional: default=$FWADAPT1)"
	echo "      -fwad2 <2nd adaptor sequence to look in forward read/R1> (optional: default=$FWADAPT2)"
	echo "      -fwad3 <3rd adaptor sequence to look in forward read/R1> (optional: default=$FWADAPT3)"
	echo "      -rvad1 <1st adaptor sequence to look in reverse read/R2> (optional: default=$RVADAPT1)"
	echo "      -rvad2 <2nd adaptor sequence to look in reverse read/R2> (optional: default=$RVADAPT2)"
	echo "      -rvad3 <3rd adaptor sequence to look in reverse read/R2> (optional: default=$RVADAPT3)"
	echo "      -minl <--minimum-length option> => see cutadapt manual (optional: default=$MINLEN)"
	echo "      -phredq <Quality trimming -q option> => see cutadapt manual (optional: default=$PHREDQUAL)"
        echo "      -job <name of job>    => Name of current job"
	echo "      -q <queue> => SGE queue: iceberg | popgenom | molecol (optional: default=$QUEUE)"
	echo "      -e <email> => Notification email address (default=$EMAIL)"
	echo "      -h => show this help"
	echo ""
	echo "  Example:"
	echo "      $(basename $0) -i raw_reads -o trimmed_reads -q popgenom"
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
			-n)shift
				NCORES=$1
				;;
			-t)	shift
				HRS=$1
				;;
			-m)	shift
				MEM=$1
				;;
			-ph) shift
				PHREDSCORE=$1
				;;
			-fwad1) shift
				FWADAPT1=$1
				;;
			-fwad2) shift
				FWADAPT2=$1
                ;;
			-fwad3) shift
				FWADAPT3=$1
				;;
			-rvad1) shift
				RVADAPT1=$1
				;;
			-rvad2) shift
				RVADAPT2=$1
				;;
			-rvad3) shift
				RVADAPT3=$1
				;;
			-minl) shift
				MINLEN=$1
				;;
			-phredq) shift
				PHREDQUAL=$1
				;;
                        -job)   shift
                                JOBNAME=$1
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

N=$(find $INDIR -maxdepth 1 -name "*R1*.fastq*" | wc -l | cut -f1 -d" ")

TIMESTAMP=$(date +%Y%m%d-%H%M%S)
TMPDIR="parallel_cutadapt_"$TIMESTAMP"_"$RANDOM

INPUT_LOCALDIR=$INDIR/"$TMPDIR""_in"

#if [[ $QUEUE == "popgenom" || $QUEUE == "molecol" ]];
#then
#	INPUT_LOCALDIR="/local/tmp/$TMPDIR""_in"
#fi

SMSJOB="$OUTDIR/parallel_cutadapt.$TIMESTAMP.smsjob.sh"
LOG="$OUTDIR/parallel_cutadapt.$TIMESTAMP.smsjob.log"

# Initialize submission script
mkdir -p $OUTDIR
echo '#!/bin/bash' > $SMSJOB

echo ""
echo "Load necessary modules"
echo "module load languages/python-2.7.10"
echo "cutadapt installed locally"
echo ""

# PBS OPTIONS
# -----------------------------------------
echo '#PBS -N '$JOBNAME'' >> $SMSJOB  
echo '#PBS -l nodes=1:ppn=1' >> $SMSJOB 
echo '#PBS -l mem='$MEM'gb' >> $SMSJOB 
echo '#PBS -l walltime='$HRS':00:00' >> $SMSJOB
echo '#PBS -j oe' >> $SMSJOB #concatenates error and output files (with prefix job1)
echo '#PBS -t 1-'$N >> $SMSJOB
echo '#PBS -o '$LOG >> $SMSJOB
# -----------------------------------------

echo ""


echo '# Run in local directory' >> $SMSJOB
echo 'cd $PBS_O_WORKDIR' >> $SMSJOB
echo ''
echo ''

cat >> $SMSJOB <<EOF

INDIR=$INDIR
OUTDIR=$OUTDIR
FQFILES1=(\$INDIR/*R1*.fastq*)
FQFILES2=(\$INDIR/*R2*.fastq*)
INDEX=\$((PBS_ARRAYID-1))
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

echo "cutadapt \$FQ1 \$FQ2 files..." > \$LOG
echo >> \$LOG 
echo "CMD: " >> \$LOG
echo "$CMD  -a $FWADAPT1 -a $FWADAPT2 -a $FWADAPT3 -A $RVADAPT1 -A $RVADAPT2 -A $RVADAPT3 \\\">> \$LOG
echo "-m $MINLEN -q $PHREDQUAL \\\">> \$LOG

EOF

if [[ $PHREDSCORE == 64 ]];
 then
     echo 'echo "--quality-base=64 \\\">> \$LOG' >> $SMSJOB;
fi

cat >> $SMSJOB <<EOF

echo "-o \${FQ1%%R1.fastq*}cutadapt_filtered_R1.fastq.gz -p \${FQ2%%R2.fastq*}cutadapt_filtered_R2.fastq.gz \\\">> \$LOG
echo "\$FQ1 \$FQ2 \\\">> \$LOG
echo >> \$LOG
echo "---------------------------------------------------" >> \$LOG
echo >> \$LOG

cd \$INPUT_TMPDIR

EOF

if [[ $PHREDSCORE == 33 ]];
	then 
		echo "$CMD -a $FWADAPT1 -a $FWADAPT2 -a $FWADAPT3 -A $RVADAPT1 -A $RVADAPT2 -A $RVADAPT3 -m $MINLEN -q $PHREDQUAL -o \${FQ1%%R1.fastq*}cutadapt_filtered_R1.fastq.gz -p \${FQ2%%R2.fastq*}cutadapt_filtered_R2.fastq.gz \$FQ1 \$FQ2 >> \$LOG 2>&1" >> $SMSJOB;
fi

 if [[ $PHREDSCORE == 64 ]];
     then 
         echo "$CMD -a $FWADAPT1 -a $FWADAPT2 -a $FWADAPT3 -A $RVADAPT1 -A $RVADAPT2 -A $RVADAPT3 -m $MINLEN -q $PHREDQUAL --quality-base=64 -o \${FQ1%%R1.fastq*}cutadapt_filtered_R1.fastq.gz -p \${FQ2%%R2.fastq*}cutadapt_filtered_R2.fastq.gz" >> $SMSJOB;
fi

cat >> $SMSJOB <<EOF

# Copy results back to output directory
# and check all files are good
PAIRED1=\${FQ1%%R1.fastq*}cutadapt_filtered_R1.fastq.gz;
PAIRED2=\${FQ2%%R2.fastq*}cutadapt_filtered_R2.fastq.gz;

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

#remove temporary files
rm -f \$INPUT_TMPDIR/\$PAIRED1 \$INPUT_TMPDIR/\$PAIRED2

EOF

chmod u+x $SMSJOB


echo "Command to submit the job to BlueCrystal p3:"
echo
echo "qsub $SMSJOB"
echo
