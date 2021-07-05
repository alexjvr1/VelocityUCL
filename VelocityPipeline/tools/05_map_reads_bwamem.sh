#!/bin/bash

# (c) Victor Soria-Carrasco
# victor.soria.carrasco@gmail.com
# Last modified: 25/04/2018 12:33:01

# Description:
# Map multiple fastq file using BWA-mem

# Changelog
#	1.2 - 08/03/2016
#	  - Updated to use bowtie 2.2.7
#	  - Updated to use samtools 1.3

#	1.3 - 05/07/2016
#	  - Updated to use bowtie 2.2.9
#     - Updated to use samtools 1.3.1
#	  - Fixed bug when sorting big files, that caused potential concurrent 
#		writing to the same temporary file  when sorting the alignment of 
#		multiple samples (when exceeding the memory allowed for sorting)
#		memory for sorting has been set to the same memory allocated for
#		the run and temporary files are saved with sample name in the same
#		output folder to avoid collision with the temporary file of other samples

#	1.3.1 - 11/08/2016
#	  - Fixed issue when any directory on the path have dots

#	1.3.2 - 04/05/2017
#	  - Fixed issue with output files when input filenames have dots

#   1.3.3 - 02/08/2017
#     - Updated to use bowtie 2.3.2
#     - Updated to use samtools 1.5
#     - Reduced memory for sorting big files with samtools to half of
#       the requested memory
#     - Fixed bug related to creation of temporary directory used 
#       by samtools sort when run out of memory

#   1.4 - 22/08/2017
#     - Added option to map unpaired reads along with paired reads
#	1.4.1 - 30/08/2017
#	  - Fixed bug when using single reads
#	1.4.2 - 1/09/2017
#	  - Fixed bug when testing paired end reads without giving directory
#       of unpaired reads
#   1.4.3 - 10/10/2017
#     - Corrected log output; added samtools merge log when mapping
#       paired and unpaired reads; added option for minimum insert size
#   1.4.4 - 23/04/2018
#     - Updated to use bowtie 2.3.4.1
#     - Updated to use samtools 1.8
#     - Made phred33 the default QS encoding

# ToDo: add summarization scripts with number of reads mapped, etc. allow output in CRAM format

VERSION='1.4.4-2018.04.25'

# External programs path and options
module load apps/bwa-0.7.15
module load apps/samtools-1.8
module load languages/perl-5.14.2

MAPPER='bwa'
MAPPER_ARGS='mem'

MININS='0' # minimum fragment length (insert size)
MAXINS='500' # maximum fragment length (insert size)

SAMTOOLS='samtools'

# Default values for optional variables
ENC='1' # Quality score encoding
EMAIL=''
QUEUE='BlueCrystal'
HRS=24
MEM=16
NTHR=1
PAIREND=0
UNPAIRED=0
ERR=0
JOBNAME='bwamem'
# Show author
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------

# Show usage
# -----------------------------------------------------------------------------
function usage {
	echo
	echo "Usage:"
	echo "  $(basename $0)"
	echo "      -i <input directory> => directory with fastq format files (.fq.gz, .fastq.gz, .fq.bz2, and .fastq.bz2 accepted)"
	echo "      -p <input directory> => directory with paired fastq format files (.fq.gz, .fastq.gz, .fq.bz2, and .fastq.bz2 accepted) (optional)"
	echo "      -u <input directory> => directory with unpaired fastq format files (.fq.gz, .fastq.gz, .fq.bz2, and .fastq.bz2 accepted) (optional)"
	echo "      -o <output directory> => Output folder"
	echo "      -r <reference file> => Genome/reference in fasta format"
	echo "      -s <0 (phred64, Illumina 1.5) | 1 (phred33, Illumina 1.9)> => Phred quality score (optional: default=1)"
	echo "      -n <number of procesors> => Number of processors (for each job) (optional: default=1)"
	echo "      -xi <fragment length> => Minimum fragment length for paired-end reads (optional: default=0)"
	echo "      -mi <fragment length> => Maximum fragment length for paired-end reads (optional: default=500)"
	echo "      -t <allocated time> => Allocated time (in hours) for the analysis (optional: default=24)"
	echo "      -m <allocated memory> => Allocated memory (in gigabytes) per processor for each analysis (optional: default=4)"
	echo "      -q <queue> => PBS queue: bluecrystal (fixed for bluecrystal usage)"
	echo "      -e <email> => Notification email address (default=none)"
	echo "      -j <job name> => name of job for bluecrystal queue
	echo "      -h => show this help"
	echo ""
	echo "  Example:"
	echo "      $(basename $0) -i input_dir -o output_dir -r reference -j job_name"
	echo ""
	echo ""
	exit 0
}
# -----------------------------------------------------------------------------

author

# Get options from the command line
# -----------------------------------------------------------------------------
if [ "$#" -ge "6" ]; # min 6 args: 2 for -i <input directory>, 2 for -o <output directory>, 2 for -r <reference file>
then 
	while [ $# -gt 0 ]; do
		case "$1" in
			-h|-help) usage
					  ;;
			-i)	shift
				INDIR=$(readlink -f $1)
				;;
			-p)	shift
				INDIRP=$(readlink -f $1)
				;;
			-u)	shift
				INDIRU=$(readlink -f $1)
				;;
			-o)	shift
				OUTDIR=$(readlink -f $1)
				;;
			-r)	shift
				REFERENCE=$(readlink -f $1)
				;;
			-s)	shift
				ENC=$1
				;;			
			-n)	shift
				NTHR=$1
				;;
			-xi)shift
				MININS=$1
				;;			
			-mi)shift
				MAXINS=$1
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
			-j)	shift
				JOBNAME=$(readlink -f $1)
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
# -----------------------------------------------------------------------------

# Check $INDIR, $OUTDIR are $REFERENCE are defined
# -----------------------------------------------------------------------------
if [[ -z $INDIR || -z $OUTDIR || -z $REFERENCE ]]; then
	if [ -z $INDIR ]; then
		echo
		echo "ERROR: You must specify an input directory"
		echo
		ERR=1
	fi
	if [ -z $OUTDIR ]; then
		echo
		echo "ERROR: You must specify an output directory"
		echo
		ERR=1
	fi
	if [ -z $REFERENCE ]; then
		echo
		echo "ERROR: You must specify a reference file"
		echo
		ERR=1
	fi
	usage
fi
# -----------------------------------------------------------------------------

# Check for paired-end reads
# -----------------------------------------------------------------------------
if [[ ! -z $INDIRP ]]; then
	PAIREND=1
	if [[ ! -z $INDIRU ]]; then
		UNPAIRED=1
	fi
fi
# -----------------------------------------------------------------------------

# Check $INDIR and $REFERENCE exist (and $INDIRP if specified)
# -----------------------------------------------------------------------------
if [ ! -d $INDIR ]; then
	echo
	echo "ERROR: I can't find the input directory $INDIR"
	echo
	ERR=1
fi

if [[ $PAIREND -eq 1 && ! -d $INDIRP ]]; then
	echo
	echo "ERROR: I can't find the pair-end reads input directory $INDIRP"
	echo
	ERR=1
fi

if [[ $PAIREND -eq 1 && ! -z $INDIRU && ! -d $INDIRU ]]; then
	echo
	echo "ERROR: I can't find the unpaired reads input directory $INDIRU"
	echo
	ERR=1
fi

if [ ! -f $REFERENCE ]; then
	echo
	echo "ERROR: I can't find the reference file $REFERENCE"
	echo
	ERR=1
fi

# -----------------------------------------------------------------------------


# Get number of fastq files in the input directory
N=$(find $INDIR -maxdepth 1 | grep -P '.+\.f(ast)?q(\.(gz|bz2))?' | wc -l | cut -f1 -d" ")

# Check if there are fastq files in the input directory
# -----------------------------------------------------------------------------
if [ "$N" -eq "0" ]; then
	echo
	echo "ERROR: No fastq files in $INDIR"
	echo
	ERR=1
fi
# -----------------------------------------------------------------------------

# Check if there are fastq files in the paired input directory
# -----------------------------------------------------------------------------
if [ $PAIREND -eq 1 ]; then
	N2=$(find $INDIRP -maxdepth 1 |  grep -P '.+\.f(ast)?q(\.(gz|bz2))?' | wc -l | cut -f1 -d" ")

	if [ "$N2" -eq "0" ]; then
		echo
		echo "ERROR: No fastq files in $INDIRP"
		echo
		ERR=1
	else
		# Check if #files reads1 = # files reads2
		if [ $N -ne "$N2" ]; then
			echo
			echo "ERROR: Different number of fastq files for reads1 ($N) and reads2 ($N2)"
			echo
			ERR=1
		fi
	fi
	if [ ! -z $INDIRU ];
	then
		NU=$(find $INDIRU -maxdepth 1 |  grep -P '.+\.f(ast)?q(\.(gz|bz2))?' | wc -l | cut -f1 -d" ")
		if [ "$NU" -eq "0" ]; then
			echo
			echo "ERROR: No fastq files in $INDIRU"
			echo
			ERR=1
		fi
	fi
fi
# -----------------------------------------------------------------------------

# Quality score encoding 
# -----------------------------------------------------------------------------
if [ $ENC -eq 0 ];
then
	QSENC='--phred64'
elif [ $ENC -eq 1 ];
then
	QSENC='--phred33'
else
	echo
	echo "ERROR: quality score incorrect"
	echo
	ERR=1
fi
# -----------------------------------------------------------------------------

# Stop and show help if there are errors
# -----------------------------------------------------------------------------
if [ $ERR -eq 1 ]; then
	echo "Aborted because of errors - ERR $ERR"
	exit
fi
# -----------------------------------------------------------------------------

# Index reference if it is not indexed
# -----------------------------------------------------------------------------
if [ ! -e $REFERENCE.fai]; then
	echo
	echo "Building indexes for reference $REFERENCE in $REFERENCEDIR..."
	echo "(this might take quite some time...)"
	echo
	cd $REFERENCEDIR
	bwa index $REFERENCE >& $REFERENCE.bwamem.log
	echo
	echo "Done. Log file saved in $REFERENCE.bwamem.log"
	echo
	echo "----"
	echo
fi
# -----------------------------------------------------------------------------

# Create output directory
# -----------------------------------------------------------------------------
if [ ! -d $OUTDIR ]; then
	mkdir -p $OUTDIR
fi
# -----------------------------------------------------------------------------

# Memory for samtools sort set to 1/2 of requested memory (in MB)
# -----------------------------------------------------------------------------
SAMMEM=$(perl -e 'printf ("%.0f",(1000*'$NTHR'*'$MEM'/2));')
# SAMMEM=5 # debug
# -----------------------------------------------------------------------------

# Temporary file for sorting files with samtools
# -----------------------------------------------------------------------------
TIMESTAMP=$(date +%Y%m%d-%H%M%S)
GLOBTMPDIR="$USER-map_reads_samtools_sort_"$TIMESTAMP"_"$RANDOM

#if [[ $QUEUE == "popgenom" || $QUEUE == "molecol" ]];
#then
#	GLOBTMPDIR="/local/tmp/$GLOBTMPDIR"
#else
#	GLOBTMPDIR="/scratch/$GLOBTMPDIR"
#fi
# -----------------------------------------------------------------------------



# Set script and log filenames
TIMESTAMP=$(date +%Y%m%d-%H%M%S)
SMSJOB="$OUTDIR/mapping.$TIMESTAMP.smsjob.sh"
LOG="$OUTDIR/mapping.$TIMESTAMP.smsjob.log"
SMSJOBSUM="$OUTDIR/mapping.$TIMESTAMP.smsjobsum.sh"
LOGSUM="$OUTDIR/mapping.$TIMESTAMP.smsjobsum.log"

# Make submission file
# -----------------------------------------------------------------------------
# Initialize submission script
echo '#!/bin/bash' > $SMSJOB

# PBS OPTIONS
# -----------------------------------------
echo '#PBS -N '$JOBNAME'' >> $SMSJOB  
echo '#PBS -l nodes=1:ppn=1' >> $SMSJOB 
echo '#PBS -l mem='$MEM'gb' >> $SMSJOB 
echo '#PBS -l walltime='$HRS':00:00' >> $SMSJOB
echo '#PBS -j oe' >> $SMSJOB 
echo '#PBS -t 1-'$N >> $SMSJOB
echo '#PBS -o '$LOG >> $SMSJOB
# -----------------------------------------

cat >> $SMSJOB <<EOF


# -----------------------------------------
# Run job in local directory 
# -----------------------------------------
echo 'Run in local directory' >> $SMSJOB
echo 'cd $PBS_O_WORKDIR' >> $SMSJOB
# -----------------------------------------



# -----------------------------------------
# Load modules 
# -----------------------------------------

module load apps/bwa-0.7.15
module load apps/samtools-1.8
module load languages/perl-5.14.2

# -----------------------------------------


INDEX=\$((PBS_ARRAYID-1))

INDIR=$INDIR
FQFILES=(\$(find \$INDIR -maxdepth 1 | grep -P '.+\.f(ast)?q(\.(gz|bz2))?' | sort -V))
FQ=\${FQFILES[\$INDEX]}

INEXT=\${FQ##*.}
READIN='cat'
if [ "\$INEXT" == "bz2" ];then
	READIN='bzip2 -dck'
	# READIN='lbzip2 -dck -n $NTHR'
elif [ "\$INEXT" == "gz" ];then 
	READIN='gzip -dc'
	# READIN='gzip -dck -p $NTHR'
fi


EOF

if [ $PAIREND -eq 1 ];then

cat >> $SMSJOB << EOF
INDIRP=$INDIRP
FQFILES2=(\$(find \$INDIRP -maxdepth 1 | grep -P '.+\.f(ast)?q(\.(gz|bz2))?' | sort -V))
FQ2=\${FQFILES2[\$INDEX]}
	
INPUT="-I $MININS -X $MAXINS -1 <(\$READIN \$FQ) -2 <(\$READIN \$FQ2)"

EOF
else # single reads
	echo 'INPUT="-U <($READIN $FQ)"' >> $SMSJOB
fi

cat >> $SMSJOB <<EOF

FQBN=\$(basename \$FQ)
OUTPUT="$OUTDIR/\${FQBN%.[fq|fastq]*}"

TMPDIR="$GLOBTMPDIR/\${FQBN%.[fq|fastq]*}"
TMPSRTPFX="\$TMPDIR/\${FQBN%.[fq|fastq]*}"

# create temporary directory
mkdir -p \$TMPDIR

IDLOG="$OUTDIR/\${FQBN%.[fq|fastq]*}.mapping.log"

# echo \$INPUT

hostname > \$IDLOG
uname -a >> \$IDLOG
date >> \$IDLOG
echo "---------------------------------------------------" >> \$IDLOG
echo >> \$IDLOG
echo "Mapping \$FQ file to $REFERENCE..." >> \$IDLOG
echo >> \$IDLOG
echo "CMD: " >> \$IDLOG
echo "    $MAPPER $MAPPER_ARGS $QSENC -p $NTHR -x $REFERENCE -q \$INPUT | \\\" >> \$IDLOG
echo "    $SAMTOOLS view -b - | \\\" >> \$IDLOG
echo "    $SAMTOOLS sort - -m "$SAMMEM"M -T \$TMPSRTPFX -o \$OUTPUT.sorted.bam" >> \$IDLOG
echo  >> \$IDLOG
echo "    $SAMTOOLS index \$OUTPUT.sorted.bam" >> \$IDLOG
echo >> \$IDLOG
echo "---------------------------------------------------" >> \$IDLOG
echo >> \$IDLOG

bash -c "$MAPPER $MAPPER_ARGS $QSENC -p $NTHR -x $REFERENCE -q \$INPUT" 2>> \$IDLOG | \\
$SAMTOOLS view -b -  2>> \$IDLOG | \\
$SAMTOOLS sort - -m "$SAMMEM"M -T \$TMPSRTPFX -o \$OUTPUT.sorted.bam >> \$IDLOG 2>&1

$SAMTOOLS index \$OUTPUT.sorted.bam >> \$IDLOG 2>&1

# find unpaired containing sample with same id	
SID=\$(\$READIN \$FQ | head -n1 | perl -pe 's/ *\-.*//g')
FQU=\$(find $INDIRU -maxdepth 1 | grep -P '.+\.f(ast)?q(\.(gz|bz2))?' | \\
	   parallel -j $NTHR 'echo -en "{}\t"; '\$READIN' {} | head -n1 | perl -pe "s/ *\-.*//g";' | \\
	   grep -P "\$SID$" | cut -f1)

if  [[ ! -z \$FQU ]]; then
	# align unpaired reads

	echo >> \$IDLOG
	date >> \$IDLOG
	echo "---------------------------------------------------" >> \$IDLOG
	INPUT="-U <(\$READIN \$FQU)"
	echo >> \$IDLOG
	echo "Mapping \$FQU file to $REFERENCE..." >> \$IDLOG
	echo >> \$IDLOG
	echo "CMD: " >> \$IDLOG
	echo "    $MAPPER $MAPPER_ARGS $QSENC -p $NTHR -x $REFERENCE -q \$INPUT | \\\" >> \$IDLOG
	echo "    $SAMTOOLS view -b - | \\\" >> \$IDLOG
	echo "    $SAMTOOLS sort - -m "$SAMMEM"M -T \$TMPSRTPFX -o \$OUTPUT-u.sorted.bam" >> \$IDLOG
	echo  >> \$IDLOG
	echo "    $SAMTOOLS index \$OUTPUT-u.sorted.bam" >> \$IDLOG
	echo >> \$IDLOG
	echo "---------------------------------------------------" >> \$IDLOG

	bash -c "$MAPPER $MAPPER_ARGS $QSENC -p $NTHR -x $REFERENCE -q \$INPUT" 2>> \$IDLOG | \\
	$SAMTOOLS view -b -  2>> \$IDLOG | \\
	$SAMTOOLS sort - -m "$SAMMEM"M -T \$TMPSRTPFX -o \$OUTPUT-u.sorted.bam >> \$IDLOG 2>&1

	$SAMTOOLS index \$OUTPUT-u.sorted.bam >> \$IDLOG 2>&1


	mv \$OUTPUT.sorted.bam \$OUTPUT-p.sorted.bam
	mv \$OUTPUT.sorted.bam.bai \$OUTPUT-p.sorted.bam.bai

	# merge bam files
	echo >> \$IDLOG
	date >> \$IDLOG
	echo "---------------------------------------------------" >> \$IDLOG
	echo >> \$IDLOG
	echo "Merging \$OUTPUT-p.sorted.bam \$OUTPUT-u.sorted.bam..." >> \$IDLOG
	echo >> \$IDLOG
	echo "CMD: " >> \$IDLOG
	echo "    $SAMTOOLS merge \$OUTPUT.sorted.bam \$OUTPUT-p.sorted.bam \$OUTPUT-u.sorted.bam" >> \$IDLOG
	echo  >> \$IDLOG
	echo "    $SAMTOOLS index \$OUTPUT.sorted.bam" >> \$IDLOG
	echo >> \$IDLOG
	echo "---------------------------------------------------" >> \$IDLOG


	$SAMTOOLS merge \$OUTPUT.sorted.bam \$OUTPUT-p.sorted.bam \$OUTPUT-u.sorted.bam >> \$IDLOG 2>&1
	$SAMTOOLS index \$OUTPUT.sorted.bam >> \$IDLOG 2>&1

	rm -f \$OUTPUT-p.sorted.bam* \$OUTPUT-u.sorted.bam*
fi

# delete temporary directory
rm -rf \$TMPDIR

echo >> \$IDLOG
echo "---------------------------------------------------" >> \$IDLOG
echo >> \$IDLOG
date >> \$IDLOG
echo >> \$IDLOG

EOF

chmod +x $SMSJOB
# -----------------------------------------------------------------------------

# Add summary script
# Make summary script (simply gather report of mapped reads)
# -----------------------------------------------------------------------------

# ls *.bam | \
# parallel -j 20 'S=$(echo {} | perl -pe 's/\.sorted.*//g'); P=$(samtools view -F12 {} | wc -l); U=$(samtools view -F8 -f4 {} | wc -l); N=$(samtools view -f4 {} | wc -l); T=$(perl -e "print $P+$U"); echo -e "$S\t$P\t$U\t$T\t$N"' | sort -V > body
# -----------------------------------------------------------------------------

echo "Command to submit the job to Iceberg ($QUEUE queue):"
echo
echo "qsub $SMSJOB"
echo
# echo "qsub $SMSJOBSUM"
# echo
