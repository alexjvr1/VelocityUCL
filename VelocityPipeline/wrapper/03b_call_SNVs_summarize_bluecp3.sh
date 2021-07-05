#!/bin/bash

# (c) Victor Soria-Carrasco
# victor.soria.carrasco@gmail.com
# Last modified: 26/04/2018 13:35:57
#########################################################################
# Modified by Alexandra Jansen van Rensburg for Blue Crystal (UoBristol)
# alexjvr@gmail.com
# Last modified: 23/10/2018 14:12:00
#########################################################################

# Description:
# Given an input directory with bcf files with variants called
# for a number of regions, it merges and filter them

# Changelog
#	1.01 - 02/04/2016
#	- Fixed bug when plotting stats
#	1.02 - 26/04/2016
#	- Fixed bug output name didn't include all filters
#	- Set genotypes of samples without reads to missing (.)
#	1.03 - 05/07/2016
#	- Updated bcftools to version 1.3.1
#	1.04 - 15/05/2017
#	- Added option to avoid producing any filtered dataset
#	- Added option to use multithreading bcftools capabilities
#   1.05 - 09/08/2017
#   - Updated bcftools to version 1.5
#   1.06 - 12/11/2017
#   - Added option to process output including invariants
#   1.07 - 23/04/2018
#   - Updated bcftools to version 1.8
#   1.08 - 26/04/2018
#   - Removed output of statistics at the end
#   
#   2.00 - 23/10/2018
#   - Modified by AJvR for bluecrystal
#   - Add conversion from bcf to vcf	
#
#
#
#
#

#Define module versions to be loaded

module load apps/bcftools-1.8
module load languages/perl-5.14.2
module load languages/java-jdk-1.8.0-66
#
#
#
#  
# hardcode variables for tools
BCFTOOLS="bcftools"
PLOTSTATS="plot-vcfstats"
VCFUTILS="/newhome/aj18951/bin/vcfutils.pl"

VERSION='2.00-2018.10.23'

# Default values for optional variables
NTHREADS=1
FLT=0
SAMPCOV="0.4"
QS=20
MINDP=3
MAXDP=100000
NOINDEL=0
NOMULTI=0
NOPRIV=0
CLEANTMP=0
MINDISGAP=0

function author {
	echo
	echo "#########################################"
	echo "  $(basename $0)"
	echo "  version $VERSION"
	echo "  (c) Victor Soria-Carrasco"
	echo "  victor.soria.carrasco@gmail.com"
	echo "#########################################"
	echo
}

function usage {
	echo
	echo "Usage:"
	echo "  $(basename $0)"
	echo "      -i  <input directory> => With files from call_SNVs_bluecrystal.sh"
	echo "      -n  <number of threads> (optional, default=1)"
	echo "      -flt => Apply filters after merging bcfs (optional, default=doesn't do any filtering)"
	echo "      -scov <sampling coverage> => Minimum sampling coverage to retain a variant (optional, default=$SAMPCOV)"
	echo "      -qs <phred quality score> => Minimum quality score to retain a variant (optional, default=$QS)"
	echo "      -mindp <depth> => Minimum depth to retain a variant (optional, default=$MINDP)"
	echo "      -maxdp <depth> => Maximum depth to retain a variant (optional, default=$MAXDP)"
	echo "      -mindisgap <no base pairs> => Minimum distance to gap (indel) to retain a variant (optional, default=$MINDISGAP)"
	echo "      -noindel => Exclude indels"
	echo "      -nomulti => Exclude multiallelic variants (>2 alleles)"
	echo "      -nopriv => Exclude private variants (fixed, but different from reference)"
	echo "      -clean => Remove temporary files"
	echo "      -h => show this help"
	echo ""
	echo "  Example:"
	echo "      $(basename $0) -i ids_reads"
	echo ""
	echo ""
	exit 0
}

author

if [ "$#" -ge "2" ]; # min 2 args: 2 for -i <input directory>
then 
	while [ $# -gt 0 ]; do
		case "$1" in
			-h|-help) 
				usage
				;;
			-i)	shift
				INDIR=$(readlink -f $1)
				;;
			-n)	shift
				NTHREADS=$1
				;;
			-flt)
				FLT=1
				;;
			-scov)	shift
				SAMPCOV=$1
				;;
			-qs) shift
				QS=$1
				;;
			-mindp) shift
				MINDP=$1
				;;
			-maxdp) shift
				MAXDP=$1
				;;
			-mindisgap) shift
				MINDISGAP=$1
				;;
			-noindel)
				NOINDEL=1
				;;
			-nomulti)
				NOMULTI=1
				;;
			-nopriv)
				NOPRIV=1
				;;
                        -module1)       shift
                                        MODULE1=$1
                                        ;;
                        -module2)       shift
                                        MODULE2=$1
                                        ;;
                        -module3)       shift
                                        MODULE3=$1
                                        ;;
                        -module4)       shift
                                        MODULE4=$1
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

# Check existence input directory

if [ ! -d $INDIR ];
then
	echo
	echo "ERROR: Input directory $INDIR does not exist or is not accessible"
	echo
	exit
fi

TMPDIR="$INDIR/tmp"

mkdir -p $TMPDIR
mv $INDIR/*.bcf* $TMPDIR/
mv $INDIR/*.log $TMPDIR/
mv $INDIR/*.report.txt $TMPDIR/

# Get list of bcf files that are not empty
ls $TMPDIR/*.bcf | \
xargs -I {} sh -c '
	if [[ $(file {} | grep -c empty) == 0 ]];
	then 
		echo {}; 
	fi' | \
sort -V > $TMPDIR/bcflist

OUTF=$(ls $TMPDIR/*.bcf | head -n 1 | perl -pe 's/.*\///g;s/job[0-9]+\.//g;s/\.bcf//g')

$BCFTOOLS concat -O b -a -d none -f $TMPDIR/bcflist --threads $NTHREADS 2> $INDIR/$OUTF.log | \
#$BCFTOOLS concat -O u -a -D -f $TMPDIR/bcflist --threads $NTHREADS 2> $INDIR/$OUTF.log | \
#$BCFTOOLS filter -O b -S . -e 'FMT/DP=0' --threads $NTHREADS 2>> $INDIR/$OUTF.log > $INDIR/$OUTF.bcf
$BCFTOOLS index $INDIR/$OUTF.bcf

# plot stats
# $BCFTOOLS stats $INDIR/$OUTBCF 2> $INDIR/$OUTF.stats.log > $INDIR/$OUTF.stats
# $PLOTSTATS -p $INDIR/$OUTF.stats.plot/ -T $OUTF $INDIR/$OUTF.stats >> $INDIR/$OUTF.stats.log 2>&1

NSAMPLES=$($BCFTOOLS view -h $INDIR/$OUTF.bcf | grep CHROM | perl -pe 's/.*FORMAT\s//g' | wc -w)


if [[ $FLT == 1 ]];
then
	# Filters
	# ---------------------------------------------------------------------------------------------------------------
	echo
	echo "Filtering:"
	echo

	SCOV=$(echo $SAMPCOV | perl -pe '$_*=100')
	FLTNAME=$SCOV"cov"$QS"QS"$MINDP"minDP"$MAXDP"maxDP"$MINDISGAP"minDisGap"
	FLT="'(AN/2)>=N_SAMPLES*$SAMPCOV & QUAL>=$QS & MIN(DP)>=$MINDP & MAX(DP)<=$MAXDP"

	echo -e "\tsampling < $SCOV% excluded"
	echo -e "\tdepth < $MINDP excluded"
	echo -e "\tdepth > $MAXDP excluded"
	echo -e "\tquality score < $QS excluded"
	echo -e "\tminimum distance to gap (indel) < $MINDISGAP excluded"

	if [[ $NOINDEL == 0 ]];
	then
		FLT="$FLT & (TYPE=\"snp\" | TYPE=\"indel\")"
	else
		echo -e "\tindels excluded"
		FLT="$FLT & TYPE=\"snp\""
		FLTNAME=$FLTNAME"noIndel"
	fi

	if [[ $NOMULTI == 1 ]];
	then
		echo -e "\tmultiallelic (>2alleles) excluded"
		FLT="$FLT & N_ALT=1"
		FLTNAME=$FLTNAME"noMulti"
	fi

	if [[ $NOPRIV == 1 ]];
	then
		echo -e "\tprivate (ALT AF=1) excluded"
		FLT="$FLT & AC/AN<1"
		FLTNAME=$FLTNAME"noPriv"
	fi
	FLT="$FLT'"

	OUTFF=$(echo $OUTF | perl -pe 's/\.raw/\.flt/g')
	eval $BCFTOOLS filter -O b -i $FLT -g $MINDISGAP $INDIR/$OUTF.bcf --threads $NTHREADS 2> $INDIR/$OUTFF.$FLTNAME.log > $INDIR/$OUTFF.$FLTNAME.bcf
	$BCFTOOLS index $INDIR/$OUTFF.$FLTNAME.bcf

	NVARRAW=$($BCFTOOLS view $INDIR/$OUTF.bcf | grep -chv "^#")
	NVARFLT=$($BCFTOOLS view $INDIR/$OUTFF.$FLTNAME.bcf | grep -chv "^#")

	echo
	echo -e "$NVARFLT out of $NVARRAW variable sites retained after filtering.\n\n";
	echo
	# ---------------------------------------------------------------------------------------------------------------
else
	echo
	echo -e "$NVARRAW variable sites in $INDIR/$OUTF.bcf\n\n";
	echo
fi

if [ $CLEANTMP == 1 ];
then
	echo
	echo "Cleaning up..."
	echo
	rm -rf $TMPDIR
fi

echo "converting bcf to vcf"
$BCFTOOLS view $INDIR/$OUTF.bcf -Ov > $INDIR/$OUTF.vcf
$BCFTOOLS view $INDIR/$OUTFF.$FLTNAME.bcf -Ov > $INDIR/$OUTFF.$FLTNAME.vcf

echo "Finished"
echo
