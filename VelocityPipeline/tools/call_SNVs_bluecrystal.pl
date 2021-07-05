#!/usr/bin/env perl

# (c) Victor Soria-Carrasco
# victor.soria.carrasco@gmail.com
# Last modified: 23/10/2018 09:23:41

# Description:
# This script will call SNPs using samtools mpileup/bcftools, and will filter them
# with vcfutils. See 'Hardcoded constants' section if you want to change bcftools options.
# See parameters file example below.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


# Changelog
# 1.2 - 09/03/2016 
#	  - Updated to use bowtie 2.2.7
#	  - Updated to use samtools/bcftools 1.3
#	  - Added option to call variants for a set of regions only
# 1.2.1 - 16/03/2016
#     - Added control for empty bcf files
#     - Added mapping quality score filter
#     - Added option to choose between consensus caller and multi-allelic caller
# 1.2.2 - 25/03/2016
#     - Cosmetic changes
# 1.2.3 - 05/07/2016
#     - Updated samtools and bcftools to version 1.3.1
# 1.2.4 - 11/05/2017
#     - Increased number of open files to 10240 (ulimit -n 10240)
#
# 1.2.5 - 03/08/2017
#     - More accurate use of the max number of open files with ulimit (=No regions + 100)
#
# 1.2.5 - 04/08/2017
#     - Updated samtools and bcftools to version 1.5
#
# 1.2.6 - 11/08/2017
#      - Set the hard limit for the number of open files 10000 
#        (maximum allowed on Iceberg and ShARC HPC clusters)
#
# 1.2.7 - 23/11/2017
#      - Added option to call also invariants
#      - It now also call indels (removed flag from samtools mpileup preventing it)
#
# 1.2.8 - 25/01/2018
#	   - Fixed bug related to using ulimit from perl forks (call of the shell required, i.e. bash -c ulimit)
#	   - Fixed subroutine for progress bar (stopped working with ForkManager)
#
# 1.2.9 - 23/04/2018
#      - Updated samtools and bcftools to version 1.8
#
#     ToDo: add option to call variants from all regions at once?
#
# 1.3 - 23/10/2018 AJvR
#	- Modify to run on BlueCrystal
#	- Callum installed perl packages Parallel and Term
#	- set path for samtools and 
#
#



use strict;
use warnings;
use File::Path qw(make_path rmtree);
use File::Basename;
use File::Spec;
use Getopt::Long;
use POSIX qw/ strftime /;
use Parallel::ForkManager;
use Term::ProgressBar;
use subs 'printR'; # Special print function, print to both STDOUT and report file

my $version='1.3-2018.10.23';

# Harcoded constants
# -------------------------------------------------------------------------------------
my $samtools="samtools"; # path to samtools
my $bcftools="bcftools"; # path to bcftools
my $vcfutils="/newhome/aj18951/bin/vcfutils.pl"; # path to vcfutils
# -------------------------------------------------------------------------------------

&author;
# Get options from the command line
# -------------------------------------------------------------------------------------
my ($inbamdir, $refseq, $outdir);
my $ncpus=1;
my @regions=();
my $pvar=0.05; # bcftools option - include SNP only if P(ref|D)<$pvar
			   # (probability of data under null model of all reads homozygous for ref allele)
# my $sampfreq=0.4; # bcftools option - minimum fraction of reads required to include a SNP
# my $minQS=20; # minimum quality score (phred)
# my $maxdepth=10000; # max SNP depth
my $caller='c'; # caller used in bcftools (c=consensus, m=multi-allelic)
my $minMQS=20; # minimum mapping quality score (phred)
my $together=0; # Call variants from all regions together
my $quiet=0; # Show/hide progress bar
my $cleantmp=0;
my $maxopenfiles=1024;
my $varonly=0;
my $indels=0;

GetOptions( 
    'b=s'       => \$inbamdir, 
    'r=s'       => \$refseq, 
    'o=s'       => \$outdir, 
    'regs=s{,}' => \@regions, 
    'n=i'       => \$ncpus, 
    # 's=f'     => \$sampfreq, 
    'c=s'       => \$caller, 
    'v=i'       => \$varonly, 
    'd=i'       => \$indels, 
    'q=i'       => \$minMQS, 
    'p=f'       => \$pvar, 
    'together'  => \$together, 
    'quiet'     => \$quiet, 
    'clean'     => \$cleantmp, 
    'help'      => \&usage 
);

&usage if (!defined($inbamdir) || !defined($refseq) || !defined($outdir));

$caller=lc($caller);
die ("\nERROR: caller must be 'c' or 'm' (option specified was -c $caller)\n\n")
	if ($caller ne "c" && $caller ne "m");

# Get absolute paths
# -------------------------------------------------------------------------------------
$inbamdir=File::Spec->rel2abs($inbamdir);
$refseq=File::Spec->rel2abs($refseq);
$outdir=File::Spec->rel2abs($outdir);
# -------------------------------------------------------------------------------------


# Output directory
# -------------------------------------------------------------------------------------
if (! -e $outdir){
	eval {make_path($outdir)}
		or die ("\nCan't create output directory: $outdir\n\n");
}
# -------------------------------------------------------------------------------------

# Temp directory
# -------------------------------------------------------------------------------------
my $tmpdir=$outdir."/tmp";
if (! -e $tmpdir){
	eval {make_path($tmpdir)}
		or die ("\nCan't create temporary directory: $tmpdir\n\n");
}
# -------------------------------------------------------------------------------------

# Regions
# -------------------------------------------------------------------------------------
if (defined($regions[0])){
	@regions=split(/,/,join(',',@regions));
}
# -------------------------------------------------------------------------------------

# Open report file
# -------------------------------------------------------------------------------------
my $repfile=$outdir."/".strftime("%Y%m%d-%H%M%S", localtime).".report.txt";
open (REPFILE, ">$repfile")
	or die ("\nCan't write report to file $repfile");

# Print execution information
# -------------------------------------------------------------------------------------
printR "-------------------------------------------------------------------------------\n\n";
printR basename($0)." executed on:\n\n";
printR "\t".strftime("%A %d %B %Y %Z %H:%M:%S", localtime)."\n\n";
printR "using the following command:\n\n";
printR "\t".File::Spec->rel2abs($0)."\n";
printR "\t\t-b $inbamdir\n";
printR "\t\t-r $refseq\n";
printR "\t\t-o $outdir\n";
printR "\t\t-regs @regions\n";
printR "\t\t-n $ncpus\n";
printR "\t\t-c $caller\n";
printR "\t\t-v $varonly\n";
printR "\t\t-q $minMQS\n";
printR "\t\t-p $pvar\n";
printR "\t\t-together\n" if ($together ==1);
printR "\t\t-quiet\n" if ($quiet ==1);
printR "\t\t-clean\n" if ($cleantmp ==1);
printR "\n";
printR "-------------------------------------------------------------------------------\n\n";


# Index reference if needed
# -------------------------------------------------------------------------------------
if (!-e "$refseq.fai"){
	printR "\nThe reference sequence is not indexed. Indexing reference file...\n\n";
	system ("$samtools faidx $refseq"); # Index fasta file
	die ("\nsamtools failed to index reference fasta file\n\n")	
		if (!-e "$refseq.fai");
	printR "Indexing finished\n\n";
}
# -------------------------------------------------------------------------------------

 
# Call variants using samtools mpileup + bcftools
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# call SNPs with samtools + bcftools for each scaffold/chromosome
# -------------------------------------------------------------------------------------
my $vartxt="variants";
$vartxt.=" and invariants" if ($varonly==0);

my $lgord=0;
if (!defined($regions[0])){
	printR "Calling $vartxt for each scaffold/chromosome using samtools/bcftools...\n\n";
	# Check if we have Timema linkage groups
	# and sort scaffolds accordingly
	$lgord=`grep -chP "lg([0-9]{1,}|NA*)_ord([0-9]{1,}|NA*)_scaf([0-9]{1,})" $refseq.fai`;
	if ($lgord > 0){
		@regions=sort {
			my ($alg) = $a =~ /lg([0-9]*)/;
			my ($blg) = $b =~ /lg([0-9]*)/;
			my ($aord) = $a =~ /ord([0-9]*)/;
			my ($bord) = $b =~ /ord([0-9]*)/;
			my ($ascaf) = $a =~ /scaf([0-9]*)/;
			my ($bscaf) = $b =~ /scaf([0-9]*)/;
			$alg <=> $blg || $aord <=> $bord || $ascaf <=> $bscaf;
		} `cut -f1 $refseq.fai | sed 's/NA/999999999/g'`;
	}
	else{
		@regions=`cut -f1 $refseq.fai`;
	}
}
else{
	printR "Calling $vartxt for regions: ".join(',',@regions)." using samtools/bcftools...\n\n";
}

# Arguments for samtools mpileup and bcftools call
# my $mpileupargs="-q $minMQS -g -t DP,SP,AD,INFO/AD -P ILLUMINA -u -f $refseq"; # platform option prevents calling indels
my $mpileupargs="-q $minMQS -g -t DP,SP,AD,INFO/AD -u -f $refseq";
# my $callargs="-p $pvar -O b -cv";
my $callargs="-p $pvar -O b -c";
if ($caller eq 'm'){
	# $callargs="-p $pvar -O b -mv -f GP,GQ";
	$callargs="-p $pvar -O b -m -f GP,GQ";
}
$callargs.=" -v" if ($varonly==1);
$callargs.=" -V indels" if ($indels==0);

my ($progress, $next_update); #required for progress bar

my $nregs=scalar(@regions);
$maxopenfiles=$nregs+100 if ($nregs > $maxopenfiles);
$maxopenfiles=10000 if ($maxopenfiles > 10000); # max hard limit in Iceberg and ShARC

if ($together==1){
	my $allreg=join(',',@regions);
	if ($varonly == 1){
		system("bash -c 'ulimit -n $maxopenfiles; cd $inbamdir && $samtools mpileup $mpileupargs -r $allreg *.bam 2> $tmpdir/var_all.raw.log | \
		$bcftools call $callargs - 2>> $tmpdir/var_all.raw.log > $tmpdir/var_all.raw.bcf'");
	}
	else{
		system("bash -c 'ulimit -n $maxopenfiles'; cd $inbamdir && $samtools mpileup $mpileupargs -r $allreg *.bam 2> $tmpdir/var_all.raw.log | \
		$bcftools call $callargs - 2>> $tmpdir/var_all.raw.log | \
		$bcftools view -O b -U - 2>> $tmpdir/var_all.raw.log > $tmpdir/var_all.raw.bcf");
	}
	bcfchkidx("$tmpdir/var_all.raw.bcf");
}
else{
	# Using ForkManager
	# my $nregs=scalar(@regions);

	# Progress bar
	# ..................................
	 my $progress=Term::ProgressBar->new ({
		name=>"\tProgress", 
		count=>$nregs, 
		ETA=>'linear', 
		term_width=>80,
		silent=>$quiet});

	 $progress->minor(0);
	 $progress->max_update_rate(1);
	 my $next_update=0;
	 # ..................................

	my $c=0;
	my $pm = Parallel::ForkManager->new(($ncpus-1));
	foreach my $reg (@regions){
		chomp($reg);
		$reg=~ s/999999999/NA/g if ($lgord > 0); # Put back NAs in Timema
		$c++;	
		# &update_bar($c) if (!$quiet && !($c % 10));
		my $pid = $pm->start and next;

		# print "$samtools mpileup $mpileupargs -r $reg\n"; # debug
		# print "$bcftools call $callargs\n"; # debug
		if ($varonly==1){
			system("bash -c 'ulimit -n $maxopenfiles; cd $inbamdir && $samtools mpileup $mpileupargs -r $reg *.bam 2> $tmpdir/var_$reg.raw.log | \
			$bcftools call $callargs - 2>> $tmpdir/var_$reg.raw.log > $tmpdir/var_$reg.raw.bcf'");
		}
		else{
			system("bash -c 'ulimit -n $maxopenfiles; cd $inbamdir && $samtools mpileup $mpileupargs -r $reg *.bam 2> $tmpdir/var_$reg.raw.log | \
			$bcftools call $callargs - 2>> $tmpdir/var_$reg.raw.log | \
			$bcftools view -O b -U - 2>> $tmpdir/var_$reg.raw.log> $tmpdir/var_$reg.raw.bcf'");
		}
		bcfchkidx("$tmpdir/var_$reg.raw.bcf");
		&update_bar($c,$next_update,$progress) if (!$quiet);
		$pm->finish();
	}
	$pm->wait_all_children;

	$progress->update($nregs) if ($nregs>$next_update && !$quiet);
}

printR "\nFinished.\n";
# -------------------------------------------------------------------------------------

# Merge results
# -------------------------------------------------------------------------------------
opendir(TMPDIR, $tmpdir)
	or die ("\nCan't open temporary directory $tmpdir\n\n");
	my @bcfs=grep(/var\_.*\.raw\.bcf$/, readdir(TMPDIR));
closedir (TMPDIR);

if (scalar(@bcfs) > 0){
	printR "\nMerging results for all scaffolds/chromosomes...\n";
	# Sort by linkage group, order, and scaffold (for Timema)
	my @list=();
	if ($lgord > 0){
		@list=sort {
			my ($alg) = $a =~ /lg([0-9]*)/;
			my ($blg) = $b =~ /lg([0-9]*)/;
			my ($aord) = $a =~ /ord([0-9]*)/;
			my ($bord) = $b =~ /ord([0-9]*)/;
			my ($ascaf) = $a =~ /scaf([0-9]*)/;
			my ($bscaf) = $b =~ /scaf([0-9]*)/;
			$alg <=> $blg || $aord <=> $bord || $ascaf <=> $bscaf;
		} `ls $tmpdir/var_*.raw.bcf | sed 's/NA/999999/g'`;
	}
	else{
		@list=sort {
			my ($an) = $a =~ /.*([0-9]+)\:?/;
			my ($bn) = $b =~ /.*([0-9]+)\:?/;
			($an) <=> ($bn);
		} `ls $tmpdir/var_*.raw.bcf`
	}
	# Put back NAs and print ordered list to temporary file
	open (TMP, ">$tmpdir/scafforder.tmp")
		or die ("\nCan't write to temporary file $tmpdir/scafforder.tmp\n\n");
		foreach (@list){
			s/999999/NA/g if ($lgord > 0);
			print TMP $_;
		}
	close (TMP);

	# Merge all files
	system("$bcftools concat -O b -a -d none -f $tmpdir/scafforder.tmp 2> $tmpdir/variants.raw.merging.log >$outdir/variants.tmp.raw.bcf");
	system("$bcftools index $outdir/variants.tmp.raw.bcf");

	# Removed path or extension from sample names
	my $smpnames=`$bcftools view $outdir/variants.tmp.raw.bcf | grep -m1 'CHROM'`;
	$smpnames=~ s/.*FORMAT\t//g;
	$smpnames=~ s/(\.sort(ed)?)?\.bam//g;
	$smpnames=~ s/\t/\n/g;
	open (FILE, ">$outdir/newsmpnames")
		or die ("\nCannot write to $outdir/newsmpnames\n\n");
		print FILE $smpnames;
	close (FILE);

	my $outbcf="variants.raw.bcf";
	$outbcf="sites.raw.bcf" if ($varonly==0);

	system("$bcftools reheader -s $outdir/newsmpnames $outdir/variants.tmp.raw.bcf > $outdir/$outbcf");
	system("$bcftools index $outdir/$outbcf");
	unlink glob ("$outdir/variants.tmp.raw.bcf");
	unlink glob ("$outdir/variants.tmp.raw.bcf.csi");
	unlink glob ("$outdir/newsmpnames");

	# Get number of samples and variants
	my $nsmp=()=$smpnames=~ /\n/g;
	my $nvar=`$bcftools view -H -v snps,indels $outdir/$outbcf | wc -l`;
	chomp($nvar);
	
	my $vartxt="$nvar variants";

	if ($varonly==0){
		my $ninvar=`$bcftools view -H $outdir/$outbcf | wc -l`;
		$ninvar=$ninvar-$nvar;
		$vartxt.=" and $ninvar invariants";
	}

	printR "$vartxt called for $nsmp samples/individuals.\n\n";
	printR "Results merged into $outdir/$outbcf\n\n";
}
else{
	my $vartxt="variants";
	$vartxt.=" nor invariants" if ($varonly==0);
	printR "Warning: No $vartxt found for any region\n\n";
}
printR "Finished.\n";
# -------------------------------------------------------------------------------------

# Clean temporary files
# -------------------------------------------------------------------------------------
rmtree($tmpdir) if ($cleantmp==1);
# -------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

printR "All processes finished.\n\n";

close (REPFILE);

# ==============================================================================
# ==============================================================================
# ============================== SUBROUTINES ===================================
# ==============================================================================
# ==============================================================================

# Show copyright
# ==============================================================================
sub author{
    print "\n";
    print "#########################################\n";
    print "  ".basename($0)."\n";
	print "  version $version     \n";
    print "  (c) Victor Soria-Carrasco             \n";
    print "  victor.soria.carrasco\@gmail.com      \n";
    print "#########################################\n";
	print "\n";
}
# ==============================================================================

# Show usage
# ==============================================================================
sub usage{
    print "\n";
	print "  Usage:\n";
    print "    ".basename($0)."\n";
	print "      -b <input directory> => directory with bam files (1 file per sample)\n";
	print "      -r <fasta file> => reference sequence used for alignment\n";
	print "      -o <output directory>\n";
	print "      -regs <regions> => restrict to regions (i.e. scaffolds or chromosomes), separated by spaces or commas\n";
	print "      -n <number of cores> => optional, default=$ncpus\n";
	# print "      -s <minimum sampling to retain a SNP (optional, default=0.4)>\n";
	print "      -c <c|m> => caller to be used with bcftools (c=consensus, m=multi-allelic) (optional, default=$caller)\n";
	print "      -v <0|1> => call only variants (optional, default=$varonly)\n";
	print "      -d <0|1> => call also indels (optional, default=$indels)\n";
	print "      -q <phred quality score> => minimum phred mapping quality score to consider an alignment for calling (optional, default=$minMQS)\n";
	print "      -p <float 0-1> => prob. of data under the hypothesis that locus is invariant (optional, default=$pvar)\n";
	print "      -together => consider all regions together to call variants (default: call variants separately, optional)\n";
	print "      -quiet => quiet output (don't show progress bar, optional)\n";
	print "      -clean => remove temporary files (kept by default, optional)\n";
    print "\n";
    exit;
}
# ==============================================================================

# Print to report file and STDOUT
# ==============================================================================
sub printR{
	foreach (@_){
		print $_;
		print REPFILE $_;  
	}
}
# ==============================================================================

# Update progress bar
# ==============================================================================
sub update_bar{
    my $c=shift(@_);
    my $next_update=shift(@_);
	my $progress=shift(@_);
    $next_update=$progress->update($c) if ($c>=$next_update);
}
# ==============================================================================

# Check empty bcf, index if exists
# ==============================================================================
sub bcfchkidx{
	my $bcf=shift;
	my $checkbcf=`file $bcf | grep -c empty`;
	if ($checkbcf==1){
		unlink glob ("$bcf");
		printR "Warning: empty bcf file was deleted\n";
	}
	else{
		system ("$bcftools index $bcf");
	}
}
# ==============================================================================
