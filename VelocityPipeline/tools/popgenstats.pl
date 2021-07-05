#!/usr/bin/env perl

# (c) Victor Soria-Carrasco
# Last modified: 27/04/2018 17:31:32

# Description:
# This script estimate fsts and a number of other population
# genetics statistics for two populations defined
# by a set of individuals using a bcf file as input (which can
# contain other populations as well)

# Mean genotype probabilities can be estimated using a uniform or
# a Hardy-Weinberg equilibrium prior
# Be aware that SNPs with >2 alleles will be excluded

# population pair file must have the following structure:
#
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop2

# windows file must have the following structure
# scaffold lower_bound upper_bound
#

# Possible improvements: 
#	MOST IMPORTANT  - Use only individuals with reads for sample size when estimating
#   Fst with Hudsonloci or Weir&Cockerham
# - Consider only samples with reads for each SNP
#	for AF and also for Fsts??
# - Use genotype probabilities as input file (bimbam)
# - Filter by sampling in each population
# - Filter by number of reads
# - sort by position in scaffold (not really necessary, because it 
#   is already done by bcftools)

# ToDo:
#  - simplify: remove all options and just calculate all statistics?
#  - add support for bcftools 1.x
#  - include Fu's Fs, ... 

# Changelog:
#  1.22 - 31/05/2016 - fixed small bug - uninitialized warning when (wrongly) 
#				       trying to sort windows when no windows specified
#  1.3  - 23/06/2016 - major changes
#					 - changed name: popgenstats.pl
#					 - added calculation of pi, dxy, da, Tajima's D
#					 - added support for bcftools 1.x (including automatic detection)
#  1.4  - 23/08/2017 - Updated to use bcftools 1.5, added option to filter by 
#                      population sampling coverage
#  1.4.1 - 05/09/2017 - Fixed several bugs: incorrect conversion from bcf to gl format,
#                       handling of regions when with no overlapping snvs (problematic
#                       for calculating Tajima D)
#  1.5 - 06/09/2017 - Changed bcftools merge command to not use threads - there is a big
#                    memory leak when combined with a list of regions
#                    Recoded to avoid loading array of genotype likelihoods and reduce 
#                    memory consumption
# 1.5.1 - 06/09/2017 - No intermediary bcf files written to disk required anymore,
#                      fixed bug when reading unsorted windows/regions file
# 1.5.2 - 27/04/2018 - Fixed help output (missing -msp description)
#

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path remove_tree);

# use Memory::Usage; # debug
# my $mu = Memory::Usage->new(); # debug
# $mu->record('starting work'); # debug

my $version='1.5.2-2018.04.27';

# Path to external programs
my $bcftools='bcftools';
my $estpEM='/newhome/aj18951/software/estpEM_2014-10-08/estpEM'; # program for estimating AFs by ML


&author;

# Read arguments and set up output directories/files
# =============================================================================
my $prior=0; # 0-> uniform, 1-> H-W
# my $samples=0;
my $mincov="0.5";
my $mincovp="0.5";
my $mmaf=0.05; # Minimum MAF allowed
my $separate=0;
my $afmethod=0; #0-> AFs estimated from empirical Bayesian genotypes, 1-> AFS estimated using EM algorithm
my $fstmethod=0; #0-> Hudson, 1->HudsonLoci, 2->Weir&Cockerham
my %fstmethod_names=(
	0=>"HudsonFst",
	1=>"HudsonFstloci",
	2=>"Weir&CockerhamFst");
my $dopi=0; # calculate pi (nucleotide diversity)
my $dodxy=0; # calculate dxy
my $doda=0; # calculate da
my $dotajimad=0; # calculate Tajima's d
my $seqlen=0; # required for calculations of pi per site, etc.
my $keeptmp=0;
my $ignorelg=0; # ignore lg and order for windows
my $nthreads=1;

&usage if (scalar(@ARGV)==0);

my ($infile, $outfile, $poppair, $windowsfile, $outfile2);
GetOptions(
    'i|I=s'     => \$infile,
    'p|P=s'     => \$poppair,  
    'w|W=s'     => \$windowsfile,  
    'o|O=s'     => \$outfile,
    'ow|OW=s'   => \$outfile2,
    'n|N=s'     => \$nthreads,
    'pr|PR=i'   => \$prior,
    'ms|MS=s'   => \$mincov,
    'msp|MSP=s' => \$mincovp,
    'mf|MF=s'   => \$mmaf,
    'am|AM=i'   => \$afmethod,
    'fm|FM=i'   => \$fstmethod,
	'sep|SEP=i' => \$separate,
    'pi=s'      => \$dopi,
    'dxy=s'     => \$dodxy,
    'da=s'      => \$doda,
    'tajd=s'    => \$dotajimad,
    'sl|SL=i'   => \$seqlen,
	'k|K=i'     => \$keeptmp,
	'ilg=i'		=> \$ignorelg,
    'h|help'    => \&usage
)
or &usage("\nERROR: Argument invalid (see above)\n\n");


# Check arguments
# ----------------------------------------------------------
my @errors=();
if (!defined($infile)){
	push (@errors, "bcf/vcf input file not specified");
}
else{
	push (@errors, "bcf/vcf input file does not exist or cannot be reached (you typed '-i $infile')")
		if (! -e $infile);
}
if (!defined($poppair)){
	push (@errors, "population pair input file not specified");
}
else{
	push (@errors, "population pair inputfile does not exist or cannot be reached (you typed '-p $poppair')")
		if (! -e $poppair);
}
if (!defined($outfile)){
	push (@errors, "output file not specified");
}
if (defined($windowsfile) && !defined($outfile2)){
	push (@errors, "output file for genomic windows not specified");
}
if ($prior !~ m/^[0|1]$/){
	push (@errors, "prior must be 0 (uniform) or 1 (Hardy-Weinberg) (you typed '-pr $prior')");
}
if ($mincov <= 0 || $mincov >= 1){
	push (@errors, "minimum sampling coverage must be a number >0 and <1 (you typed '-ms $mincov')");
}
if ($mincovp <= 0 || $mincovp >= 1){
	push (@errors, "minimum sampling coverage per population must be a number >0 and <1 (you typed '-msp $mincovp')");
}
if ($mmaf < 0 || $mmaf >=0.5){
	push (@errors, "minimum minor allele frequency be a number >=0 and <0.5 (you typed '-mf $mmaf')");
}
if ($afmethod !~ m/^[0|1]$/){
	push (@errors, "use EM algorithm to estimate AF must be 0 (no) or 1 (yes) (you typed '-am $afmethod')");
}
if ($fstmethod !~ m/^[0|1|2]$/){
	push (@errors, "method for Fst estimation must be 0 (Hudson), 1 (Hudson loci), or 2 (Weir-Cockerham) (you typed '-fm $fstmethod')");
}
if ($separate !~ m/^[0|1]$/){
	push (@errors, "estimate the AF separately for each population must be 0 (no) or 1 (yes) (you typed '-sep $separate')");
}
if ($dopi !~ m/^[0|1]$/){
	push (@errors, "calculate pi (nucleotide diversity) must be 0 (no) or 1 (yes) (you typed '-pi $dopi')");
}
if ($dodxy !~ m/^[0|1]$/){
	push (@errors, "calculate dxy must be 0 (no) or 1 (yes) (you typed '-dxy $dodxy')");
}
if ($doda !~ m/^[0|1]$/){
	push (@errors, "calculate da must be 0 (no) or 1 (yes) (you typed '-da $doda')");
}
if ($dotajimad !~ m/^[0|1]$/){
	push (@errors, "calculate da must be 0 (no) or 1 (yes) (you typed '-d $dotajimad')");
}
if ($seqlen !~ m/^[0-9]+$/){
	push (@errors, "sequence length must be a positive integer (you typed '-sl $seqlen')");
}
if ($keeptmp !~ m/^[0|1]$/){
	push (@errors, "keep the temporary files must be 0 (no) or 1 (yes) (you typed '-k $keeptmp')");
}
# if ($dotajimad == 1 && $seqlen==0){
# 	push (@errors, "sequence length must be provided for Tajima's D test");
# }

if (scalar(@errors)>0){
	my $out="\n\nThere were ERRORS:\n";
	foreach my $e (@errors){
		$out.="\t$e\n";
	}
	$out.="\n\n";
	&usage($out);
}
# ----------------------------------------------------------

# correct number of threads for compressing (as they are extra)
# ----------------------------------------------------------
$nthreads--;
# ----------------------------------------------------------

# Adjust options
# ----------------------------------------------------------
if ($afmethod==1 && $prior==0 && $separate==0){
	print "\nWarning: using uniform priors and EM algorithm implies AFs are estimated separately\n\n";
	$separate=1;
}
# ----------------------------------------------------------

# Absolute paths to files
# ----------------------------------------------------------
$infile=File::Spec->rel2abs($infile);
$poppair=File::Spec->rel2abs($poppair);
$windowsfile=File::Spec->rel2abs($windowsfile) if (defined($windowsfile));
$outfile=File::Spec->rel2abs($outfile);
$outfile2=File::Spec->rel2abs($outfile2) if (defined($windowsfile));
# ----------------------------------------------------------

# temporary files
# ----------------------------------------------------------
my $randpart=$ENV{USER}."-".time()."-".int(rand(1000000));
my $tmp=dirname($outfile)."/tmp_".$randpart;
my $sampfile="$tmp/samp_".$randpart;
my $scafile="$tmp/sca_".$randpart;
my $glfile="$tmp/gl_".$randpart; # genotype likelihood
my $afemfile="$tmp/afem_".$randpart;
my $bcffile="$tmp/bcf_".$randpart;
# my $estpEM="$tmp/estpEM_".$randpart;
# ----------------------------------------------------------


# Read population pair file and create file with samples
# ----------------------------------------------------------
my %sampop=();
my %samplesp=();
my @samples=();
my @pops=();
open (FILE, "$poppair")
	or die ("\nCan't open populations file $poppair\n\n");
	while (<FILE>){
		next if (/^(\#|\s*\n)/);
		my @aux=split(/\s|\,/,$_); 
		$aux[0]=~ s/(\.sorted)?\.[b|s]am//g; # remove any possible suffixes
		push (@samples, $aux[0]);
		push (@{$samplesp{$aux[1]}},$aux[0]);
		$sampop{$aux[0]}=$aux[1];
		push (@pops, $aux[1]) if (!grep(/^$aux[1]$/,@pops));
	}
close (FILE);
#----------------------------------------------------------

# Read windows file
# ----------------------------------------------------------
# Format must be scaffold\tlb\tub (with header)
my %windows;
my $i=0;
my $sca='';
if (defined($windowsfile)){
	open (FILE, "$windowsfile")
		or die ("\nCan't open genomic windows file $windowsfile\n\n");
		while (<FILE>){
			next if ($.==1 && !/[0-9]/); # header
			next if (/^(\#|\s*\n)/); # comments and empty lines
			my @aux=split(/\s|\,/,$_);
			$i=$#{$windows{$aux[0]}}+1;
			# print "line $_ -- $aux[0]-$i-lb $aux[1]\n"; # debug
			$windows{$aux[0]}[$i]{'lb'}=$aux[1];
			$windows{$aux[0]}[$i]{'ub'}=$aux[2];
			$i++;
		}
	close (FILE);
}
# ----------------------------------------------------------

# Calculate number of samples per population
my %nsamples=();
foreach my $s (keys %sampop){
	$nsamples{$sampop{$s}}++;;
}

# Exit if # pops > 2 in poppair file
die ("\nERROR: More than 2 populations in the population file $poppair\n\n")
	if (scalar(@pops)>2);

# Exit if # samples is <2 for any population
foreach my $p (@pops){
	die ("\nERROR: Fewer than 2 individuals in population $p\n\n")
		if ($nsamples{$p} < 2);
}
# ----------------------------------------------------------

# Exit if format is not vcf or bcf
# ----------------------------------------------------------
die ("\nERROR: $infile is not a bcf/vcf file\n\n")
	if ($infile !~ m/(\.bcf|\.vcf|\.vcf\.gz)$/);
# ----------------------------------------------------------

# Create temporary directory and file with samples
# ----------------------------------------------------------
make_path("$tmp");
open (TMP, ">$sampfile")
	or die ("\nCan't write to file $sampfile\n\n");
	foreach my $s (@samples){
		print TMP "$s\n";
	}
close (TMP);

foreach my $p (keys %samplesp){
	open (TMP, ">$sampfile.$p")
		or die ("\nCan't write to file $sampfile.$p\n\n");
		foreach my $s (@{$samplesp{$p}}){
			print TMP "$s\n";
		}
	close (TMP);
}
# ----------------------------------------------------------

# # Create estpEM executable from embedded data
# # ----------------------------------------------------------
# my $bin=pack 'H*', map s{\s+}{}gr, do { local $/; scalar <DATA> };
# open (EXEC, ">$estpEM");
# 	print EXEC "$bin";
# close(EXEC);
# system("chmod +x $estpEM");
# # ----------------------------------------------------------

# Filter bcf by population coverage separately and merge if need be
# -------------------------------------------------------------------
my $in=""; # command line to read data afterwards
if ($mincovp > 0){
	# Get list of SNVs with sampling coverage >= $mincovp
	foreach my $p (keys %nsamples){
		system("$bcftools view --threads $nthreads -O v $infile 2> /dev/null | \\
			   perl -pe 's/(\.sorted)?\.[b|s]am//g if (/^\#CHROM/)' | \\
			   $bcftools view -S $sampfile.$p -O v 2> /dev/null | \\
			   $bcftools filter -e '(AN/2)<($nsamples{$p}*$mincovp)' -O v | \\
			   $bcftools view -H -O v | cut -f1,2 > $bcffile.$p.sites 2> /dev/null");
	}
	# Get list of snvs shared in the two populations
	system("cat $bcffile.*.sites | sort | uniq -d > $bcffile.incsites");
	
	# data stream filtering also by global sampling coverage 
	$in=("$bcftools view --threads $nthreads -O v -R $bcffile.incsites $infile 2> /dev/null | \\
		  perl -pe 's/(\.sorted)?\.[b|s]am//g if (/^\#CHROM/)' | \\
		  $bcftools view -S $sampfile -O v | \\
		  $bcftools +fill-tags -O v -- -t AF 2> /dev/null | \\
		  $bcftools filter -O v -e '(AN/2)<(".scalar(@samples)."*$mincov)' -O v 2> /dev/null |");
}
else{
	# data stream filtering only by global sampling coverage
	$in=("$bcftools view --threads $nthreads -O v $infile 2> /dev/null| \\
	      perl -pe 's/(\.sorted)?\.[b|s]am//g if (/^\#CHROM/)' | \\
		  $bcftools view -S $sampfile -O v | \\
          $bcftools +fill-tags -O v -- -t AF 2> /dev/null | \\
	      $bcftools filter -e '(AN/2)<(".scalar(@samples)."*$mincov)' -O v 2> /dev/null |");
}
# -------------------------------------------------------------------

# Convert from bcf to gl 
# ----------------------------------------------------------
# print "IN: $in\n"; # debug

my ($x, $y)=&bcf2gl($in, $glfile);
my @bcf_af=@$x;
my @ids=@$y;


# Read glfile and split by population
# This is all assuming samples will be in the same order as in the samples file
if ($separate==1){
	my %ids_files;
	foreach my $i (0..$#pops){
		open (my $fh, ">$glfile.$i")
			or die ("\nCan't write to file $glfile.$i\n\n");
		$ids_files{$pops[$i]}=$fh;
		print {$ids_files{$pops[$i]}} "$nsamples{$pops[$i]} ".scalar(@bcf_af)."\n";
	}

	open (GL, $glfile)
		or die "\nCan't read $glfile\n\n";
		# my $header=<GL>;
		<GL>;

		while (<GL>){
			chomp;
			my @aux=split(/\s+/,$_);
			# print SNP id 
			my $snpid=shift(@aux);

			foreach my $p (@pops){
				print {$ids_files{$p}} "$snpid";
			}
			# print the genotype likelihoods to the corresponding files
			# OJO GL van de 3 en 3
			my $i=0;
			while ($i < $#aux){
				print {$ids_files{$sampop{$samples[$i/3]}}} " $aux[$i] $aux[$i+1] $aux[$i+2]";
				$i+=3;
			}
			# print linebreak
			foreach my $p (@pops){
				print {$ids_files{$p}} "\n";
			}
		}
	close (GL);
		
	foreach my $p (keys %ids_files){
		close ($ids_files{$p});
	}
}
# ----------------------------------------------------------

# Calculate AFs with EM algorithm
# ----------------------------------------------------------
my %afl=(); # AF for each locus
if ($afmethod==1){ # use EM algorithm to estimate allele frequencies
	print "Estimating allele frequencies by maximum likelihood...\n\n";
	
	# Note: global AF estimation is always needed to filter out by MAF
	system("$estpEM -i $glfile -o $afemfile -m 20 -e 0.001 >& $afemfile.log"); 
	open (FILE, "$afemfile")
		or die ("\nCan't open file $afemfile\n\n");
		while (<FILE>){
			chomp;
			my @aux=split(/\s+/,$_);
			$afl{$aux[0]}=$aux[2];
		}
	close (FILE);
	
	# if ($separate==1 && $dopi==1){
	if ($separate==1){
		system("$estpEM -i $glfile.0 -o $afemfile.0 -m 20 -e 0.001 >& $afemfile.0.log"); 
		system("$estpEM -i $glfile.1 -o $afemfile.1 -m 20 -e 0.001 >& $afemfile.1.log");

		open (FILE, "$afemfile.0")
			or die ("\nCan't open file $afemfile.0\n\n");
			while (<FILE>){
				chomp;
				my @aux=split(/\s+/,$_);
				$afl{$pops[0]}{$aux[0]}=$aux[2];
			}
		close (FILE);

		open (FILE, "$afemfile.1")
			or die ("\nCan't open file $afemfile.1\n\n");
			while (<FILE>){
				chomp;
				my @aux=split(/\s+/,$_);
				$afl{$pops[1]}{$aux[0]}=$aux[2];
			}
		close (FILE);
	}

	print "Done\n\n";
}
# ----------------------------------------------------------

# Calculate Fst for each locus
# =============================================================================
my $nsnps=0;
my @statsloci=(); # loci statistics: fst,pi,dxy,da
my @nums=(); # numerators for calculating mean Fst
my %numsloci=(); # numerators for calculating mean Fst per window
my @denoms=(); # denominator for calculating mean Fsts
my %denomsloci=(); # denominators for calculating mean Fst per window
my %pi1loci=(); # pop1 pi
my %pi2loci=(); # pop2 pi
my %piloci=(); # pooled pops pi
my %dxyloci=(); # dxy
my %daloci=(); # da
my %tajimad1loci=(); # pop1 Tajima D
my %tajimad2loci=(); # pop2 Tajima D
my %tajimadloci=(); # pooled pops Tajima D

my %afsp1=(); # pop1 allele frequencies for pi, dxy, da
my %afsp2=(); # pop1 allele frequencies for pi, dxy, da
my %afspavg=(); # pooled pops allele frequencies for pi, dxy, da

open (GL, $glfile)
	or die "\nCan't read $glfile\n\n";
		# my $header=<GL>;
	<GL>;

	# foreach my $i (0..$#genolhl){ # For each locus (i.e. SNP)
	while (<GL>){
		chomp;
		my $genolhl=$_;
		print "SNPs processed (used for Fst estimation): $nsnps\r" if ($nsnps%100==0);
		$|++; # flush buffer
		
		my @aux=split(/\s+|\:/,$genolhl);
		# my $snpid=shift(@aux);
		my $sca=shift(@aux);
		my $pos=shift(@aux);
		my $snpid=$sca.':'.$pos;

		# Calculate AFs
		# ------------------------------------------------------------------
		my $af='';
		my %afs=(); # allele frequencies of each population
		my %genotypes=(); # required when using uniform prior
		if ($afmethod==0 && $prior==0){ # use empirical Bayesian genotypes with uniform prior
			# Estimate allele frequency using genotype probabilities 
			# with uniform priors (doesn't require knowing AFs beforehand)
			# --------------------------------------------------------
			# Note: @aux only contains genotype likelihoods, first position storing snp id was removed before
			my $j=0;
			while ($j < $#aux){
				# print "$snpid - $ids[$j/3] - ".$sampop{$ids[$j/3]}." - ".$aux[$j].",".$aux[$j+1].",".$aux[$j+2]."\n";
				push (@{$genotypes{$sampop{$ids[$j/3]}}},&pl2cgp_uni($aux[$j],$aux[$j+1],$aux[$j+2]));
				$j+=3;
			}

			# Calculate allele frequencies for each population
			# ------------------------------------------------
			foreach my $p (keys %genotypes){
				$afs{$p}=sum(@{$genotypes{$p}})/(2*scalar(@{$genotypes{$p}}))
			}
			# ------------------------------------------------

			# AF of 2 pops together (required only for MAF filter)
			$af=sum((@{$genotypes{$pops[0]}},@{$genotypes{$pops[1]}}))/(2*(scalar(@{$genotypes{$pops[0]}})+scalar(@{$genotypes{$pops[1]}})));
		}
		elsif ($afmethod==0 && $prior==1){ # Hardy-Weinberg priors using AFs coming from bcf file
			$af=$bcf_af[$nsnps]; # AF from bcftools
			if ($separate==0){
				# Estimate genotypes, which will be used for estimating
				# AFs for each population
				# Note: @aux only contains genotype likelihoods, first position storing snp id was removed before
				my $j=0;
				while ($j < $#aux){
					my @gt=($aux[$j],$aux[$j+1],$aux[$j+2]);
					push (@{$genotypes{$sampop{$ids[$j/3]}}},&pl2cgp_hw(\@gt,\$af));
					$j+=3;
				}
				# Calculate allele frequencies for each population
				# ------------------------------------------------
				foreach my $p (keys %genotypes){
					$afs{$p}=sum(@{$genotypes{$p}})/(2*scalar(@{$genotypes{$p}}))
				}
				# ------------------------------------------------
			}
			else{
				# ToDo - extract population specific AF from bcf file
				# $afs{$pops[0]}=$afl{$pops[0]}{$snpid};	
				# $afs{$pops[1]}=$afl{$pops[1]}{$snpid};
			}
		}
		elsif ($afmethod==1 && $prior==0){ # uniform prior implies AFs must be separated
			$af=$afl{$snpid};
			$afs{$pops[0]}=$afl{$pops[0]}{$snpid};
			$afs{$pops[1]}=$afl{$pops[1]}{$snpid};
		}
		elsif ($afmethod==1 && $prior==1){ # Hardy-Weinberg prior using AFs estimated with EM algorithm
			$af=$afl{$snpid};
			# print "AF EM: $snpid - $af\n"; # debug
			if ($separate==0){
				# Estimate genotypes, which will be used for estimating
				# AFs for each population
				# Note: @aux only contains genotype likelihoods, first position storing snp id was removed before
				my $j=0;
				while ($j < $#aux){
					my @gt=($aux[$j],$aux[$j+1],$aux[$j+2]);
					push (@{$genotypes{$sampop{$ids[$j/3]}}},&pl2cgp_hw(\@gt,\$af));
					$j+=3;
				}
				# Calculate allele frequencies for each population
				# ------------------------------------------------
				foreach my $p (keys %genotypes){
					$afs{$p}=sum(@{$genotypes{$p}})/(2*scalar(@{$genotypes{$p}}));
					# print "\t$p - $afs{$p}\n";
				}
				# ------------------------------------------------
			}
			else{
				# Estimate genotypes, which will be used for estimating
				# AFs for each population
				# Note: @aux only contains genotype likelihoods, first position storing snp id was removed before
				my $j=0;
				while ($j < $#aux){
					my @gt=($aux[$j],$aux[$j+1],$aux[$j+2]);
					push (@{$genotypes{$sampop{$ids[$j/3]}}},&pl2cgp_hw(\@gt,\$afl{$sampop{$ids[$j/3]}}{$snpid}));
					$j+=3;
				}
				# Calculate allele frequencies for each population
				# ------------------------------------------------
				foreach my $p (keys %genotypes){
					$afs{$p}=sum(@{$genotypes{$p}})/(2*scalar(@{$genotypes{$p}}));
				}
				# ------------------------------------------------
			}
		}
		# ------------------------------------------------------------------

		# minor allele frequency (for filtering)
		# --------------------------------------------------------
		my $maf=$af;
		$maf=1-$af if ($af>0.5);
		# -------------------------------------------------------

		if ($af < 1 &&     # exclude 'invariant' SNPs: all individuals have the same allele, but it is different from reference
			$maf > $mmaf){ # exclude SNPs with MAF < $mmaf

			# Estimate dAF and Fst for each locus
			# ------------------------------------------------
			# print "$snpid\t$pops[0]\n" if (!exists($afs{$pops[0]})); # debug
			# print "$snpid\t$pops[1]\n" if (!exists($afs{$pops[1]})); # debug

			my $dAF=sprintf("%.10f", $afs{$pops[0]}-$afs{$pops[1]});
			my ($fst,$num,$den);
			if ($fstmethod == 0){
				($fst,$num,$den)=HudsonFst($afs{$pops[0]},$afs{$pops[1]});
			}
			elsif ($fstmethod == 1){ # Hudson loci
				($fst,$num,$den)=HudsonFstloci($afs{$pops[0]},$afs{$pops[1]},$nsamples{$pops[0]},$nsamples{$pops[1]});
			}
			elsif ($fstmethod == 2){ # Weir & Cockerham
				($fst,$num,$den)=WeirCockerhamFst($afs{$pops[0]},$afs{$pops[1]},$nsamples{$pops[0]},$nsamples{$pops[1]});
			}

			$fst=sprintf("%.10f", $fst);
			
			my $locusline="$sca\t$pos\t$dAF\t$fst";
			# print "SNP $aux[0]-$aux[1] - AF1 $afs{$pops[0]} - AF2 $afs{$pops[1]}- FST $fst\n";
			# ------------------------------------------------
			# ------------------------------------------------

			# Store numerators and denominators for estimating
			# genome-wide Fst, and window Fst at the end
			# ------------------------------------------------
			push (@nums, $num);
			$sca=~ s/.*scaf//g if ($ignorelg==1);
			$numsloci{$sca}{$pos}=$num;
			push (@denoms, $den);
			$denomsloci{$sca}{$pos}=$den;
			# ------------------------------------------------
		
			# Store allele frequencies
			# ------------------------------------------------
			$afsp1{$sca}{$pos}=$afs{$pops[0]};
			$afsp2{$sca}{$pos}=$afs{$pops[1]};
			$afspavg{$sca}{$pos}=mean($afs{$pops[0]},$afs{$pops[1]});
			# ------------------------------------------------

			# Calculate loci pi
			# ------------------------------------------------
			if ($dopi==1){
				$pi1loci{$sca}{$pos}=estpi(1, $afs{$pops[0]});
				$pi2loci{$sca}{$pos}=estpi(1, $afs{$pops[1]});
				$piloci{$sca}{$pos}=estpi(1, $afspavg{$sca}{$pos});# use average allele frequency
				$locusline.="\t".sprintf("%.10f",$pi1loci{$sca}{$pos});
				$locusline.="\t".sprintf("%.10f",$pi2loci{$sca}{$pos});
				$locusline.="\t".sprintf("%.10f",$piloci{$sca}{$pos});
			}
			# ------------------------------------------------
			# Calculate loci dxy
			# ------------------------------------------------
			if ($dodxy==1){
				$dxyloci{$sca}{$pos}=estdxy(1,\$afs{$pops[0]},\$afs{$pops[1]});
				$locusline.="\t".sprintf("%.10f",$dxyloci{$sca}{$pos});
			}
			# ------------------------------------------------
			# Calculate loci da	
			# ------------------------------------------------
			if ($doda==1){
				my $dxy;
				if (defined($dxyloci{$sca}{$pos})){
					$dxy=$dxyloci{$sca}{$pos};
				}
				else{ 
					$dxy=estdxy(1,\$afs{$pops[0]},\$afs{$pops[1]}); 
				}
				my $pi1; my $pi2;
				if (defined($pi1loci{$sca}{$pos})){
					$pi1=$pi1loci{$sca}{$pos};
					$pi2=$pi2loci{$sca}{$pos};
				}
				else{ 
					$pi1=estpi(1,$afs{$pops[0]}); 
					$pi2=estpi(1,$afs{$pops[1]}); 
				}

				$daloci{$sca}{$pos}=estda($dxy,$pi1,$pi2);
				$locusline.="\t".sprintf("%.10f",$daloci{$sca}{$pos});
			}
			# ------------------------------------------------


			# Calculate Tajima D
			# ------------------------------------------------
			if ($dotajimad==1){
				$tajimad1loci{$sca}{$pos}=estTajimaD($nsamples{$pops[0]},$afs{$pops[0]});
				$tajimad2loci{$sca}{$pos}=estTajimaD($nsamples{$pops[1]},$afs{$pops[1]});
				$tajimadloci{$sca}{$pos}=estTajimaD(($nsamples{$pops[0]}+$nsamples{$pops[1]}),($afs{$pops[0]},$afs{$pops[1]}));
				$locusline.="\t".sprintf("%.10f",$tajimad1loci{$sca}{$pos});
				$locusline.="\t".sprintf("%.10f",$tajimad2loci{$sca}{$pos});
				$locusline.="\t".sprintf("%.10f",$tajimadloci{$sca}{$pos});
			}
			# ------------------------------------------------
			
			push (@statsloci, "$locusline");
			
			$nsnps++;
		}
	}
close (GL);

print "SNPs processed (used for Fst estimation): $nsnps\r";
print "\n\n";

# Genome-wide statistics
# ----------------------------------------------------------
# ----------------------------------------------------------

# Estimate genome-wide Fst
# ----------------------------------------------------------
my $meanfst=0;
if ($fstmethod==0 || $fstmethod==2){
	$meanfst=sprintf("%.10f", (1-mean(@nums)/mean(@denoms)));
}
elsif ($fstmethod==1){
	$meanfst=sprintf("%.10f", (mean(@nums)/mean(@denoms)));
}
my $header="pop1\tpop2\tmean".$fstmethod_names{$fstmethod};
my $line="$pops[0]\t$pops[1]\t$meanfst";
# print "pop1\tpop2\tmean".$fstmethod_names{$fstmethod}."\n";
# print "$pops[0]\t$pops[1]\t$meanfst\n\n";
# ----------------------------------------------------------

# Estimate global (i.e. genome-wide) pi
# ----------------------------------------------------------
if ($dopi==1){
	my $globalpi1=my $globalpi2=my $globalpi=0;

	my $sl=scalar(dumpHoH(\%afsp1)); # use no snps as seq length
	
	$globalpi1=estpi($sl, dumpHoH(\%afsp1));
	$globalpi2=estpi($sl, dumpHoH(\%afsp2));
	$globalpi=estpi($sl, dumpHoH(\%afspavg));
	$header.="\tpi1\tpi2\tpi";
	$line.="\t".sprintf("%.10f",$globalpi1);
	$line.="\t".sprintf("%.10f",$globalpi2);
	$line.="\t".sprintf("%.10f",$globalpi);

	if ($seqlen > 0){
		my $globalpi1site=estpi($seqlen, dumpHoH(\%afsp1));
		my $globalpi2site=estpi($seqlen, dumpHoH(\%afsp2));
		my $globalpisite=estpi($seqlen, dumpHoH(\%afspavg));
		$header.="\tpi1_site\tpi2_site\tpi_site";
		$line.="\t".sprintf("%.10f",$globalpi1site);
		$line.="\t".sprintf("%.10f",$globalpi2site);
		$line.="\t".sprintf("%.10f",$globalpisite);
	}
}
# ----------------------------------------------------------

# Estimate global (i.e. genome-wide) dxy
# ----------------------------------------------------------
if ($dodxy==1){
	my $sl=scalar(dumpHoH(\%afsp1)); # use no snps as seq length
	my $globaldxy=estdxy($sl,\dumpHoH(\%afsp1),\dumpHoH(\%afsp2));
	$header.="\tdxy";
	$line.="\t".sprintf("%.10f",$globaldxy);

	if ($seqlen > 0){
		my $globaldxysite=estdxy($seqlen,\dumpHoH(\%afsp1),\dumpHoH(\%afsp2));
		$header.="\tdxy_site";
		$line.="\t".sprintf("%.10f",$globaldxysite);
	}
}
# ----------------------------------------------------------

# Estimate global (i.e. genome-wide) da
# ----------------------------------------------------------
if ($doda==1){
	my $sl=scalar(dumpHoH(\%afsp1)); # use no snps as seq length
	my $globaldxy=estdxy($sl,\dumpHoH(\%afsp1),\dumpHoH(\%afsp2));
	my $globalpi1=estpi($sl, dumpHoH(\%afsp1));
	my $globalpi2=estpi($sl, dumpHoH(\%afsp2));

	my $globalda=estda($globaldxy,$globalpi1,$globalpi2);
	$header.="\tda";
	$line.="\t".sprintf("%.10f",$globalda);

	if ($seqlen > 0){
		my $globaldxysite=estdxy($seqlen, \dumpHoH(\%afsp1),\dumpHoH(\%afsp2));
		my $globalpi1site=estpi($seqlen, dumpHoH(\%afsp1));
		my $globalpi2site=estpi($seqlen, dumpHoH(\%afsp2));
	
		my $globaldasite=estda($globaldxysite,$globalpi1site,$globalpi2site);

		$header.="\tda_site";
		$line.="\t".sprintf("%.10f",$globaldasite);
	}
}
# ----------------------------------------------------------

# Estimate global (i.e. genome-wide) Tajima D
# ----------------------------------------------------------
if ($dotajimad==1){
	$header.="\ttajimad1\ttajimad2\ttajimad";
	if ($seqlen>0){
		my $globaltajimad1=estTajimaD($nsamples{$pops[0]}, dumpHoH(\%afsp1));
		my $globaltajimad2=estTajimaD($nsamples{$pops[1]}, dumpHoH(\%afsp2));
		my $globaltajimad=estTajimaD(($nsamples{$pops[0]}+$nsamples{$pops[1]}),(dumpHoH(\%afsp1),dumpHoH(\%afsp2)));
		$line.="\t".sprintf("%.10f",$globaltajimad1);
		$line.="\t".sprintf("%.10f",$globaltajimad2);
		$line.="\t".sprintf("%.10f",$globaltajimad);
	}
	else{
		# print warning
		$line.="\tNA\tNA\tNA";
	}
}
# ----------------------------------------------------------

print "$header\n";
print "$line\n";
# ----------------------------------------------------------
# ----------------------------------------------------------

# Sort loci by linkage group -> scaffold -> position
# ----------------------------------------------------------
# Remove NAs with a number for numerical sorting
if ($statsloci[0]=~ /^lg/){
	foreach my $fl (@statsloci){
		$fl=~ s/NA/999999999/g;
	}
	@statsloci=sort {
			my ($alg) = $a =~ /lg([0-9]+)/;
			my ($blg) = $b =~ /lg([0-9]+)/;
			my ($aord) = $a =~ /ord([0-9]+)/;
			my ($bord) = $b =~ /ord([0-9]+)/;
			my ($ascaf) = $a =~ /scaf*([0-9]+)/;
			my ($bscaf) = $b =~ /scaf*([0-9]+)/;
			$alg <=> $blg || $aord <=> $bord || $ascaf <=> $bscaf;
			} @statsloci;
	# Put back NAs
	foreach my $fl (@statsloci){
		$fl=~ s/999999999/NA/g;
	}
}
elsif ($statsloci[0]=~ /^scaffold/){
	@statsloci=sort {
			my ($ascaf) = $a =~ /scaffold_([0-9]+)/;
			my ($bscaf) = $b =~ /scaffold_([0-9]+)/;
			$ascaf <=> $bscaf;
			} @statsloci;
}
elsif ($statsloci[0]=~ /^pseudoscaff/){
	my @sortstatsloci=sort {
		my ($ascaf) = $a =~/pseudoscaff_([0-9]+)/;
		my ($bscaf) = $b =~/pseudoscaff_([0-9]+)/;
		$ascaf <=> $bscaf;
	} @statsloci;
}
# ----------------------------------------------------------


# window statistics
# ----------------------------------------------------------
# ----------------------------------------------------------

# Estimate windows Fst
# ----------------------------------------------------------
my @fstwindows=();
if (defined($windowsfile)){
	foreach my $sca (keys %windows){
		foreach my $i (0..$#{$windows{$sca}}){
			my $lb=$windows{$sca}[$i]{'lb'};
			my $ub=$windows{$sca}[$i]{'ub'};
			my $wsize=$ub-$lb;
			my @nums=();
			my @denoms=();
			my @af1=();
			my @af2=();
			# print "window $sca: $lb-$ub\n";
			my $meanfst=my $meanpi1=my $meanpi2=my $meandxy=my $meanda='NA';
			my $statsline="$sca\t$lb\t$ub";

			# Get allele frequencies for SNPs in the window
			# ----------------------------------------------------------
			foreach my $pos (sort {$a <=> $b} keys %{$pi1loci{$sca}}){
				last if ($pos > $ub);
				if ($pos > $lb && $pos < $ub){
					push (@af1,$pi1loci{$sca}{$pos});
					push (@af2,$pi2loci{$sca}{$pos});
				}
			}
			# ----------------------------------------------------------

			if (scalar(@af1)>0){ # there are no snps within this window
				# Estimate windows fst
				# ----------------------------------------------------------
				if (defined ($numsloci{$sca})){
					foreach my $pos (sort {$a <=> $b} keys %{$numsloci{$sca}}){
						last if ($pos > $ub);
						if ($pos > $lb && $pos < $ub){
							push (@nums, $numsloci{$sca}{$pos});
							push (@denoms, $denomsloci{$sca}{$pos});
						}
					}
					if (scalar(@nums)>0){
						if ($fstmethod==0 || $fstmethod==2){
							$meanfst=sprintf("%.10f", (1-mean(@nums)/mean(@denoms)));
						}
						elsif ($fstmethod==1){
							$meanfst=sprintf("%.10f", (mean(@nums)/mean(@denoms)));
						}
					}
					my $nsnps=scalar(@nums);
					$statsline.="\t$nsnps\t$meanfst";
				}
				# ----------------------------------------------------------
				
				# Estimate windows pi
				# ----------------------------------------------------------
				if ($dopi==1){
					my $windowpi1=estpi($wsize, @af1);
					my $windowpi2=estpi($wsize, @af2);
					my $windowpi=estpi($wsize, (@af1,@af2));

					$statsline.="\t".sprintf("%.10f", $windowpi1);
					$statsline.="\t".sprintf("%.10f", $windowpi2);
					$statsline.="\t".sprintf("%.10f", $windowpi);
				}
				# ----------------------------------------------------------

				# Estimate windows dxy
				# ----------------------------------------------------------
				if ($dodxy==1){
					my $windowdxy=estdxy($wsize,\@af1,\@af2);
					$statsline.="\t".sprintf("%.10f", $windowdxy);
				}		
				# ----------------------------------------------------------
				
				# Estimate windows da
				# ----------------------------------------------------------
				if ($doda==1){
					my $windowdxy=estdxy($wsize,\@af1,\@af2);
					my $windowpi1=estpi($wsize, @af1);
					my $windowpi2=estpi($wsize, @af2);
					
					my $windowda=estda($windowdxy,$windowpi1,$windowpi2);
					$statsline.="\t".sprintf("%.10f", $windowda);
				}
				# ----------------------------------------------------------

				# Estimate windows Tajima's D
				# ----------------------------------------------------------
				if ($dotajimad==0){
					my $windowtajimad1=estTajimaD($nsamples{$pops[0]},@af1);
					my $windowtajimad2=estTajimaD($nsamples{$pops[1]},@af2);
					my $windowtajimad=estTajimaD(($nsamples{$pops[0]}+$nsamples{$pops[1]}),@af1,@af2);

					$statsline.="\t".sprintf("%.10f", $windowtajimad1);
					$statsline.="\t".sprintf("%.10f", $windowtajimad2);
					$statsline.="\t".sprintf("%.10f", $windowtajimad);
				}
				# ----------------------------------------------------------
			}
			else{
				# print STDERR "WARNING: No data for window $sca:$lb-$ub \n";
				$statsline.="\t0\tNA";
				$statsline.="\tNA\tNA\tNA" if ($dopi==1);
				$statsline.="\tNA" if ($dodxy==1);
				$statsline.="\tNA" if ($doda==1);
				$statsline.="\tNA\tNA\tNA" if ($dotajimad==1);
			}
			# push (@fstwindows,"$sca\t$lb\t$ub\t$nsnps\t$meanfst");
			push (@fstwindows,"$statsline");
		}
	}
}
# ----------------------------------------------------------

# ----------------------------------------------------------
# ----------------------------------------------------------

# Sort windows by linkage group -> scaffold -> position
# ----------------------------------------------------------
# Remove NAs with a number for numerical sorting
if (defined($windowsfile)){
	if ($fstwindows[0]=~ /^lg/){
		foreach my $fl (@fstwindows){
			$fl=~ s/NA/999999999/g;
		}
		@fstwindows=sort {
				my ($alg) = $a =~ /lg([0-9]+)/;
				my ($blg) = $b =~ /lg([0-9]+)/;
				my ($aord) = $a =~ /ord([0-9]+)/;
				my ($bord) = $b =~ /ord([0-9]+)/;
				my ($ascaf) = $a =~ /scaf*([0-9]+)/;
				my ($bscaf) = $b =~ /scaf*([0-9]+)/;
				$alg <=> $blg || $aord <=> $bord || $ascaf <=> $bscaf;
				} @fstwindows;
		# Put back NAs
		foreach my $fl (@fstwindows){
			$fl=~ s/999999999/NA/g;
		}
	}
	elsif ($fstwindows[0]=~ /^scaffold/){
		@fstwindows=sort {
				my ($ascaf) = $a =~ /scaffold_([0-9]+)/;
				my ($bscaf) = $b =~ /scaffold_([0-9]+)/;
				$ascaf <=> $bscaf;
				} @fstwindows;
	}
	elsif ($fstwindows[0]=~ /^pseudoscaff/){
		my @sortfstwindows=sort {
			my ($ascaf) = $a =~/pseudoscaff_([0-9]+)/;
			my ($bscaf) = $b =~/pseudoscaff_([0-9]+)/;
			$ascaf <=> $bscaf;
		} @fstwindows;
	}
}
# ----------------------------------------------------------

# =============================================================================

# Write output to files
# =============================================================================
open (FILE, ">$outfile")
	or die ("\nCan't write to file $outfile\n\n");
	$header="scaffold\tposition\tdAF\t".$fstmethod_names{$fstmethod};
	$header.="\tpi1\tpi2\tpi" if ($dopi==1);
	$header.="\tdxy" if ($dodxy==1);
	$header.="\tda" if ($doda==1);
	$header.="\ttajimaD1\ttajimaD2\ttajimaD" if ($dotajimad==1);

	print FILE "$header\n";

	foreach my $l (@statsloci){
		print FILE "$l\n";
	}
close (FILE);

if (defined($windowsfile)){
	open (FILE, ">$outfile2")
		or die ("\nCan't write to file $outfile\n\n");
		$header="scaffold\tlb\tub\tnloci\t".$fstmethod_names{$fstmethod};
		$header.="\tpi1\tpi2\tpi" if ($dopi==1);
		$header.="\tdxy" if ($dodxy==1);
		$header.="\tda" if ($doda==1);
		$header.="\ttajimaD1\ttajimaD2\ttajimaD" if ($dotajimad==1);

		print FILE "$header\n";

		foreach my $w (@fstwindows){
			print FILE "$w\n";
		}
	close (FILE);
}
# =============================================================================

# Remove temporary files
remove_tree($tmp) if ($keeptmp==0);


# debug
# ................................................
# Record amount in use afterwards
# $mu->record('after something_memory_intensive()');
#
# # Spit out a report
# $mu->dump();
# ................................................
# ==============================================================================
# ==============================================================================
# ============================== SUBROUTINES ===================================
# ==============================================================================
# ==============================================================================


# Convert Phred-scaled genotype likelihoods (pl) to
# composite genotype probabilities (cgp)
# ==============================================================================

# using uniform priors 
sub pl2cgp_uni{
	my @pls=@_;
	
	my $sum=0;
	my @gps=();
	foreach my $pl (@pls){
		my $gl=10**(-$pl/10);
		push (@gps, $gl);
		$sum+=$gl;
	}
	foreach my $pr (@gps){
		$pr=$pr/$sum;
	}
	my $cgp=$gps[2]*2+$gps[1]; # alternate allele dosage

	return($cgp);
}
# -----------------------------------------------------------------------------

# using Hardy-Weinberg priors
# -----------------------------------------------------------------------------
sub pl2cgp_hw{
	my ($a,$b)=@_;
	my @pls=@$a;
	my $af=$$b;

	# if ref allele is not the major allele
	# likelihoods must be swapped
	# minor allele frequency is then 1-af
	# if ($af>0.5){
	# 	my @aux=@pls;
	# 	$pls[0]=$aux[2];
	# 	$pls[2]=$aux[0];
	# 	$af=1-$af;
	# }

	my @gps=();
	# 1st term => Phred-descalation
	# 2nd term => Hardy-Weinberg prior
	$gps[0]=(10**(-$pls[0]/10)) * ((1-$af)**2);
	$gps[1]=(10**(-$pls[1]/10)) * (2*$af*(1-$af));
	$gps[2]=(10**(-$pls[2]/10)) * ($af**2);

	my $sum=0;
	foreach my $pr (@gps){
		$sum+=$pr
	}
	foreach my $pr (@gps){
		$pr=$pr/$sum;
	}
	my $cgp=$gps[2]*2+$gps[1]; # alternate allele dosage 

	return($cgp);
}
# -----------------------------------------------------------------------------
# ==============================================================================

# Calculate Fst
# ==============================================================================
# original Hudson Fst=1-Hw/Hb, see SI in Bhatia et al. 2013, doi:10.1101/gr.154831.113
sub HudsonFst{
	my $af1=shift;
	my $af2=shift;
	my $numerator=$af1*(1-$af1)+$af2*(1-$af2); # Hw
	my $denominator=$af1*(1-$af2)+$af2*(1-$af1); # Hb
	my $fst=(1-$numerator/$denominator);
	return($fst,$numerator,$denominator);
}

# equation 10 in Bhatia et al. 2013, doi:10.1101/gr.154831.113
sub HudsonFstloci{
	my $af1=shift;
	my $af2=shift;
	my $n1=shift;
	my $n2=shift;
	my $numerator=($af1-$af2)**2-(($af1*(1-$af1))/($n1-1))-(($af2*(1-$af2))/($n2-1));
	my $denominator=$af1*(1-$af2)+$af2*(1-$af1);
	my $fst=$numerator/$denominator;
	return($fst,$numerator,$denominator);
}

# equation 6 in Bhatia et al. 2013, doi:10.1101/gr.154831.113
sub WeirCockerhamFst{
	my $af1=shift;
	my $af2=shift;
	my $n1=shift;
	my $n2=shift;
	my $numerator=(2*(($n1*$n2)/($n1+$n2))*(1/($n1+$n2-2)))*($n1*$af1*(1-$af1)+$n2*$af2*(1-$af2));
	my $denominator=(($n1*$n2)/($n1+$n2))*($af1-$af2)**2+(2*(($n1*$n2)/($n1+$n2))-1)*(1/($n1+$n2-2))*($n1*$af1*(1-$af1)+$n2*$af2*(1-$af2));
	my $fst=(1-$numerator/$denominator);
	return($fst,$numerator,$denominator);
}
# ==============================================================================

# Calculate pi (nucleotide diversity)
# ==============================================================================
# sub estpi{
# 	my $af=shift;
# 	my $pi=2*$af*(1-$af);
# 	return ($pi)
# }
sub estpi{
	my $n=shift; # sequence length
	my @afs=@_; # allele frequencies
	my $pi=0;
	foreach my $af (@afs){
		$pi+=2*$af*(1-$af);
	}
	$pi/=$n;
	return ($pi)
}
# ==============================================================================

# Calculate dxy
# ==============================================================================
# sub estdxy{
# 	my $af1=shift;
# 	my $af2=shift;
# 	my $dxy=$af1*(1-$af2)+(1-$af1)*$af2;
# 	return($dxy);
# }
sub estdxy{
	my $n=shift; # sequence length
	my $a=shift; # allele frequencies pop 1
	my $b=shift; # allele frequencies pop 2

	my @afs1=my @afs2=();
	if (ref($a) eq 'ARRAY'){
		@afs1=@$a;
		@afs2=@$b;
	}
	else{
		push (@afs1,$$a);
		push (@afs2,$$b);
	}

	my $dxy=0;
	foreach my $i (0..$#afs1){
		$dxy+=$afs1[$i]*(1-$afs2[$i])+(1-$afs1[$i])*$afs2[$i];
	}
	$dxy/=$n;

	return($dxy);
}
# ==============================================================================

# Calculate da
# ==============================================================================
sub estda{
	my $dxy=shift;
	my $pi1=shift;
	my $pi2=shift;

	my $da=$dxy-($pi1+$pi2)/2;
	return($da);
}

# ==============================================================================

# Calculate Tajima's D
# Derived partially from vcftools code
# Based on Tajima 1989 and Carlston et al 2005 doi:10.1101/gr.4326505
# ==============================================================================
sub estTajimaD{
	my $nsamp=shift; # number of samples=individuals
	my @afs=@_; # allele frequencies
	
	my $nsnps=scalar(@afs); # no segregating sites
	my $n=2*$nsamp; # number of alleles

	# estimate pi (eq 2 in Carlston et al 2005)
	my $pisum=0;
	foreach my $af (@afs){
		$pisum+=2*$af*(1-$af);
	}
	my $pi=$pisum*$n/($n-1); # pi

	# estimate var(pi-thetas)
	# multiple equations from Tajima 1989
	my $a1=my $a2=0;
	foreach my $i (1..$n){
		$a1+=1/$i;
		$a2+=1/$i**2;
	}
	my $b1=($n+1)/(3*($n-1));
	# print "NSAMP: $nsamp - N: $n - A1:$a1\n";
	my $c1=$b1-(1/$a1);
	my $e1=$c1/$a1;
	my $b2=2*($n**2+$n+3)/(9*$n*($n-1));
	my $c2=$b2-($n+2)/($a1*$n)+($a2/$a2**2);
	my $e2=$c2/(($a1**2)+$a2);

	my $var=$e1*$nsnps+$e2*$nsnps*($nsnps-1);

	# estimate theta s (eq 1 in Carlston et al 2005)
	my $thetaS=$nsnps/$a1; # thetas

	# Tajima's D
	my $tajd=($pi-$thetaS)/sqrt($var);

	return($tajd);
}
# ==============================================================================

# Basic statistics
# ==============================================================================
sub mean{
	my $avg=0;
	foreach my $v (@_){
		$avg+=$v;
	}
	$avg/=scalar(@_);

	return($avg);
}

sub sum{
	my $s=0;
	foreach (@_){
		$s+=$_;
	}

	return($s)
}

# mean of variable stored in hash
# with scaffolds and positions
sub meanscapos{
	my $a=shift;
	my %var=%$a;
	my $meanvar=0;
	foreach my $sca (keys %var){
		foreach my $pos (keys %{$var{$sca}}){
			$meanvar+=$var{$sca}{$pos};
			$nsnps++;
		}
	}
	$meanvar/=$nsnps;

	return($meanvar);
}
sub meanscapos2var{
	my $a=shift;
	my $b=shift;
	my %var1=%$a;
	my %var2=%$b;
	my $meanvar=0;
	foreach my $sca (keys %var1){
		foreach my $pos (keys %{$var1{$sca}}){
			$meanvar+=($var1{$sca}{$pos}+$var2{$sca}{$pos})/2;
			$nsnps++;
		}
	}
	$meanvar/=$nsnps;

	return($meanvar);
}

# dump hash of hashes
sub dumpHoH{
	my $a=shift;
	my %var=%$a;
	my @out=();
	foreach my $sca (keys %var){
		foreach my $pos (keys %{$var{$sca}}){
			push (@out,$var{$sca}{$pos});
		}
	}
	return (@out);
}

# ==============================================================================

# Convert bcf to gl
# ==============================================================================
sub bcf2gl{
	my $in=shift;
	my $out=shift;
	my @af_bcf=();
	my @ids = ();
	my $readgendata=0; # read genetic data

	open (OUT, "> $out.tmp") 
		or die "\nCan't write to $out.tmp\n\n";
	open (IN, "$in") 
		or die "\nCan't read bcf file with $in\n\n";
		my $nsnps=0;
		while (<IN>){
			chomp;
			## get individual ids
			if (m/^#CHROM/){
				my @aux=split(m/\s+/, $_);
				foreach my $i (9..$#aux){
					$aux[$i]=~ s/(\.sorted)?(\.[s|b]am)?//g;
					push (@ids, $aux[$i]);
				}
				# print OUT join (" ", @ids)."\n";
				$readgendata=1;
				next;
			}
			## read genetic data lines, write gl
			elsif ($readgendata && (!m/[AGCT],[AGCT]/)){ # second condition discards multiallelic snps
				print "SNPs read from input file: $nsnps\r" if ($nsnps%100==0);
				$|++; # flush buffer
				# last if ($nsnps>1000);

				my @aux=split(/\s+/,$_);
				my $id = "$aux[0]".":"."$aux[1]";
				my $line = "$id ";
				@aux = split(m/\s+/, $_);
				next if ($aux[4] =~ m/\,/); # discard multiallelic snps

				my @afs=();
				while ($aux[7]=~ m/AF\=([0-1]\.?[0-9]*e?\-?[0-9]*)\;/g){
					push (@afs, $1);
				}
				# Store last estimated allele frequency
				push (@af_bcf, $afs[$#afs]);

				# get column index for genotype likelihoods from FORMAT field	
				my @aux2=split(/\:/,$aux[8]);
				my $plcol=0;
				while ($aux2[$plcol] ne "PL"){
					$plcol++;
				}
				foreach my $a (@aux[9..$#aux]){
					my @aux2=split(/\:/,$a);
					my $genotype=$aux2[0];
					my $genolhl=$aux2[$plcol];
					# print "$id - $a - genolhl $genolhl\n";
					if ($genolhl =~ /[0-9]+,[0-9]+,[0-9]+/){
						$genolhl =~ s/,/ /g;
						$line.=" $genolhl";
					}
					elsif ($genotype eq '0/0' && $genolhl eq '0'){ # invariant, this used to happen with bcftools 0.1x, should not happen with bcftools 1.x
						$line.=" 0 999 999";	
					}
					elsif ($genotype eq './.' && $genolhl eq '.'){ # missing data
						$line.=" 0 0 0";
					}
					else{
						die ("\nERROR: genotype ($genotype) and genotype likelihoods ($genolhl) format not recognised\n\n");
					}
				}
				print OUT "$line\n";
				$nsnps++;
			}	
		}
	close (IN);

	close (OUT);
	open (OUT, "> $out") 
		or die "\nCan't write to $out\n\n";
	open (IN, "$out.tmp")
		or die "\nCan't read $out.tmp\n\n";
		print OUT scalar(@ids)." $nsnps\n";
		while (<IN>){
			print OUT $_;
		}
	close (IN);
	close (OUT);
	unlink glob "$out.tmp";

	print "SNPs read from input file: $nsnps\r";
	print "\n\n";
	
	return (\@af_bcf, \@ids);
}
# ==============================================================================

# Show copyright
# ==============================================================================
sub author{
    print "\n";
	print "#########################################\n";
	print "  ".basename($0)." version $version     \n";
    print "  (c) Victor Soria-Carrasco             \n";
    print "  victor.soria.carrasco\@gmail.com      \n";
	print "#########################################\n";
}
# ==============================================================================

# Show usage
# ==============================================================================
sub usage{
	my $txt='';
	$txt=shift if (scalar(@_)>0);
	print "\n";
    print "  Usage:\n";
    print "    ".basename($0)."\n";
    print "      -i    <SNPs input file (bcf/vcf format)>\n";
	print "      -p    <population pair file>\n";
	print "      -w    <genome windows file (optional)>\n";
    print "      -o    <output file>\n";
    print "      -ow   <output file for genome windows (optional)>\n";
    print "      -n    <number of compression threads (optional, default=1)>\n";
    print "      -ms   <minimum sampling coverage> (optional, default=0.5)\n";
    print "      -msp  <minimum sampling coverage per population> (optional, default=0.5)\n";
    print "      -mf   <minimum minor allele frequency> (optional, default=0.05)\n";
	print "      -pr   <0|1> prior for calculating empirical posterior probabilities of genotypes\n";
	print "            (0->uniform | 1->Hardy-Weinberg*; optional, default=0)\n";
	print "      -am   <0|1> allele frequency estimation method\n";
	print "            (0->using empirical Bayesian genotypes | 1->using EM algorithm*, default=0)\n";
	print "      -fm   <0|1|2> Fst estimation method\n";
	print "            (0->Hudson | 1->Hudson loci | 2->Weir & Cockerham; optional, default=0)\n";
	print "      -sep  <0|1> estimate allele frequency of each population separately\n";
	print "            (0->no | 1->yes; optional, default=0)\n";
	print "      -pi   <0|1> estimate pi (nucleotide frequency)\n";
	print "            (0->no | 1->yes; optional, default=0)\n";
	print "      -dxy  <0|1> estimate dxy\n";
	print "            (0->no | 1->yes; optional, default=0)\n";
	print "      -da   <0|1> estimate da\n";
	print "            (0->no | 1->yes; optional, default=0)\n";
	print "      -tajd <0|1> estimate Tajima's D\n";
	print "            (0->no | 1->yes; optional, default=0)\n";
	print "      -sl   <integer> (sequence length in bp)\n";
	print "      -k    <0|1> keep temporary files\n";
	print "            (0->no | 1->yes; optional, default=0)\n";
    print "\n";
    print "      * Using a H-W prior will enable maximum-likelihood estimation of allele frequencies\n";
	print "        using the EM algorithm described in Li 2011 doi:10.1093/bioinformatics/btr509\n";
    print "\n";
	print "  Example:\n";
	print "      ".basename($0)." -i input.bcf -p poppairs.csv -o fst_loci.csv\n";
    print "\n\n";
	print $txt;
    exit;
}
# ==============================================================================
