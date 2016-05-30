#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use Cwd;
use Cwd 'abs_path';

# check the number of argument
if ($#ARGV+1 != 1) {
	print STDERR "Usage: ./Run_RACA.pl <parameter file>\n";
	exit(1);
}

my $params_f = $ARGV[0];

# parse parameter file
my %params = ();
open(F,"$params_f");
while(<F>) {
	chomp;
	my $line = trim($_);
	if ($line =~ /^#/ || $line eq "") { next; }
	my ($name, $value) = split(/=/);
	$name = trim($name);
	$value = trim($value);
	if (-f $value || -d $value) {
		$params{$name} = abs_path($value);
	} else {
		$params{$name} = $value;
	}
}
close(F);

check_parameters(\%params);

my $log_dir = $params{"OUTPUTDIR"}."/logs";
my $sepmap_dir = $params{"OUTPUTDIR"}."/mapping_eachscf";
my $cov_dir = $params{"OUTPUTDIR"}."/coverage";
my $sf_dir = $params{"OUTPUTDIR"}."/SFs";

# collect paired-end reads mapping for each target scaffold
print STDERR "Processing paired-end read mappping data...\n";
`$Bin/script/collect_readmapping.pl $params{"READMAPPINGLIB"} $params{"READMAPPINGDIR"} $sepmap_dir`;

# compute coverage for each target scaffold
print STDERR "Computing paired-end read coverage...\n";
`$Bin/script/wrap_compute_coverage.pl $params{"INSERTSIZETHR"} $params{"NCPUS"} $params{"INSERTLIBFILE"} $params{"SCFSIZEFILE"} $params{"SCFPREFIX"} $sepmap_dir $cov_dir $log_dir`; 

`mkdir -p $sf_dir`;
`mkdir -p $params{"OUTPUTDIR"}`;

# make blocks
print STDERR "Constructing syntenic fragments...\n"; 
my $cwd = getcwd();
`mkdir -p $sf_dir`;
`cp $params{"CONFIGSFSFILE"} $sf_dir/config.file`; 
#`cp $params{"MAKESFSFILE"} $sf_dir/Makefile`; 
`sed -e 's:<needtobechanged>:$Bin/code/makeBlocks:' $params{"MAKESFSFILE"} > $sf_dir/Makefile`;
chdir($sf_dir);
`make`;
chdir($cwd);
`$Bin/script/create_blocklist.pl $params{"REFSPC"} $params{"TARSPC"} $sf_dir`; 

# estimate coverage threshold for intra-blocks
print STDERR "Statistical test for blocks within a scaffold...\n";
my $out2 = `$Bin/script/estimate_intrablock_thr.pl $params{"SCFPREFIX"} $sf_dir/block_list.txt $params{"WINDOWSIZE"} $cov_dir $params{"OUTPUTDIR"}/intrablock_thrs.txt`;
my ($intrathr_5p, $intrathr_1p) = split(/\s+/, $out2);
#print STDERR "\tEstimated cutoffs for intra blocks: $intrathr_5p(5%), $intrathr_1p(1%)\n";

# get threshold value
my $intrablock_thr = `$Bin/script/get_thr.pl $params{"MIN_INTRACOV_PERC"} $params{"OUTPUTDIR"}/intrablock_thrs.txt`;
chomp($intrablock_thr);
my $thr_perc = $params{"MIN_INTRACOV_PERC"};
print STDERR "\tEstimated cutoffs for intra blocks: $intrablock_thr ($thr_perc %)\n";

# reconstruct PCFs
`$Bin/script/wrap_recon_pcf.pl $params{"BENADJFILE"} $params{"TREEFILE"} $params{"INSERTSIZETHR"} $params{"RESOLUTION"} $params{"REFSPC"} $params{"TARSPC"} $params{"SCFPREFIX"} $sf_dir $params{"SCFSEQFILE"} $params{"INSERTLIBFILE"} $params{"SCFSIZEFILE"} $intrablock_thr $params{"WINDOWSIZE"} $params{"READMAPPINGDIR"} $cov_dir $params{"OUTPUTDIR"} $params{"IGNORE_ADJS_WO_READS"}`; 

###############################################################
sub check_parameters {
	my $rparams = shift;
	my $flag = 0;
	my $out = "";
	my @parnames = ("INSERTLIBFILE","READMAPPINGDIR","READMAPPINGLIB","NCPUS","SCFSIZEFILE","SCFPREFIX","INSERTSIZETHR","REFSPC","TARSPC","WINDOWSIZE","OUTPUTDIR","RESOLUTION","SCFSEQFILE","TREEFILE","BENADJFILE","CONFIGSFSFILE","MAKESFSFILE","MIN_INTRACOV_PERC","IGNORE_ADJS_WO_READS"); 

	foreach my $pname (@parnames) {
		if (!defined($$rparams{$pname})) {
			$out .= "$pname "; 
			$flag = 1;
		}
	}

	if ($flag == 1) {
		print STDERR "missing parameters: $out\n";
		exit(1);
	}	
}

sub trim {
	my $str = shift;
	$str =~ s/^\s+//;
	$str =~ s/\s+$//;
	return $str;
}
