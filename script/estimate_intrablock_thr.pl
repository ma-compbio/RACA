#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use Statistics::Descriptive::Weighted; 

my $scfprefix = shift;
my $block_f = shift;
my $windowsize = shift;
my $cov_dir = shift;
my $out_f = shift; 

my $flanksize = 10000;	
my $stepsize = int($windowsize/2);

my %hs_blocks = ();
open(F,"$block_f");
while(<F>) {
	chomp;
	my ($scfnum, $start, $end) = split(/\s+/);
	$hs_blocks{$scfnum}{$start} = $end;
}
close(F);

my @values = ();
my @weights = ();
foreach my $scfnum (sort {$a<=>$b} keys %hs_blocks) {
	my $rhs = $hs_blocks{$scfnum};
	my @starts = sort {$a<=>$b} keys %$rhs;
	if (scalar(@starts) > 1) { next; }

	# read coverage data
	my $f = "$cov_dir/$scfprefix${scfnum}.cov";
	if (!(-f $f)) {
		print STDERR "File $f doesn't exist...\n"; next; 
	}

	my %hs_cov = ();
	my $scfsize = 0;
	my $gsum = 0;
	my $gcnt = 0;
	open(F,"$f");
	while(<F>) {
		chomp;
		my ($pos, $cov) = split(/\s+/);
		$hs_cov{$pos} = $cov;
		$scfsize = $pos;
		$gsum += $cov;
		$gcnt++;
	}
	close(F);
	my $gavg = $gsum/$gcnt;

	for (my $p = $flanksize+1; $p < $scfsize - $flanksize; $p += $stepsize) {
		my ($wstart, $wend) = ($p, $p+$windowsize-1);
		if ($wend >= $scfsize-$flanksize) { $wend = $scfsize-$flanksize-1; $p = $scfsize-$flanksize; }
		
		my $sum = 0;
		my $cnt = 0;
		for (my $pp = $wstart; $pp <= $wend; $pp++) {
			my $cov = $hs_cov{$pp};
			$sum += $cov;
			$cnt++;
		}
		my $avg = $sum/$cnt;
		my $pct = $avg/$gavg;
		push(@values, $pct);
		push(@weights, 1);
	}
 	
}

my $stat = Statistics::Descriptive::Weighted::Full->new();
$stat->add_data(\@values, \@weights);
my $mean = $stat->mean();
my $median = $stat->median();
my $cnt = $stat->count();
my $min = $stat->min();
my $max = $stat->max();
my $val01 = $stat->percentile(1);
my $val05 = $stat->percentile(5);
print "$val05\t$val01\n";

# Output to a file
open(O,">$out_f");
print O "0\t0\n";
for (my $pthr = 1; $pthr <= 99; $pthr++) {
	my $val = $stat->percentile($pthr);
	print O "$pthr\t$val\n";
}
close(O);
