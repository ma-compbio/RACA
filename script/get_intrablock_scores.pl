#!/usr/bin/perl

use strict;
use warnings;

my $scfprefix = shift;
my $block_f = shift;
my $windowsize = shift;
my $cov_dir = shift;
my $out_f = shift;

my $flanksize = 50000;	
my $stepsize = int($windowsize/2);

my %hs_blocks = ();
my %hs_blockschr = ();
open(F,"$block_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	my ($scfnum, $start, $end, $rchr) = ($ar[0], $ar[1], $ar[2], $ar[5]);
	$hs_blocks{$scfnum}{$start} = $end;
	$hs_blockschr{$scfnum}{$start}{$end} = $rchr;
}
close(F);

my @finals = ();
foreach my $scfnum (sort {$a<=>$b} keys %hs_blocks) {
	my $rhs = $hs_blocks{$scfnum};
	my @starts = sort {$a<=>$b} keys %$rhs;
	if (scalar(@starts) == 1) { next; }

	# read coverage data
	my $f = "$cov_dir/$scfprefix${scfnum}.cov";
	if (!(-f $f)) {
		print STDERR "File $f doesn't exist...\n"; next; 
	}

	my %hs_cov = ();
	open(F,"$f");
	my $gsum = 0;
	my $gcnt = 0;
	while(<F>) {
		chomp;
		my ($pos, $cov) = split(/\s+/);
		$hs_cov{$pos} = $cov;
		$gsum += $cov;
		$gcnt++;
	}
	close(F);
	my $gavg = $gsum/$gcnt;

	for (my $i = 0; $i < $#starts; $i++) {
		my $start1 = $starts[$i];
		my $end1 = $$rhs{$start1};
		my $start2 = $starts[$i+1];
		my $end2 = $$rhs{$start2};

		my $rchr1 = $hs_blockschr{$scfnum}{$start1}{$end1};
		my $rchr2 = $hs_blockschr{$scfnum}{$start2}{$end2};

		my ($brstart, $brend) = ($end1, $start2);
		if ($end1 > $start2) { ($brstart, $brend) = ($start2, $end1); }

		$brstart -= $flanksize;
		$brend += $flanksize;

		my $minavg = -1;
		my ($min_wstart, $min_wend) = (0,0);
		for (my $p = $brstart; $p < $brend; $p += $stepsize) {
			my ($wstart, $wend) = ($p, $p+$windowsize);
			if ($wend > $brend) { $wend = $brend; $p = $brend; }

			my $sum = 0;
			my $cnt = 0;
			for (my $pp = $wstart+1; $pp<=$wend; $pp++) {
				my $cov = $hs_cov{$pp};
				$sum += $cov;
				$cnt++;
			}
			my $avgcov = $sum/$cnt;
			if ($minavg == -1 || $avgcov < $minavg) {
				$minavg = $avgcov;
				($min_wstart, $min_wend) = ($wstart, $wend);
			}	
		}
		my $type = "INTER";
		if ($rchr1 eq $rchr2) { $type = "INTRA"; }
		my $minpct = $minavg/$gavg;
		my $out = "$scfnum:$start1-$end1\t$scfnum:$start2-$end2\t$brstart-$brend\t$min_wstart-$min_wend\tMIN=$minavg\t$gavg\t$minpct\t$type\t$rchr1 $rchr2";
		push(@finals, $out);
	} 	
}

open(O,">$out_f");
foreach my $out (@finals) {
	my @ar = split(/\s+/, $out);
	print O "$ar[0] +\t$ar[1] +\t$ar[6]\n";
}
close(O);
