#!/usr/bin/perl

use strict;
use warnings;

my $blist_f = shift;
my $prob_f = shift;
my $out_f = shift;

my %hs_blocks = ();
my %hs_blocks_dir = ();
open(F, "$blist_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	my ($scf, $start, $end, $dir, $bid) = ($ar[0], $ar[1], $ar[2], $ar[3], $ar[4]);
	$hs_blocks{$bid} = "$scf:$start-$end";
	$hs_blocks_dir{$bid} = $dir;
}
close(F);

my %hs_probs = ();
open(F, "$prob_f");
while(<F>) {
	chomp;
	if ($_ =~ /^#/ || length($_) == 0) { next; }
	my ($bid1, $bid2, $prob) = split(/\s+/);
	
	my $pval = $hs_probs{$bid1}{$bid2};
	if (defined($pval)) {
		# check consistency
		die if ($prob != $pval);
		next;
	}

	$pval = $hs_probs{-$bid2}{-$bid1};	
	if (defined($pval)) {
		# check consistency
		die if ($prob != $pval);
		next;
	}

	$hs_probs{$bid1}{$bid2} = $prob;
}
close(F);

my %hs_used = ();
open(O, ">$out_f");
foreach my $bid1 (keys %hs_probs) {
	my $rhs = $hs_probs{$bid1};
	foreach my $bid2 (keys %$rhs) {
		my $prob = $$rhs{$bid2};

		my $bstr1 = $hs_blocks{abs($bid1)};
		my $bdir1 = $hs_blocks_dir{abs($bid1)};
		if ($bid1 < 0) {
			if ($bdir1 eq "+") { $bdir1 = "-"; }
			else { $bdir1 = "+"; }
		}	
	
		my $bstr2 = $hs_blocks{abs($bid2)};
		my $bdir2 = $hs_blocks_dir{abs($bid2)};
		if ($bid2 < 0) {
			if ($bdir2 eq "+") { $bdir2 = "-"; }
			else { $bdir2 = "+"; }
		}	

		if (defined($hs_used{"$bstr1 $bdir1"}{"$bstr2 $bdir2"})) { next; }
		my $rbdir1 = "+";
		if ($bdir1 eq "+") { $rbdir1 = "-"; }
		my $rbdir2 = "+";
		if ($bdir2 eq "+") { $rbdir2 = "-"; }
		if (defined($hs_used{"$bstr2 $rbdir2"}{"$bstr1 $rbdir1"})) { next; }
	
		print O "$bstr1 $bdir1\t$bstr2 $bdir2\t$prob\n";		

		$hs_used{"$bstr1 $bdir1"}{"$bstr2 $bdir2"} = $prob;	
	}
}
close(O);
