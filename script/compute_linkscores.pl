#!/usr/bin/perl

use strict;
use warnings;

my $minavg = shift;
my $blist_f = shift;
my $inter_f = shift;
my $intra_f = shift;

#$minavg = 0.5;
print STDERR "MINAVG=$minavg\n";

my $gMaxScore = 1.0;
my $gMinScore = 0.0000001;
my $gDiff = $gMaxScore - $gMinScore;

my %hs_bid = ();
open(F,"$blist_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	my ($scfnum, $start, $end, $dir, $bid) = ($ar[0], $ar[1], $ar[2], $ar[3], $ar[4]);
	$hs_bid{"$scfnum:$start-$end"} = $bid;
}
close(F);

my %hs_scores = ();

my $sum = 0;
my $cnt = 0;
my %hs_inters = ();
open(F,"$inter_f");
while(<F>) {
	chomp;
	my ($block1, $dir1, $block2, $dir2, $score) = split(/\s+/);
	$hs_inters{"$block1 $dir1"}{"$block2 $dir2"} = $score;
	$sum += $score;
	$cnt++;
}
close(F);

my ($maxinter, $mininter) = (-1,-1);
my $avginter = $sum/$cnt;
foreach my $block1 (keys %hs_inters) {
	my $rhs = $hs_inters{$block1};
	foreach my $block2 (keys %$rhs) {
		my $score = $$rhs{$block2};
		my $pct = $score/$avginter;
		$hs_inters{$block1}{$block2} = $pct;
		if ($maxinter == -1 || $pct > $maxinter) {
			$maxinter = $pct;
		}
		if ($mininter == -1 || $pct < $mininter) {
			$mininter = $pct;
		}
	}
}

my $diff = $maxinter - $mininter;
foreach my $block1 (keys %hs_inters) {
	my $rhs = $hs_inters{$block1};
	foreach my $block2 (keys %$rhs) {
		my $score = $$rhs{$block2};
		my $newscore = $gMinScore + ($score - $mininter)/$diff*$gDiff;
		$hs_scores{$block1}{$block2} = $newscore;	
	}
}

$sum = 0;
$cnt = 0;
my %hs_intras = ();
open(F,"$intra_f");
while(<F>) {
	chomp;
	my ($block1, $dir1, $block2, $dir2, $score) = split(/\s+/);
	#if ($score < $minavg) { next; }

	$hs_intras{"$block1 $dir1"}{"$block2 $dir2"} = $score;
	$sum += $score;
	$cnt++;
}
close(F);

if ($cnt > 0) {
my ($maxintra, $minintra) = (-1,-1);
my $avgintra = $sum/$cnt;
foreach my $block1 (keys %hs_intras) {
	my $rhs = $hs_intras{$block1};
	foreach my $block2 (keys %$rhs) {
		my $score = $$rhs{$block2};
		my $pct = $score/$avgintra;
		$hs_intras{$block1}{$block2} = $pct;
		if ($maxintra == -1 || $pct > $maxintra) {
			$maxintra = $pct;
		}
		if ($minintra == -1 || $pct < $minintra) {
			$minintra = $pct;
		}
	}
}

$diff = $maxintra - $minintra;
foreach my $block1 (keys %hs_intras) {
	my $rhs = $hs_intras{$block1};
	foreach my $block2 (keys %$rhs) {
		my $score = $$rhs{$block2};
		if ($score < $minavg) { next; }
		my $newscore = $gMinScore + ($score - $minintra)/$diff*$gDiff;
		$hs_scores{$block1}{$block2} = $newscore;	
	}
}
}

# print
foreach my $block1 (keys %hs_scores) {
	my $rhs = $hs_scores{$block1};
	foreach my $block2 (keys %$rhs) {
		my $score = $$rhs{$block2};
		print "$block1\t$block2\t$score\n";	
	}
} 
