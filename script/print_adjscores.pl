#!/usr/bin/perl

use strict;
use warnings;

my $scf_prefix = shift;
my $chr_f = shift;
my $score_f = shift;

my %hs_scores = ();
open(F,"$score_f");
while(<F>) {
	chomp;
	if ($_ =~ /^>/ || length($_) == 0) { next; }

	my ($block1, $dir1, $block2, $dir2, $weight) = split(/\s+/);
	$hs_scores{"$scf_prefix$block1"}{"$scf_prefix$block2"} = $weight;
	$hs_scores{"$scf_prefix$block2"}{"$scf_prefix$block1"} = $weight;
}
close(F);

my @lines = ();
open(F,"$chr_f");
while(<F>) {
	chomp;
	if ($_ =~/GAPS/) { next; }
	push(@lines, $_);
}
close(F);
chomp(@lines);

for (my $i = 0; $i < $#lines; $i++) {
	my $line1 = $lines[$i];
	my $line2 = $lines[$i+1];
	my ($chr1, $start1, $end1, $scf1, $tscfstart1, $tscfend1, $dir1, $tstart1, $tend1) = split(/\s+/, $line1);
	my ($chr2, $start2, $end2, $scf2, $tscfstart2, $tscfend2, $dir2, $tstart2, $tend2) = split(/\s+/, $line2);

	if ($chr1 ne $chr2) { next; }
	my $weight = $hs_scores{"$scf1:$tstart1-$tend1"}{"$scf2:$tstart2-$tend2"};
	print "$chr1\t$end1\t$start2\t$weight\n";
}
