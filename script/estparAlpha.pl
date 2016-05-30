#!/usr/bin/perl

use strict;
use warnings;

my $reladj_f = shift;
my $cons_f = shift;
my $link_f = shift;

my %hs_cons = ();
open(F,"$cons_f");
while(<F>) {
	chomp;
	my ($bid1, $dir1, $bid2, $dir2, $score) = split(/\s+/);
	$hs_cons{"$bid1 $dir1"}{"$bid2 $dir2"} = $score;
}
close(F);

my %hs_link = ();
open(F,"$link_f");
while(<F>) {
	chomp;
	my ($bid1, $dir1, $bid2, $dir2, $score) = split(/\s+/);
	$hs_link{"$bid1 $dir1"}{"$bid2 $dir2"} = $score;
}
close(F);

my $sse_cons = 0;
my $sse_link = 0;

open(F,"$reladj_f");
while(<F>) {
	chomp;
	my ($bid1, $dir1, $bid2, $dir2) = split(/\s+/);
	my $key1 = "$bid1 $dir1";
	my $key2 = "$bid2 $dir2";

	my $rdir1 = "+";
	if ($dir1 eq "+") { $rdir1 = "-"; }	
	my $rdir2 = "+";
	if ($dir2 eq "+") { $rdir2 = "-"; }	
	my $rkey1 = "$bid2 $rdir2";
	my $rkey2 = "$bid1 $rdir1";

	my $cons_score = $hs_cons{$key1}{$key2};
	if (!defined($cons_score)) {
		$cons_score = $hs_cons{$rkey1}{$rkey2};
	}
	if (!defined($cons_score)) {
		$sse_cons += 1;
	} else {
		$sse_cons += (1 - $cons_score)**2;
	}	
	
	my $link_score = $hs_link{$key1}{$key2};
	if (!defined($link_score)) {
		$link_score = $hs_link{$rkey1}{$rkey2};
	}
	if (!defined($link_score)) {
		$sse_link += 1;
	} else {
		$sse_link += (1 - $link_score)**2;
	}	
}
close(F);

my $alpha = $sse_link/($sse_cons + $sse_link);
print "$alpha\n";
