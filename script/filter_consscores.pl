#!/usr/bin/perl

use strict;
use warnings;

my $ref = shift;
my $blist_f = shift;
my $score_f = shift;
my $src_dir = shift;
my $outgroup_f = shift;

my %hs_blocks = ();
my %hs_dirs = ();
open(F,"$blist_f");
while(<F>) {
	chomp;
	my ($scf, $start, $end, $dir, $bid) = split(/\s+/);
	$hs_blocks{"$scf:$start-$end"} = $bid;	
	$hs_dirs{"$scf:$start-$end"} = $dir;
}
close(F);

my %hs_refjoins = ();
open(F,"$src_dir/$ref.joins");
while(<F>) {
	chomp;
	if ($_ =~ /^#/) { next; }
	my @ar = split(/\s+/);
	my ($bid1, $bid2) = ($ar[0], $ar[1]);
	if ($bid1 eq "") {
		($bid1, $bid2) = ($ar[1], $ar[2]);
	}

	$hs_refjoins{"$bid1:$bid2"} = 1;
	my $rbid1 = -1 * $bid1;	
	my $rbid2 = -1 * $bid2;	
	$hs_refjoins{"$rbid2:$rbid1"} = 1;
}
close(F);


open(F,"$outgroup_f");
my @outgroups = <F>;
close(F);
chomp(@outgroups);

my %hs_outjoins = ();
foreach my $outgroup (@outgroups) {
	my $outgroup_f = "$src_dir/$outgroup.joins";
	open(F,"$outgroup_f");
	while(<F>) {
		chomp;
		if ($_ =~ /^#/) { next; }
		my @ar = split(/\s+/);
		my ($bid1, $bid2) = ($ar[0], $ar[1]);
		if ($bid1 eq "") {
			($bid1, $bid2) = ($ar[1], $ar[2]);
		}
		
		$hs_outjoins{$outgroup}{"$bid1:$bid2"} = 1;
		my $rbid1 = -1 * $bid1;	
		my $rbid2 = -1 * $bid2;	
		$hs_outjoins{$outgroup}{"$rbid2:$rbid1"} = 1;
	}
	close(F);
}

open(F,"$score_f");
while(<F>) {
	chomp;
	my ($blk1, $dir1, $blk2, $dir2, $score) = split(/\s+/);
	my $bid1 = $hs_blocks{$blk1};
	my $odir1 = $hs_dirs{$blk1};
	if ($dir1 ne $odir1) { $bid1 = -1 * $bid1; }	

	my $bid2 = $hs_blocks{$blk2};
	my $odir2 = $hs_dirs{$blk2};
	if ($dir2 ne $odir2) { $bid2 = -1 * $bid2; }	

	my $flag = 0;

	if (defined($hs_refjoins{"$bid1:$bid2"})) {
		$flag = 1;
	} else {
		foreach my $outgroup (@outgroups) {
			if (defined($hs_outjoins{$outgroup}{"$bid1:$bid2"})) {
				$flag = 1;
				last;
			}		
		}
	}

	if ($flag == 1) {
		print "$_\n";
	} else {
		#print STDERR "Filtered $bid1 $bid2 $_\n";
	}
}
close(F);
