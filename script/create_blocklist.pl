#!/usr/bin/perl

use strict;
use warnings;

my $rspc = shift;
my $spc = shift;
my $data_dir = shift;

my %data = ();
my %data_blockid = ();
my %data_dir = ();
my %data_rchr = ();
my $f = "$data_dir/Conserved.Segments";
my $o = "$data_dir/block_list.txt";

my $blockid = -1;
my $rchr = "";
open(F,"$f");
while(<F>) {
	chomp;

	if ($_ =~ /^>/) {
		$blockid = substr($_,1);
	} elsif ($_ =~ /^$rspc\.(\S+):(\S+)\-(\S+) (\S+)/) {
		$rchr = $1;
	} elsif ($_ =~ /^$spc\.(\S+):(\S+)\-(\S+) (\S+)/) {
		my ($scf, $start, $end, $dir) = ($1,$2,$3,$4);
		$data{$scf}{$start} = $end;	
		$data_blockid{$scf}{$start} = $blockid;	
		$data_dir{$scf}{$start} = $dir;	
		$data_rchr{$scf}{$start} = $rchr;	
	}
}
close(F);

open(O,">$o");
foreach my $scf (sort keys %data) {
	my $rhs = $data{$scf};
	$scf =~ /\D+(\d+)/;
	my $scfnum = $1;
	foreach my $start (sort {$a<=>$b} keys %$rhs) {
		my $end = $$rhs{$start};
		my $blockid = $data_blockid{$scf}{$start};
		my $dir = $data_dir{$scf}{$start};
		my $rchr = $data_rchr{$scf}{$start};
		print O "$scfnum\t$start\t$end\t$dir\t$blockid\t$rchr\n";
	} 
}
close(O);
