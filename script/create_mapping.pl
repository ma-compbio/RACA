#!/usr/bin/perl

use strict;
use warnings;

my $refspc = shift;
my $tarspc = shift;
my $src_f = shift;

my %hs_outs = ();
open(F,"$src_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	if ($ar[3] eq "GAPS") { next; }

	my ($chr, $start, $end, $tscf, $tstart, $tend, $tdir, $tblkstart, $tblkend) = split(/\s+/);

	# estimate the coordinates of a block
	my ($newstart, $newend) = (0,0);
	my $rsize = $end - $start;
	my $scfsize = $tend - $tstart;
	
	if ($tblkstart <= $tstart) { $newstart = $start; }
	else {
		my $diff = $tblkstart - $tstart;
		$newstart = $start + $diff;
	}

	if ($tblkend >= $tend) { $newend = $end; } 
	else {
		my $diff = $tend - $tblkend;
		$newend = $end - $diff;
	}

	my $out = "$refspc.$chr:$newstart-$newend +\n";
	$out .= "$tarspc.$tscf:$tblkstart-$tblkend $tdir\n";
	$out .= "\n";
	
	my $key = "$chr:$newstart-$newend";
	$hs_outs{$key} = $out;
}
close(F);

my $segid = 1;
foreach my $key (sort my_sort keys %hs_outs) {
	my $out = $hs_outs{$key};
	print ">$segid\n";
	print "$out";
	$segid++;
}

sub my_sort {
    $a =~ /(\S+):(\S+)\-(\S+)/;
    my ($scf1, $start1, $end1) = ($1,$2,$3);
    $b =~ /(\S+):(\S+)\-(\S+)/;
    my ($scf2, $start2, $end2) = ($1,$2,$3);

    return -1 if ($scf1 < $scf2);
    return 1 if ($scf1 > $scf2);
    return -1 if ($start1 < $start2);
    return 1 if ($start2 < $start1);
    return 0;
}

