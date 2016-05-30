#!/usr/bin/perl

use strict;
use warnings;

my $perc = shift;
my $data_f = shift;

my %hs_thrs = ();
open(F,"$data_f");
while(<F>) {
	chomp;
	my ($perc, $thr) = split(/\s+/);
	$hs_thrs{$perc} = $thr;
}
close(F);

my $thr = $hs_thrs{int($perc)};
print "$thr\n";
