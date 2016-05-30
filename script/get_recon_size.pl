#!/usr/bin/perl

use strict;
use warnings;

my $scfprefix = shift;
my $src_f = shift;

my $total_size = 0;
my $pchr = "";
my $psize = "";
my %pscfcnt = ();
open(F,"$src_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	my $chr = $ar[0];
	my $size = $ar[2];

	if ($pchr ne "" && $pchr ne $chr) {
		my $pcnt = scalar(keys %pscfcnt);
		print "$pchr\t$psize\t$pcnt\n";
		$total_size += $psize;
		%pscfcnt = ();
	}
	
	if ($ar[3] =~ /^$scfprefix/) { $pscfcnt{$ar[3]} = 1; }

	$pchr = $chr;
	$psize = $size;
}
close(F);
my $pcnt = scalar(keys %pscfcnt);
print "$pchr\t$psize\t$pcnt\n";
$total_size += $psize;
print "Total\t$total_size\n";
