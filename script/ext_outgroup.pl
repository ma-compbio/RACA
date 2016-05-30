#!/usr/bin/perl 

use strict;
use warnings;

my $data_dir = shift;

my $src_f = "$data_dir/Conserved.Segments";
my $out_f = "$data_dir/outgroup.txt";

open(F,"$src_f");
my @lines = <F>;
close(F);
chomp(@lines);

my %hs_outgroup = ();
for (my $i = 0; $i <= $#lines; $i++) {
	my $line = $lines[$i];
	if ($line !~ /^>/) { next; }

	my $si = $i+3;
	while(1) {
		my $sline = $lines[$si];
		if (length($sline) == 0) { last; }
		
		$sline =~ /^(\S+)\.\S+:/;
		my $spc = $1;
		$hs_outgroup{$spc} = 1;
		$si++;
	}
}

open(O,">$out_f");
foreach my $spc (sort keys %hs_outgroup) {
	print O "$spc\n";
}
close(O);
