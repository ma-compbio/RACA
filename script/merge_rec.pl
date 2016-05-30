#!/usr/bin/perl

use strict;
use warnings;

my $src_f = shift;

open(F,"$src_f");
my @lines = <F>;
close(F);
chomp(@lines);

my ($psid, $pstart, $pend, $pscf, $pscfstart, $pscfend, $pdir) = ("",-1,-1,"",-1,-1,-1,-1,""); 
for (my $i = 0; $i <= $#lines; $i++) {
	my $line = $lines[$i];
	my @ar = split(/\s+/, $line);
	my $sid = $ar[0];

	if ($ar[3] eq "GAPS") {
		print "$psid\t$pstart\t$pend\t$pscf\t$pscfstart\t$pscfend\t$pdir\n";
		($psid, $pstart, $pend, $pscf, $pscfstart, $pscfend, $pdir) = ("",-1,-1,"",-1,-1,-1,-1,""); 
		print "$line\n";	# print gap line
	} else {
		if ($psid eq "") {
			($psid, $pstart, $pend, $pscf, $pscfstart, $pscfend, $pdir) = ($ar[0],$ar[1],$ar[2],$ar[3],$ar[4],$ar[5],$ar[6]); 
		} else {
if ($ar[0] ne $psid || $ar[3] ne $pscf || $ar[6] ne $pdir || ($pdir eq "+" && !($ar[4] > $pscfstart && $ar[5] > $pscfend)) || ($pdir eq "-" && !($ar[4] < $pscfstart && $ar[5] < $pscfend))) {
				print "$psid\t$pstart\t$pend\t$pscf\t$pscfstart\t$pscfend\t$pdir\n";
				($psid, $pstart, $pend, $pscf, $pscfstart, $pscfend, $pdir) = ($ar[0],$ar[1],$ar[2],$ar[3],$ar[4],$ar[5],$ar[6]); 
			} else {
				$pend = $ar[2];
				if ($ar[5] > $pscfend) { $pscfend = $ar[5]; }
				else { $pscfstart = $ar[4]; }
			}
		}
	} 
}
		
print "$psid\t$pstart\t$pend\t$pscf\t$pscfstart\t$pscfend\t$pdir\n";
