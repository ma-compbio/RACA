#!/usr/bin/perl

use strict;
use warnings;

my $str_f = shift;	# rec_chrs.refined.txt
my $genome_f = shift;	# genome.scf.fasta

my $gaps = "N"x100;

# read genome sequences
my $totalbp = 0;
my %hs_seqs = ();
my $spc = "";
my $seq = "";
open(F,"$genome_f");
while(<F>) {
	chomp;
	if ($_ =~ /^>/) {
		if (length($seq) > 0) {
			$hs_seqs{$spc} = $seq;
			$totalbp += length($seq);	
		}

		$seq = "";
		my @ars = split(/\s+/);
		$spc = substr($ars[0],1);
	} else {
		$seq .= $_;
	}
}
close(F);
		
if (length($seq) > 0) {
	$hs_seqs{$spc} = $seq;
	$totalbp += length($seq);	
}

my %hs_used_start = ();
my %hs_used_end = ();

# read chr structure file
my $usedracabp = 0;
my $chrcnt = 1;
my $chrseq = "";
my $pchrid = "";
open(F,"$str_f");
while(<F>) {
	chomp;
	my @ar = split(/\s+/);
	my ($chrid,$start,$end) = ($ar[0],$ar[1],$ar[2]);
	if ($pchrid ne "" && $chrid ne $pchrid) {
		$chrcnt++;
		$chrseq = "";
	}

	if ($ar[3] eq "GAPS") {
		$chrseq .= $gaps;
	} else {
my $csize = $end - $start;
my $ssize = $ar[5] - $ar[4];
if ($csize != $ssize) { print STDERR "DIFF $csize vs $ssize\n"; }

		my ($scfid,$scfstart,$scfend,$scfdir) = ($ar[3],$ar[7],$ar[8],$ar[6]);
		my $scffasta = $hs_seqs{$scfid};
		if (defined($hs_used_start{$scfid})) {
			if ($hs_used_start{$scfid} > $scfstart) {
				$hs_used_start{$scfid} = $scfstart;
			}	
			if ($hs_used_end{$scfid} < $scfend) {
				$hs_used_end{$scfid} = $scfend;
			}	
		} else {
			$hs_used_start{$scfid} = $scfstart;
			$hs_used_end{$scfid} = $scfend;
		} 
		my $scfseq = substr($scffasta,$scfstart,$scfend-$scfstart);
	
		$usedracabp += length($scfseq);
		$chrseq .= $scfseq;	
	}	

	$pchrid = $chrid;
}
close(F);
	
$chrcnt++;
$chrseq = "";

my $usedorgbp = 0;
foreach my $scfid (sort keys %hs_seqs) {
	if (!defined($hs_used_start{$scfid})) {
		next;
	}
	my $scffasta = $hs_seqs{$scfid};
	my $scfstart = $hs_used_start{$scfid};
	my $scfend = $hs_used_end{$scfid};
	my $scfseq = substr($scffasta, $scfstart, $scfend-$scfstart);

	$usedorgbp += length($scfseq);
	$chrcnt++;
}

my $org_coverage = $usedorgbp/$totalbp;
my $raca_coverage = $usedracabp/$totalbp;
print "#ref_length\torg_length\torg_covearge\traca_length\traca_coverage\n";
print "$totalbp\t$usedorgbp\t$org_coverage\t$usedracabp\t$raca_coverage\n";
