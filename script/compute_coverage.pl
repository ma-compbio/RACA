#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use List::Util qw (min max);

my $insertsize_limit = shift;
my $fstart = shift;
my $fend = shift;
my $insert_f = shift;
my $size_f = shift;
my $scfprefix = shift;
my $data_dir = shift;
my $out_dir = shift;

system("mkdir -p $out_dir");

my $sdtimes = 2;

my @files = <$data_dir/*.mapping>;

my %insertsizes_exp = ();
my %insertsizes_cal = ();
my %insertsizes_negsd = ();
my %insertsizes_possd = ();
open(F,"$insert_f");
while(<F>) {
    chomp;
    my ($insert, $size, $calsize, $neg_sd, $pos_sd) = split(/\s+/);
    $insertsizes_exp{$insert} = $size;
    $insertsizes_cal{$insert} = $calsize;
    $insertsizes_negsd{$insert} = $neg_sd;
    $insertsizes_possd{$insert} = $pos_sd;
}
close(F);

my %scf_size = ();
open(F,"$size_f");
while(<F>) {
    chomp;
    my ($scf, $size) = split(/\s+/);
	$scf =~ /$scfprefix(\S+)/;
    $scf_size{$1} = $size;
}
close(F);

for (my $fi = $fstart; $fi <= $fend; $fi++) {
	my $tmpfile = $files[$fi];
	my $fname = basename($tmpfile, ".mapping");
	my $scfnum = substr($fname, length($scfprefix));
	
	my %used = ();
	my %mappings = ();
	my %sm_mappings = ();
	my $cnt_dup = 0;
	my $cnt_total = 0;
	my $cnt_passed = 0;
	my $cnt_samescf_samedir = 0;
	my $cnt_samescf_disterr = 0;
	my $cnt_diffscf_disterr = 0;

	my $cntfiles = 0;

my $file = "$data_dir/$scfprefix$scfnum.mapping";
	$cntfiles++;
	print STDERR "Processing $file\n";

	my $scf = $scfnum;
	my $scfsize = $scf_size{$scf};
	
	# check the existence of the output file
    my $of = "$out_dir/$scfprefix$scf.cov";
	my @ar_finaldata = ();

	my $pline = "";
    open(F,"$file");
    my $linecnt = 0;
    while(<F>) {
        $linecnt++;
        chomp;
        my ($insertname, $readstr, $readlen, $dir, $scffull, $pos) = split(/\s+/);
		my $readid = substr($readstr, -1, 1);
		
		if ($readid != 1 && $readid != 2) { die; }
        
		if ($readid == 1) {
            $pline = $_;
        } else {
			if (length($pline) == 0) { next; }

			$cnt_total++;
            if (defined($used{"$pline $_"})) { $cnt_dup++; next; }

            my ($pinsertname, $preadid, $preadlen, $pdir, $pscffull, $ppos) = split(/\s+/,$pline);

			if ($pinsertname ne $insertname) { die; }
			my $insertsize = $insertsizes_exp{$insertname};
			my $insertsize_cal = $insertsizes_cal{$insertname};
			my $insertsize_negsd = $insertsizes_negsd{$insertname};
			my $insertsize_possd = $insertsizes_possd{$insertname};

			my $maxsize = $insertsize_cal + ($sdtimes * $insertsize_possd);
			my $minsize = $insertsize_cal - ($sdtimes * $insertsize_negsd);

			my ($pestart1, $peend1) = ($ppos, $ppos+$preadlen-1);
			my ($pestart2, $peend2) = ($pos, $pos+$readlen-1);

			# check distance between reads
			if ($pdir eq $dir) { $cnt_samescf_samedir++; next; }

			my $dist = 0;
			if ($insertsize < $insertsize_limit) {
				if ($pdir eq "+" && $dir eq "-") {
					$dist = $pos - $ppos + $readlen; 
				} else {
					$dist = $ppos - $pos + $preadlen;
				}
			} else {
				if ($pdir eq "-" && $dir eq "+") {
					$dist = $pos - $ppos + $readlen; 
				} else {
					$dist = $ppos - $pos + $preadlen;
				}
			}

			if ($dist < $minsize || $dist > $maxsize) {
				$cnt_samescf_disterr++;	
				next;
			}

			my ($prstart, $prend) = ($pestart1, $peend2);
			if ($pestart1 > $pestart2) {
				$prstart = $pestart2;
				$prend = $peend1;
			}

			my $range = "$prstart $prend";
			push(@ar_finaldata, $range);
				
            $used{"$pline $_"} = 1;
			$pline = "";
        }
    }
	close(F);

	my @ar_cov = map { 0 } 1..$scfsize;
	foreach my $range (@ar_finaldata) {
		my ($start, $end) = split(/\s+/, $range);
		for (my $p = $start; $p <= $end; $p++) {
			$ar_cov[$p-1]++;
		}	
	}

	my $outstr = "";	
	for (my $p = 1; $p <= $scfsize; $p++) {
		my $cov = $ar_cov[$p-1];
		$outstr .= "$p\t$cov\n";
	}	
	
	open(O,">$of");
	print O "$outstr";	
	close(O); 
}
