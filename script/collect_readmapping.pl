#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $filelib_f = shift;
my $src_dir = shift;
my $out_dir = shift;

# check output directory
if (-d $out_dir) {
	CHECKDIR:
	print STDERR "\nThe output directory $out_dir already exists. New data will be appended to existing files. Do you want to remove and recreate it (recommended)? (y/n)";
	my $ans = <>;
	chomp($ans);
	if ($ans eq "y") {
		system("rm -rf $out_dir");
	} elsif ($ans eq "n") {
		; 
	} else {
		print STDERR "Unrecognized answer.\n";
		goto CHECKDIR; 
	} 
}

system("mkdir -p $out_dir");

# read insert library names for each mapping file
my %hs_fname2lib = ();
open(F,"$filelib_f");
while(<F>) {
	chomp;
	if ($_ =~ /^#/) { next; }
	my ($fname, $libname) = split(/\s+/);
	$hs_fname2lib{$fname} = $libname;
}
close(F);

my @files = <$src_dir/*>;
my $numfiles = scalar(@files);

my $fcnt = 1;
foreach my $file (@files) {
	my $bname = basename($file);
	my $libname = $hs_fname2lib{$bname};
	
	my %hs_outs = ();
	my $pline = "";
	my $pscf = "";
	open(F,"$file");
	while(<F>) {
		chomp;
		my @ar = split(/\s+/);
		my $readnum = $ar[0];
		if (length($readnum) != 1) { $readnum = substr($readnum, -1, 1); }
		my $scf = $ar[3];

		if ($readnum == 1) {
			$pscf = $scf;
			$pline = $_;
		} elsif ($readnum == 2) {
			if ($scf eq $pscf) { 
				my $out1 = "$libname $pline";
				my $out2 = "$libname $_";
				if (defined($hs_outs{$scf})) {
					my $rar = $hs_outs{$scf};
					push(@$rar, $out1);
					push(@$rar, $out2);
				} else {
					$hs_outs{$scf} = [$out1,$out2];
				}
			}
			
			$pline = "";
			$pscf = "";
		}
	}
	close(F);
	
	foreach my $scf (keys %hs_outs) {
		my $f = "$out_dir/$scf.mapping";
		my $rar = $hs_outs{$scf};
		my $outstr = "";	
		foreach my $line (@$rar) {
			$outstr .= "$line\n";
		}
		open(O,">>$f");
		print O "$outstr";
		close(O);
	}
	$fcnt++;
}
