#!/usr/bin/perl

use strict;
use warnings;

my $rec_spc = shift;
my $ref_spc = shift;
my $tar_spc = shift;
my $rec_f = shift;
my $ref_f = shift;
my $tar_f = shift;

open(F,"$ref_f");
my @lines = <F>;
close(F);

my %hs_recchrs_new = ();
my %hs_refchrs_new = ();
my %hs_refchrpos = ();
for (my $i = 0; $i <= $#lines; $i++) {
    my $line = $lines[$i];
    if ($line !~ /^>/) { next; }

    $lines[$i+1] =~ /^$rec_spc\.(\S+):/;
    my $rec_chr = $1;
    $lines[$i+2] =~ /^$ref_spc\.(\S+):(\S+)\-/;
    my $ref_chr = $1;
	my $ref_start = $2;

    $hs_recchrs_new{$rec_chr}{$ref_chr} = 1;
    $hs_refchrs_new{$ref_chr}{$rec_chr} = 1;

	my $pstart = $hs_refchrpos{$rec_chr}{$ref_chr};
	if (defined($pstart)) {
		if ($ref_start < $pstart) {
			$hs_refchrpos{$rec_chr}{$ref_chr} = $ref_start;
		}
	} else {
		$hs_refchrpos{$rec_chr}{$ref_chr} = $ref_start;
	}
}

my %hs_recchrpos = ();
foreach my $rec (sort {$a<=>$b} keys %hs_refchrpos) {
	my $rhs = $hs_refchrpos{$rec};
	foreach my $ref (sort my_sort keys %$rhs) {
		my $start = $$rhs{$ref};
		$hs_recchrpos{$ref}{$start} = $rec;
	}
}

my %hs_refused_num = ();
foreach my $ref_chr (sort my_sort keys %hs_refchrs_new) {
    my $rhs = $hs_refchrs_new{$ref_chr};
    my @keys = keys %$rhs;
    my $keynum = scalar(@keys);
    $hs_refused_num{$ref_chr} = $keynum;
}

my %hs_newrecchrs = ();
my %hs_rec_used = ();

my %used_recchr = ();
foreach my $ref (sort my_sort keys %hs_recchrpos) {
	my $rhs_out = $hs_recchrpos{$ref};
	foreach my $start_out (sort {$a<=>$b} keys %$rhs_out) {
		my $rec_chr = $$rhs_out{$start_out};

		if (defined($used_recchr{$rec_chr})) { 
			next; 
		}

		$used_recchr{$rec_chr} = 1;
		my $rhs = $hs_recchrs_new{$rec_chr};
		my $chrstr = "";
		my $cnt = 0;
		foreach my $ref_chr (sort my_sort keys %$rhs) {
			my $chrnum = substr($ref_chr, 3);
			my $totalusednum = $hs_refused_num{$ref_chr};
			if ($cnt > 0) { $chrstr .= "_"; }

			my $usednum = $hs_rec_used{$chrnum};
			if (defined($usednum)) {
				my $cnum = 97 + $usednum;
				if ($cnum > 122) {
					my $diff = $cnum - 122;
					$chrstr .= $chrnum.chr(122).$diff;
				} else {
					$chrstr .= $chrnum.chr(97+$usednum);
				}
			} else {
				$chrstr .= $chrnum;
				if ($totalusednum > 1) {
					$chrstr .= chr(97);
				}
			}

			if (defined($hs_rec_used{$chrnum})) {
				$hs_rec_used{$chrnum} = $hs_rec_used{$chrnum} + 1;
			} else {
				$hs_rec_used{$chrnum} = 1;
			}

			$cnt++;
		}

		$hs_newrecchrs{$rec_chr} = $chrstr;
	}
}

open(F,"$rec_f");
@lines = <F>;
close(F);
chomp(@lines);

my %hs_newlines = ();
foreach my $line (@lines) {
	my @ar = split(/\s+/, $line);
	my ($rec_chr, $start, $end) = ($ar[0], $ar[1], $ar[2]);
	my $newrec_chr = $hs_newrecchrs{$rec_chr};
	my $newline = "$newrec_chr\t$start\t$end";
	for (my $i = 3; $i <= $#ar; $i++) {
		$newline .= "\t$ar[$i]";
	}
	$hs_newlines{$newrec_chr}{$start} = $newline;	
}

`mv $rec_f $rec_f.org`;
open(O,">$rec_f");
foreach my $rec_chr (sort my_sort2 keys %hs_newlines) {
	my $rhs_starts = $hs_newlines{$rec_chr};
	foreach my $start (sort {$a<=>$b} keys %$rhs_starts) {
		my $outline = $$rhs_starts{$start};
		print O "$outline\n";
	}
}
close(O);

open(F,"$ref_f");
@lines = <F>;
close(F);
chomp(@lines);

%hs_newlines = ();
for (my $i = 0; $i <= $#lines; $i++) {
	my $line = $lines[$i];
	if ($line !~ /^>/) { next; }

	$lines[$i+1] =~ /$rec_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($rec_chr, $rec_start, $rec_end, $rec_dir) = ($1,$2,$3,$4);
	$lines[$i+2] =~ /$ref_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($ref_chr, $ref_start, $ref_end, $ref_dir) = ($1,$2,$3,$4);

	my $newrec_chr = $hs_newrecchrs{$rec_chr};
	my $newline = "$rec_spc.$newrec_chr:$rec_start-$rec_end $rec_dir\n";
	$newline .= "$lines[$i+2]\n"; 
	$hs_newlines{$newrec_chr}{$rec_start} = $newline;	
}

`mv $ref_f $ref_f.org`;
open(O,">$ref_f");
my $newsegid = 1;
foreach my $rec_chr (sort my_sort2 keys %hs_newlines) {
	my $rhs_starts = $hs_newlines{$rec_chr};
	foreach my $start (sort {$a<=>$b} keys %$rhs_starts) {
		my $outline = $$rhs_starts{$start};
		print O ">$newsegid\n";	
		print O "$outline\n";
		$newsegid++;
	}
}
close(O);

open(F,"$tar_f");
@lines = <F>;
close(F);
chomp(@lines);

%hs_newlines = ();
for (my $i = 0; $i <= $#lines; $i++) {
	my $line = $lines[$i];
	if ($line !~ /^>/) { next; }

	$lines[$i+1] =~ /$rec_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($rec_chr, $rec_start, $rec_end, $rec_dir) = ($1,$2,$3,$4);
	$lines[$i+2] =~ /$tar_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($tar_chr, $tar_start, $tar_end, $tar_dir) = ($1,$2,$3,$4);

	my $newrec_chr = $hs_newrecchrs{$rec_chr};
	my $newline = "$rec_spc.$newrec_chr:$rec_start-$rec_end $rec_dir\n";
	$newline .= "$lines[$i+2]\n"; 
	$hs_newlines{$newrec_chr}{$rec_start} = $newline;	
}

`mv $tar_f $tar_f.org`;
open(O,">$tar_f");
$newsegid = 1;
foreach my $rec_chr (sort my_sort2 keys %hs_newlines) {
	my $rhs_starts = $hs_newlines{$rec_chr};
	foreach my $start (sort {$a<=>$b} keys %$rhs_starts) {
		my $outline = $$rhs_starts{$start};
		print O ">$newsegid\n";
		print O "$outline\n";
		$newsegid++;
	}
}
close(O);

sub rec_sort {
	my $va = $a;
	my $vb = $b;
	
	use vars qw(%hs_refchrpos);
	my $rhs = \%hs_refchrpos;

print STDERR "rec_sort: $va $vb\n";

	my $rhs_arefs = $$rhs{$va};
	my $rhs_brefs = $$rhs{$vb};

	my $flag = 0;
	my $order = "aa";
	foreach my $aref (sort my_sort keys %$rhs_arefs) {
		my $astart = $$rhs_arefs{$aref};

		my $bstart = $$rhs_brefs{$aref};
		if (defined($bstart)) {
print STDERR "==> $aref $astart $bstart\n";
			$flag = 1;
			if ($astart < $bstart) { $order = "ab"; } 
			elsif ($astart > $bstart) { $order = "ba"; }
		}
	}
if (($va == 19 && $vb == 28) || ($va == 28 && $vb == 19)) {
print STDERR "==>$flag\t$order\n";

}
	if ($flag == 0) {
		return -1 if ($va < $vb);
		return 1 if ($va > $vb);
		return 0;
	} else {
		if ($order eq "ab") { return -1; }
		elsif ($order eq "ba") { return 1; }
		else { return 0; }
	}
}

sub my_sort {
    $a =~ /chr(\S+)/;
    my $chr1 = $1;
    $b =~ /chr(\S+)/;
    my $chr2 = $1;

	return 0 if ($chr1 eq $chr2);
	return 1 if ($chr1 eq "X" && $chr2 ne "X");
	return -1 if ($chr1 ne "X" && $chr2 eq "X");
 
    return -1 if ($chr1 < $chr2);
    return 1 if ($chr1 > $chr2);
    return 0;
}

sub my_sort2 {
    my $ca = substr($a, 0, 1);
    my $cb = substr($b, 0, 1);
    if ($ca =~ /[Xx]/ && $cb !~ /[Xx]/) {
        return 1;
    } elsif ($ca !~ /[Xx]/ && $cb =~ /[Xx]/) {
        return -1;
    } elsif ($ca =~ /[Xx]/ && $cb =~ /[Xx]/) {
        my $ca2 = substr($a, 1, 1);
        my $cb2 = substr($b, 1, 1);
        if ($ca2 lt $cb2) { return -1; }
        if ($ca2 gt $cb2) { return 1; }

        my $ca3 = substr($a, 2);
        my $cb3 = substr($b, 2);
        if (length($ca3) == 0) { return -1; }
        if (length($cb3) == 0) { return 1; }
        if ($ca3 < $cb3) { return -1; }
        if ($ca3 > $cb3) { return 1; }
        return 0;
    }

    $a =~ /^(\d+)/;
    my $chr1 = $1;
    $b =~ /^(\d+)/;
    my $chr2 = $1;

    return -1 if ($chr1 < $chr2);
    return 1 if ($chr1 > $chr2);

    my $chr1n = substr($a, length($chr1), 1);
    my $chr2n = substr($b, length($chr2), 1);
    if ($chr1n lt $chr2n) { return -1; }
    if ($chr1n gt $chr2n) { return 1; }
    return 0;
}

