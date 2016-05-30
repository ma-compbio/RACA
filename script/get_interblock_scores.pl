#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use List::Util qw (min max);

my $insertsize_limit = shift;
my $insert_f = shift;
my $blist_f = shift;
my $size_f = shift;
my $scfprefix = shift;
my $data_dir = shift;
my $out_f = shift;

my $sdtimes = 2;

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

# read synteny blocks
my %hs_blocks = ();
open(F,"$blist_f");
while(<F>) {
	chomp;
	my ($scfnum, $start, $end) = split(/\s+/);
	$hs_blocks{$scfnum}{$start} = $end;
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

my %used = ();
my %mappings = ();
my %sm_mappings = ();
my $cnt_dup = 0;
my $cnt_total = 0;
my $cnt_passed = 0;
my $cnt_diffscf_disterr = 0;

# process paired-ends mappings
my @files = <$data_dir/*>;
my $numfiles = scalar(@files);
my $cntfiles = 0;
foreach my $file (@files) {
	$cntfiles++;

	my $bname = basename($file);
	my @artmp = split(/\./, $bname);
    my @ar = split(/_/, $artmp[0]);

    my $insertname = $ar[-1];
    my $insertsize = $insertsizes_exp{$insertname};
    my $insertsize_cal = $insertsizes_cal{$insertname};
    my $insertsize_negsd = $insertsizes_negsd{$insertname};
    my $insertsize_possd = $insertsizes_possd{$insertname};

	my $maxsize = $insertsize_cal + ($sdtimes * $insertsize_possd);

	my $pline = "";
    open(F,"$file");
    my $linecnt = 0;
    while(<F>) {
        $linecnt++;
        chomp;
        my ($readstr, $readlen, $dir, $scffull, $pos) = split(/\s+/);
		my $readid = substr($readstr, -1, 1);
		$scffull =~ /$scfprefix(\S+)/;
		my $scf = $1;
		my $scfsize = $scf_size{$scf};
		if ($readid != 1 && $readid != 2) { die; }
        
		if ($readid == 1) {
            $pline = $_;
        } else {
			if (length($pline) == 0) { next; }

			$cnt_total++;
            if (defined($used{"$pline $_"})) { $cnt_dup++; next; }

            my ($preadid, $preadlen, $pdir, $pscffull, $ppos) = split(/\s+/,$pline);
			$pscffull =~ /$scfprefix(\S+)/;
			my $pscf = $1;
			my $pscfsize = $scf_size{$pscf};

			# check the existency of scaffolds in the list
			if (!defined($hs_blocks{$pscf} || !defined($hs_blocks{$scf}))) { next; } 	

			my ($scf1len, $scf2len) = ($scf_size{$pscf}, $scf_size{$scf});
			my ($pestart1, $peend1) = ($ppos, $ppos+$preadlen-1);
			my ($pestart2, $peend2) = ($pos, $pos+$readlen-1);
	
			if ($pscf == $scf) { next; }
				
			# search for the first read	
			my ($block_start1, $block_end1) = read_search($pestart1, $peend1, $pscf, \%hs_blocks);
			if ($block_start1 == -1) { next; }	
	
			# search for the second read	
			my ($block_start2, $block_end2) = read_search($pestart2, $peend2, $scf, \%hs_blocks);
			if ($block_start2 == -1) { next; }

			# check distance between reads
			my $dist = 0;		

			if ($insertsize < $insertsize_limit) {
				if ($pdir eq "+") { $dist += ($pscfsize - $ppos + 1); }
				else { $dist += ($ppos + $preadlen - 1); }
					
				if ($dir eq "-") { $dist += ($pos + $readlen - 1); }
				else { $dist += ($scfsize - $pos + 1); }
			} else {
				if ($pdir eq "-") { $dist += ($pscfsize - $ppos + 1); }
				else { $dist += ($ppos + $preadlen - 1); }
					
				if ($dir eq "+") { $dist += ($pos + $readlen - 1); }
				else { $dist += ($scfsize - $pos + 1); }
			}
				
			if ($dist > $maxsize) {
				$cnt_diffscf_disterr++;	
				next;
			}

			$cnt_passed++;

			# check block directions
			my $key1 = "$pscf:$block_start1-$block_end1 ";
			my $key2 = "$scf:$block_start2-$block_end2 ";
			my $keydir1 = "+";
			my $keydir2 = "+";
			
			if ($insertsize < $insertsize_limit) {	# org: +/-
				if ($pdir eq "-") { $keydir1 = "-"; }
				if ($dir eq "+") { $keydir2 = "-"; }
			} else {	# org: -/+
				if ($pdir eq "+") { $keydir1 = "-"; }
				if ($dir eq "-") { $keydir2 = "-"; }
			}	
			$key1 .= $keydir1;
			$key2 .= $keydir2;

			# adjust keys
			if ($pscf > $scf) {
				if ($keydir1 eq "+") { $keydir1 = "-"; }
				else { $keydir1 = "+"; }
					
				if ($keydir2 eq "+") { $keydir2 = "-"; }
				else { $keydir2 = "+"; }
			
				$key1 = "$scf:$block_start2-$block_end2 $keydir2";
				$key2 = "$pscf:$block_start1-$block_end1 $keydir1";
			}

			# store
			if (defined($mappings{$key1}{$key2})) {
				$mappings{$key1}{$key2}++;
			} else {
				$mappings{$key1}{$key2} = 1;
			}
			
			$used{"$pline $_"} = 1;
			$pline = "";
		}
    }
	close(F);
}

# print results
open(O,">$out_f");
foreach my $key1 (sort my_sort keys %mappings) {
	my $rhs = $mappings{$key1};
	foreach my $key2 (sort my_sort keys %$rhs) {
		my $cnt = $$rhs{$key2};
		print O "$key1\t$key2\t$cnt\n";
	}
}

##################################################################
sub read_search {
	my $pestart = shift;
	my $peend = shift;
	my $scfnum = shift;
	my $rhs_blocks = shift;
			
	my $rhs = $$rhs_blocks{$scfnum};
	my ($block_start, $block_end) = (-1,-1);
	foreach my $start (sort {$a<=>$b} keys %$rhs) {
		my $end = $$rhs{$start};
		if ($pestart >= $start && $peend <= $end) {
			$block_start = $start;
			$block_end = $end;
			last;
		}
	}

	return ($block_start, $block_end);
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
