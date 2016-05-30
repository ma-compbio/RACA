#!/usr/bin/perl

use strict;
use warnings;
use List::Util qw[min max];

my $anc_spc = shift; 
my $resolution = shift;
my $ref_spc = shift;
my $out_spc = shift;
my $src_f = shift;
my $tar_f = shift;

open(F,"$src_f");
my @lines = <F>;
close(F);
chomp(@lines);

my %refpos = ();
my %tardirs = ();
for (my $i = 0; $i <= $#lines; $i++) {
	if ($lines[$i] !~ /^>/) { next; }

	my $refline = $lines[$i+1];
	my $tarline = $lines[$i+2];
	$tarline =~ /$ref_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($tchr, $tstart, $tend, $tdir) = ($1, $2, $3, $4);
	my $key = "$ref_spc.$tchr:$tstart-$tend"; 
	$refpos{$key} = $refline;
	$tardirs{$key} = $tdir;
}

open(F,"$tar_f");
@lines = <F>;
close(F);
chomp(@lines);

# pre-process conserved segments
## remove ambiguous outgroup lines
my @lines_new = ();
for (my $i = 0; $i <= $#lines; $i++) {
	if ($lines[$i] !~ /^>/) { next; }

	push(@lines_new, $lines[$i]);
	push(@lines_new, $lines[$i+1]);	# ref_line
	
	my $si = $i+2;
	my $sline = $lines[$si];
	my @artmp = ();
	while (length($sline) > 0) {
		if ($sline =~ /^$out_spc/) {
			push(@artmp, $sline);
		}
		$si++;
		$sline = $lines[$si];	
	}
	
	my $asize = scalar(@artmp);
	if ($asize <= 2) {
		push(@lines_new, @artmp);
	} else {
		my %used = ();
		my %used_index = ();

		for (my $tmpi = 0; $tmpi <= $#artmp; $tmpi++) {	
			my $tmpline = $artmp[$tmpi];
			$tmpline =~ /$out_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
			my ($chr,$start,$end,$dir) = ($1,$2,$3,$4);
			
			if ($tmpi == 0 || $tmpi == $#artmp) { 
				push(@lines_new, $tmpline); 
				$used{$chr}{$start} = $end;
				$used_index{$chr}{$start} = $#lines_new;
				next; 
			} 

			# check block size
			if ($end-$start <= $resolution) { next; }
			
			# check overlap with previous blocks
			my $rhs = $used{$chr};
			if (defined($rhs)) {
				foreach my $tmpstart (sort {$a<=>$b} keys %$rhs) {
					my $tmpend = $$rhs{$tmpstart};

					if ($start < $tmpstart && $end > $tmpstart && $end < $tmpend) {
						my $newend = $tmpstart;
						my $newline = "$out_spc.$chr:$start-$newend $dir";
						push(@lines_new, $newline);
						$used{$chr}{$start} = $newend;
						$used_index{$chr}{$start} = $#lines_new;
					} elsif ($end > $tmpend && $start > $tmpstart && $start < $tmpend) {
						my $newstart = $tmpend;
						my $newline = "$out_spc.$chr:$newstart-$end $dir";
						push(@lines_new, $newline);
						$used{$chr}{$newstart} = $end;
						$used_index{$chr}{$start} = $#lines_new;
					} elsif ($end < $tmpend && $start > $tmpstart) {
						my $len1 = $tmpend - $tmpstart;
						my $len2 = $end - $start;	
						next;	
					} elsif ($end > $tmpend && $start < $tmpstart) {
						my $pindex = $used_index{$chr}{$tmpstart};
						delete $used{$chr}{$tmpstart};
						%used_index = ();
						#delete entry
						my @tmp_lines_new = ();
						for (my $ttti = 0; $ttti <= $#lines_new; $ttti++) {
							if ($ttti == $pindex) { next; }
							my $tttmp_line = $lines_new[$ttti];
							push(@tmp_lines_new, $tttmp_line);
					
							$tttmp_line =~ /$out_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
							my ($ttchr,$ttstart,$ttend,$ttdir) = ($1,$2,$3,$4);
							$used_index{$ttchr}{$ttstart} = $ttti;
						}
						@lines_new = @tmp_lines_new;		
						push(@lines_new, $tmpline);
						$used{$chr}{$start} = $end;
						$used_index{$chr}{$start} = $#lines_new;
					} else {
						push(@lines_new, $tmpline);
						$used{$chr}{$start} = $end;
						$used_index{$chr}{$start} = $#lines_new;
					}
				}
			} else {
				push(@lines_new, $tmpline);
				$used{$chr}{$start} = $end;
				$used_index{$chr}{$start} = $#lines_new;
			}
		}
	}

	push(@lines_new, "");
}

@lines = @lines_new;
my $segid = 1;
my %outlines = ();
for (my $i = 0; $i <= $#lines; $i++) {
	if ($lines[$i] !~ /^>/) { next; }
	my $refsid = substr($lines[$i],1);
	
	my $tarline = $lines[$i+1];
	my $objline = $lines[$i+2];
	if (length($objline) == 0) { next; }
	
	$tarline =~ /$ref_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($tchr, $tstart, $tend, $tdir) = ($1, $2, $3, $4);
	my $key = "$ref_spc.$tchr:$tstart-$tend"; 

	my $refline = $refpos{$key};
	my $dir = $tardirs{$key};
	if (!defined($refline)) {
		next;
	}

	$refline =~ /$anc_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
	my ($rchr, $rstart, $rtend, $rtdir) = ($1, $2, $3, $4);
	my $total_rbsize = $rtend - $rstart;

	# For multiple outgroup blocks
	my $si = $i+2;
	my $sline = $lines[$si];
	while (length($sline) > 0) {
		$si++;
		$sline = $lines[$si];	
	}

	my $num_outblocks = $si - ($i+2);
	if ($num_outblocks == 1) {
		# With one outgroup block

		my $outstr = ">$segid\n";
		$outstr .= "$refline\n";
		if ($dir eq $tdir) {
			$outstr .= "$objline\n";	
		} else {
			$objline =~ /$out_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
			my ($ochr, $ostart, $oend, $odir) = ($1, $2, $3, $4);
			my $finaldir = "";
			if ($odir eq "+") { $finaldir = "-"; }
			else { $finaldir = "+"; }
			$outstr .= "$out_spc.$ochr:$ostart-$oend $finaldir\n";
		}	
		$outstr .= "\n";
		$outlines{$rchr}{$rstart} = $outstr;
		$segid++;
	} else {
		# With more than one outgroup blocks
		my %oblockpos = ();
		my %oblock_size = ();
		my %oblock_start = ();
		my ($pochr, $postart, $poend, $podir) = ("", -1, -1, "");
		my $obstart = 0;
		if ($dir eq $tdir) {
			for (my $ti = $i+2; $ti < $si; $ti++) {
				my $sobjline = $lines[$ti];
				$sobjline =~ /$out_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
				my ($ochr, $ostart, $oend, $odir) = ($1, $2, $3, $4);
				$oblockpos{$ochr}{$ostart} = $oend;
				$oblock_size{$ochr}{$ostart} = $oend - $ostart;

				if ($ochr ne $pochr) { 
					$oblock_start{$ochr}{$ostart} = $obstart;
					$obstart += ($oend - $ostart);
				} elsif ($ostart < $poend && $oend > $postart) {
					# overlap
					# ignore inclusion case
					if ($ostart <= $postart && $poend <= $oend) {
						$oblock_start{$ochr}{$postart} = -1;
						$obstart -= ($poend - $postart);
						$oblock_start{$ochr}{$ostart} = $obstart;
						$obstart += ($oend - $ostart);
					} elsif ($postart <= $ostart && $oend <= $poend) {
						$oblock_start{$ochr}{$ostart} = -1;
					} else {
						my $overlapsize = min($oend,$poend) - max($ostart,$postart);
						$obstart -= $overlapsize;
						$oblock_start{$ochr}{$ostart} = $obstart;
						$obstart += ($oend - $ostart);
					}
				} else {
					$oblock_start{$ochr}{$ostart} = $obstart;
					$obstart += ($oend - $ostart);
				}
				($pochr, $postart, $poend, $podir) = ($ochr, $ostart, $oend, $odir);
			}
		} else {
			for (my $ti = $si-1; $ti >= $i+2; $ti--) {
				my $sobjline = $lines[$ti];
				$sobjline =~ /$out_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
				my ($ochr, $ostart, $oend, $odir) = ($1, $2, $3, $4);
				$oblockpos{$ochr}{$ostart} = $oend;
				$oblock_size{$ochr}{$ostart} = $oend - $ostart;

				if ($ochr ne $pochr) { 
					$oblock_start{$ochr}{$ostart} = $obstart;
					$obstart += ($oend - $ostart);
				} elsif ($ostart < $poend && $oend > $postart) {
					# overlap
					# ignore inclusion case
					if ($ostart <= $postart && $poend <= $oend) {
						$oblock_start{$ochr}{$postart} = -1;
						$obstart -= ($poend - $postart);
						$oblock_start{$ochr}{$ostart} = $obstart;
						$obstart += ($oend - $ostart);
					} elsif ($postart <= $ostart && $oend <= $poend) {
						$oblock_start{$ochr}{$ostart} = -1;
					} else {
						my $overlapsize = min($oend,$poend) - max($ostart,$postart);
						$obstart -= $overlapsize;
						$oblock_start{$ochr}{$ostart} = $obstart;
						$obstart += ($oend - $ostart);
					}
				} else {
					$oblock_start{$ochr}{$ostart} = $obstart;
					$obstart += ($oend - $ostart);
				}
				($pochr, $postart, $poend) = ($ochr, $ostart, $oend);
			}
		}

		# Compute total length
		my $total_obsize = 0;
		foreach my $ochr (keys %oblockpos) {
			my $hs = $oblockpos{$ochr};
			my ($postart, $poend) = (-1,-1);
			foreach my $ostart (sort {$a<=>$b} keys %$hs) {
				my $oend = $$hs{$ostart};
				if ($ostart >= $poend) {
					$total_obsize += ($oend - $ostart);
				} elsif ($ostart < $poend && $oend > $poend) {
					# overlap
					$total_obsize += ($oend - $poend);	
				} elsif ($ostart < $poend && $oend <= $poend) { 
					# inclusion: do nothing	
				} else {
					die;
				}	
				($postart, $poend) = ($ostart, $oend);	
			}
		}

		my $new_start = $rstart;
		$poend = -1;
		$postart = -1;
		$pochr = "chr";
		if ($dir eq $tdir) {
			for (my $ti = $i+2; $ti < $si; $ti++) {
				my $sobjline = $lines[$ti];
				$sobjline =~ /$out_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
				my ($ochr, $ostart, $oend, $odir) = ($1, $2, $3, $4);
				my $obsize = $oblock_size{$ochr}{$ostart};

				my $obstart_rel = $oblock_start{$ochr}{$ostart};
				if ($obstart_rel == -1) { next; }

				my $osfrac = $obstart_rel/$total_obsize;
				my $new_start = $rstart + int($total_rbsize*$osfrac + 0.5); 
				my $ofrac = $obsize/$total_obsize;
				my $rbsize = int($total_rbsize*$ofrac + 0.5);
				my $new_end = $new_start + $rbsize;

				if ($new_end > $rtend) { $new_end = $rtend; }
				my $new_refline = "$anc_spc.$rchr:$new_start-$new_end $rtdir";
				if ($new_start > $new_end) {
					$new_refline = "$anc_spc.$rchr:$new_end-$new_start $rtdir";
				}
	
				my $outstr = ">$segid\n";
				$outstr .= "$new_refline\n";
				if ($dir eq $tdir) {
					$outstr .= "$sobjline\n";	
				} else {
					my $finaldir = "";
					if ($odir eq "+") { $finaldir = "-"; }
					else { $finaldir = "+"; }
					$outstr .= "$out_spc.$ochr:$ostart-$oend $finaldir\n";
				}	
				$outstr .= "\n";
				$outlines{$rchr}{$new_start} = $outstr;
				$segid++;
			}
		} else {
			for (my $ti = $si-1; $ti >= $i+2; $ti--) {
				my $sobjline = $lines[$ti];
				$sobjline =~ /$out_spc\.(\S+):(\S+)\-(\S+) (\S+)/;
				my ($ochr, $ostart, $oend, $odir) = ($1, $2, $3, $4);
				my $obsize = $oblock_size{$ochr}{$ostart};

				my $obstart_rel = $oblock_start{$ochr}{$ostart};
				if ($obstart_rel == -1) { next; }

				my $osfrac = $obstart_rel/$total_obsize;
				my $new_start = $rstart + int($total_rbsize*$osfrac + 0.5); 
				my $ofrac = $obsize/$total_obsize;
				my $rbsize = int($total_rbsize*$ofrac + 0.5);
				my $new_end = $new_start + $rbsize;

				if ($new_end > $rtend) { $new_end = $rtend; }
				my $new_refline = "$anc_spc.$rchr:$new_start-$new_end $rtdir";
				if ($new_start > $new_end) {
					$new_refline = "$anc_spc.$rchr:$new_end-$new_start $rtdir";
				}
	
				my $outstr = ">$segid\n";
				$outstr .= "$new_refline\n";
				if ($dir eq $tdir) {
					$outstr .= "$sobjline\n";	
				} else {
					my $finaldir = "";
					if ($odir eq "+") { $finaldir = "-"; }
					else { $finaldir = "+"; }
					$outstr .= "$out_spc.$ochr:$ostart-$oend $finaldir\n";
				}	
				$outstr .= "\n";
				$outlines{$rchr}{$new_start} = $outstr;
				$segid++;
			}		
		}
	}
}

foreach my $rchr (sort my_sort keys %outlines) {
	my $hs = $outlines{$rchr};
	foreach my $rstart (sort {$a<=>$b} keys %$hs) {
		my $outline = $$hs{$rstart};
		print "$outline";
	}
}

sub my_sort {
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

