#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib/perl";
use Parallel::ForkManager;

my $insertsize_thr = shift;
my $num_p = shift;
my $insert_f = shift;
my $size_f = shift;
my $scfprefix = shift;
my $data_dir = shift;
my $out_dir = shift;
my $log_dir = shift;

system("mkdir -p $log_dir");

my @files = <$data_dir/*.mapping>;
my $fcnt = scalar(@files);
my $cbin = int($fcnt/$num_p);

my $pm = new Parallel::ForkManager($num_p);
for (my $pi = 0; $pi < $num_p; $pi++) {
	my $fstart = $cbin*$pi;
	my $fend = $fstart + $cbin - 1;
	if ($pi == $num_p - 1) { $fend = $fcnt; }

	$pm->start and next;

	my $cmd = "$Bin/compute_coverage.pl $insertsize_thr $fstart $fend $insert_f $size_f $scfprefix $data_dir $out_dir > $log_dir/compute_coverage.p$pi.log 2>&1";
	system($cmd);

	$pm->finish; 
}
$pm->wait_all_children;
