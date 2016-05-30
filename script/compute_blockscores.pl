#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);

my $insertsize_thr = shift;
my $insertsize_f = shift;
my $window_size = shift;
my $scfprefix = shift;
my $size_f = shift;
my $reads_dir = shift;
my $cov_dir = shift;
my $data_dir = shift;


`$Bin/get_interblock_scores.pl $insertsize_thr $insertsize_f $data_dir/block_list.txt $size_f $scfprefix $reads_dir $data_dir/inter_block_scores.txt`;

`$Bin/get_intrablock_scores.pl $scfprefix $data_dir/block_list.txt $window_size $cov_dir $data_dir/intra_block_scores.w$window_size.txt`;
