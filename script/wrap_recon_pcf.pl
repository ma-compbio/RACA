#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use FindBin qw($Bin);

my $benadj_f = shift;
my $tree_f = shift;
my $insertsizelimit = shift;
my $resolution = shift;
my $ref_spc = shift;
my $tar_spc = shift;
my $scf_prefix = shift;
my $src_dir = shift;
my $fa_f = shift;
my $insertsize_f = shift;
my $size_f = shift;
my $intrablock_thr = shift;
my $window_size = shift;
my $reads_dir = shift;
my $cov_dir = shift;
my $out_dir = shift;
my $ignore_adjs_wo_reads = shift;

system("mkdir -p $out_dir");
# create mapping data
print STDERR "Computing syntenic fragment mapping score...\n";
`$Bin/compute_blockscores.pl $insertsizelimit $insertsize_f $window_size $scf_prefix $size_f $reads_dir $cov_dir $src_dir`; 

my $inter_f = "$src_dir/inter_block_scores.txt";
my $intra_f = "$src_dir/intra_block_scores.w$window_size.txt"; 
# extract outgroup species
`$Bin/ext_outgroup.pl $src_dir`;

# estimate breakpoint distance
if (-f "$src_dir/bpdist.txt") { `rm $src_dir/bpdist.txt`; } 
`$Bin/estimate_bpdist_ref.pl $tar_spc $ref_spc $src_dir >> $src_dir/bpdist.txt`;
`$Bin/estimate_bpdist_out.pl $resolution $tar_spc $src_dir >> $src_dir/bpdist.txt`;

# check tree.txt file
if (!(-f "$tree_f")) {
	print STDERR "File doesn't exist: $tree_f\n";
	die;
}

my $tout = `cat $tree_f`;

# eatimate JK model parameter
print STDERR "Estimating the JK model parameter... ";
open(F,"$src_dir/$ref_spc.joins");
my $tmp = <F>;
close(F);
chomp($tmp);
my $numblocks = substr($tmp, 1);
my $jkalpha = `$Bin/estparJK.pl $numblocks $tar_spc $tree_f $src_dir/bpdist.txt`;
chomp($jkalpha);
print STDERR "$jkalpha\n";

# create new genome file
`$Bin/create_gfile.pl $src_dir`;

print STDERR "Computing adjacency probabilities...\n";
# compute adjacency probabilities
my $curdir = getcwd;
chdir($src_dir);
`$Bin/../code/inferAdjProb $jkalpha $tree_f Genomes.Order.new`; 
chdir($curdir);

# refine adjacency probabilities
`$Bin/refine_adjprob.pl $src_dir/block_list.txt $src_dir/adjacencies.prob $src_dir/block_consscores.txt`;

# remove adjacencies without ref and outgroup support
`mv $src_dir/block_consscores.txt $src_dir/block_consscores.txt.org`;
`$Bin/filter_consscores.pl $ref_spc $src_dir/block_list.txt $src_dir/block_consscores.txt.org $src_dir $src_dir/outgroup.txt > $src_dir/block_consscores.txt`;

# compute link scores
`$Bin/compute_linkscores.pl $intrablock_thr $src_dir/block_list.txt $inter_f $intra_f > $src_dir/block_linkscores.txt`;

# estimate the controlling parameter alpha
print STDERR "Estimating the controlling parameter alpha... ";
my $reladj_f = $benadj_f; 
# check tree.txt file
if (!(-f "$reladj_f")) {
	print STDERR "File doesn't exist: $reladj_f\n";
	die;
}
my $alpha = `$Bin/estparAlpha.pl $reladj_f $src_dir/block_consscores.txt $src_dir/block_linkscores.txt`;
chomp($alpha);
print STDERR "$alpha\n";

# run RACA
print STDERR "Building chromosome fragments...\n";
my $cmd = "$Bin/../code/raca $ignore_adjs_wo_reads $alpha $src_dir/block_list.txt $src_dir/block_consscores.txt $src_dir/block_linkscores.txt > $out_dir/raca.out";
`$cmd`;

# create chromosome data
$cmd = "$Bin/create_pcfdata.pl $scf_prefix $resolution RACA $tar_spc $ref_spc $src_dir/outgroup.txt $src_dir/block_list.txt $out_dir/raca.out $src_dir/Conserved.Segments $fa_f $out_dir"; 
`$cmd`;
