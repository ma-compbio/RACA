#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);

my $scf_prefix = shift;
my $resolution = shift;
my $rec_spc = shift;
my $tar_spc = shift;
my $ref_spc = shift;
my $out_f = shift;
my $block_f = shift;
my $src_f = shift;
my $seg_f = shift;
my $fa_f = shift;	# new
my $out_dir = shift;

open(F,"$out_f");
my @outspcs = <F>;
close(F);
chomp(@outspcs);
system("mkdir -p $out_dir");
system("mkdir -p $out_dir/tmp");

my $rec_chrs_f = "$out_dir/rec_chrs.txt";
`$Bin/create_chrs.pl $scf_prefix $block_f $fa_f $src_f > $rec_chrs_f`;

my $tar_seg_f = "$out_dir/rec_chrs.$tar_spc.segments.txt";
`$Bin/create_mapping.pl $rec_spc $tar_spc $rec_chrs_f > $tar_seg_f`;

my $ref_seg_f = "$out_dir/rec_chrs.$ref_spc.segments.txt";
`$Bin/infer_mapping_ref.pl $rec_spc $tar_spc $ref_spc $tar_seg_f $seg_f > $ref_seg_f`; 

# reorder segments based on a reference
`$Bin/reorder.simple.pl $rec_spc $ref_spc $tar_spc $rec_chrs_f $ref_seg_f $tar_seg_f`;

my $ref_seg_refined_f = "$out_dir/rec_chrs.$ref_spc.segments.refined.txt";
`$Bin/merge_pos2ex.pl $rec_spc $ref_spc $ref_seg_f > $ref_seg_refined_f`; 

my $tar_seg_refined_f = "$out_dir/rec_chrs.$tar_spc.segments.refined.txt";
`$Bin/merge_pos2ex.pl $rec_spc $tar_spc $tar_seg_f > $tar_seg_refined_f`; 

foreach my $out_spc (@outspcs) {
	my $out_seg_f = "$out_dir/rec_chrs.$out_spc.segments.txt";
	`$Bin/infer_mapping.pl $rec_spc $resolution $ref_spc $out_spc $ref_seg_f $seg_f > $out_seg_f`; 

	my $out_seg_refined_f = "$out_dir/rec_chrs.$out_spc.segments.refined.txt";
	`$Bin/merge_pos2ex.pl $rec_spc $out_spc $out_seg_f > $out_seg_refined_f`;
}

my $rec_chrs_refined_f = "$out_dir/rec_chrs.refined.txt";
`$Bin/merge_rec.pl $rec_chrs_f > $rec_chrs_refined_f`; 

`$Bin/get_recon_size.pl $scf_prefix $rec_chrs_refined_f > $out_dir/rec_chrs.size.txt`;

`$Bin/print_adjscores.pl $scf_prefix $out_dir/rec_chrs.txt $out_dir/raca.out > $out_dir/rec_chrs.adjscores.txt`; 

`mv $ref_seg_f $out_dir/tmp/`;
`mv $ref_seg_f.org $out_dir/tmp/`;
`mv $tar_seg_f $out_dir/tmp/`;
`mv $tar_seg_f.org $out_dir/tmp/`;
foreach my $out_spc (@outspcs) {
	my $out_seg_f = "$out_dir/rec_chrs.$out_spc.segments.txt";
	`mv $out_seg_f $out_dir/tmp/`;
}
`mv $rec_chrs_f $out_dir/tmp/`;
`mv $rec_chrs_f.org $out_dir/tmp/`;
