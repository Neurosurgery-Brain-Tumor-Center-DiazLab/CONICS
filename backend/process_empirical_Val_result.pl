#!/usr/bin/perl -w

use strict;

my $base = $ARGV[0];
my $list = $ARGV[1];
my $box_input = $ARGV[2];

open(O, "> $box_input") or die;
open(I, $list) or die "Cannot opne $list\n";
my @file_list = ();
while(<I>){
	my $line = $_; chomp $line;
	if($line =~ /incidenceMatrix.csv$/){
		push @file_list, $line;
	}
}
close(I);

my $dir = "output_".$base."/incidenceMatrix/";
foreach my $file (@file_list){
	my $result_file = $dir.$file;
	open(R, $result_file) or die "$result_file\n";
	my $dummy = <R>;
	my $sum_all = 0;
	my $sum_1 = 0;
	while(<R>){
		my $line = $_; chomp $line;
		my @ele = split /\t/, $line;
		for(my$j=1;$j<scalar@ele;$j++){
			if($ele[$j] == 1){
				$sum_all ++;
				$sum_1 ++;
			}
			elsif($ele[$j] == 0){
				$sum_all ++;
			}
		}
	}
	my $fdr = $sum_1/$sum_all;
	print O "$fdr\n";
	close(R);
}


close(O);
