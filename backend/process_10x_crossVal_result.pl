#!/usr/bin/perl -w

use strict;

my $t = 10;
my $base = $ARGV[0];
my $path_to_rscript = $ARGV[1];
my $fdr_t = $ARGV[2];
my $box_input = $ARGV[3];

open(O, "> $box_input") or die;

for(my$i=0;$i<$t;$i++){
	my $test = $base."/cpmMatrix/Test_".$i."_cpmMatrix.csv";
	my $train = $base."/cpmMatrix/Train_".$i."_cpmMatrix.csv";
	my $result = $base."/incidenceMatrix/Result_".$i;

	my $str = "$path_to_rscript  ./backend/plot_bean_CPM.R $test $train $result $fdr_t";
	system "$str\n";
	my $result_file = $result.".incidenceMatrix.csv";
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
