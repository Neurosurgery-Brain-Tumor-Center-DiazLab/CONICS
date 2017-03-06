#!/usr/bin/perl -w

use strict;


open(BAM, $ARGV[0]) or die;
my @sample_arr = ();
while(<BAM>){
	my $line = $_; chomp $line;
	if($line =~ /bam$/){
		push @sample_arr, $line;
	}
}
close(BAM);

my $sample_num = scalar@sample_arr;

my $t = 10;
my $t_num = int($sample_num/$t);
my $int = $sample_num % $t;

my @sample_num_arr = ();
for(my$i=0;$i<$t;$i++){
    my $a = $t_num;
    if($int > 0){
        $a++;
        $int--;
    }
    push @sample_num_arr, $a;
    print "sample $i : $a files\n";
}

my %sub_sample = ();
my %all_sample = ();
srand(time|$$); #initialize the random number
for(my$j=0;$j<$t;$j++){
    while(scalar(keys%{$sub_sample{$j}}) < $sample_num_arr[$j]){
        my $k = int(rand $sample_num); #rand_num
        if(!exists$all_sample{$sample_arr[$k]}){
            $all_sample{$sample_arr[$k]} = "x";
            $sub_sample{$j}{$sample_arr[$k]} = "x";
        }
    }
}

my $dir = "$ARGV[1]"."/sample_list/";
for(my$l=0;$l<$t;$l++){


    my $test_output = $dir."Test_".$l;
    my $train_output = $dir."Train_".$l;

    open(TE, "> $test_output") or die;
    open(TR, "> $train_output") or die;

    foreach my $sample (keys %all_sample){
        if(exists $sub_sample{$l}{$sample}){
            print TE "$sample\t\n";
        }
        else{
            print TR "$sample\t\n";
        }

    }

    close(TE);
    close(TR);

}


