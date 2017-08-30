#!/usr/bin/perl -w

use strict;

if(scalar@ARGV != 2){die "Usage: MkBedFromSequenza.pl [_segments.txt file from sequenza][output_file]\n";}

open(I, $ARGV[0]) or die "cannot open segments.txt file: $ARGV[0]\n"; #open segments.txt files from Sequenza
print "reading $ARGV[0] ...\n";
my $header = <I>;
my %segment = (); #key1=chromosome, key2=start.position, key3=end.position, value=lcopy number
my %depth_ratio = (); #key1=chromosome, key2=start.position, key3=end.position, value=depth ratio
my %c_number = ();
#my %c_xy = ();
my $t = 1500000;
while(<I>){
	my $line = $_; chomp $line;
	$line =~ s/\"//g;
	my @ele = split /\t/, $line;
	if($ele[2]-$ele[1] > $t){
#		if($ele[9] != 2){
			$ele[0] =~ s/chr//gi; 
			$ele[0] =~ s/_//g; 
			$ele[0] =~ s/-//g; 
			$segment{$ele[0]}{$ele[1]}{$ele[2]} = $ele[9];
			$depth_ratio{$ele[0]}{$ele[1]}{$ele[2]} = $ele[6];
			if($ele[0] =~ m/[0-9]/){$c_number{$ele[0]} = "x";}
#			else{$c_xy{$ele[0]} = "x";}
#		}
	}
}
close(I);

my @chromosome = (); 
foreach my $chrm_num (sort {$a<=>$b} keys %c_number){
	push @chromosome, $chrm_num;
}
#foreach my $chrm_xy (sort keys %c_xy){
#	push @chromosome, $chrm_xy;
#}

open(O, "> $ARGV[1]") or die "cannot open output file: $ARGV[1]\n";
print "writing results to $ARGV[1]\n";
foreach my $chrm (@chromosome){
	my $tag = 0;
	my $start_prev = 0; my $end_prev = 0;	my $end_dist = 0; my $status = "";  
	my $length = 0; my $depth_sum = 0; my $depth_length = 0; my $depth_avg = 0; 
	
	foreach my $start (sort {$a<=>$b} keys %{$segment{$chrm}}){
		foreach my $end (sort {$a<=>$b} keys %{$segment{$chrm}{$start}}){
			
			
			my $cn = $segment{$chrm}{$start}{$end};
			my $depth = $depth_ratio{$chrm}{$start}{$end};
			$length = $end - $start +1;
			
			if($tag==0){
				if($cn == 2){next;}
				else{
					$start_prev = $start; $end_prev=$end; $end_dist=$end;
					$tag=1;
					if($cn > 2){$status = "A";}
					elsif($cn<2){$status = "D";}
					$depth_sum = ($depth*$length);
					$depth_length = $length;
				}
			}
			
			else{
				if($start-$end_dist > $t){
					$depth_avg = $depth_sum/$depth_length;
					print O "$chrm\t$start_prev\t$end_prev\t$chrm\:$start_prev\:$end_prev\:$depth_avg\:$status\n";
					if($cn==2){$tag=0;}
					else{
						$start_prev = $start; $end_prev=$end; $end_dist=$end;
						if($cn > 2){$status = "A";}
						elsif($cn<2){$status = "D";}
						$depth_sum = ($depth*$length);
						$depth_length = $length;
					}
				}
				else{
					if($cn == 2){$end_dist=$end;}
					elsif($cn > 2){
						if($status eq "A"){
							$end_prev = $end; $end_dist=$end;
							$depth_sum += ($depth*$length);
							$depth_length += $length;
						}
						else{
							$depth_avg = $depth_sum/$depth_length;
							print O "$chrm\t$start_prev\t$end_prev\t$chrm\:$start_prev\:$end_prev\:$depth_avg\:$status\n";
							if($cn==2){$tag=0;}
							else{
								$start_prev = $start; $end_prev=$end; $end_dist=$end;
								if($cn > 2){$status = "A";}
								elsif($cn<2){$status = "D";}
								$depth_sum = ($depth*$length);
								$depth_length = $length;
							}
						}
					}
					elsif($cn < 2){
						if($status eq "A"){
							$depth_avg = $depth_sum/$depth_length;
							print O "$chrm\t$start_prev\t$end_prev\t$chrm\:$start_prev\:$end_prev\:$depth_avg\:$status\n";
							if($cn==2){$tag=0;}
							else{
								$start_prev = $start; $end_prev=$end; $end_dist=$end;
								if($cn > 2){$status = "A";}
								elsif($cn<2){$status = "D";}
								$depth_sum = ($depth*$length);
								$depth_length = $length;
							}
						}
						else{
							$depth_sum += ($depth*$length);
							$depth_length += $length;
							$end_prev = $end; $end_dist=$end;
						}
					}
				}
			}
		}
	}
	if($tag ==1){
		$depth_avg = $depth_sum/$depth_length;
		print O "$chrm\t$start_prev\t$end_prev\t$chrm\:$start_prev\:$end_prev\:$depth_avg\:$status\n";
		
	}
}
close(O);
