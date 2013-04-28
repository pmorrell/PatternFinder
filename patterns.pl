#!/usr/bin/perl

##########################################################################
#																									#	
#  Copyright (C) 2007 originally by Donna Toleno and Peter Morrell			#
#                                                                         	#
#  Modifications (C) 2010 by Jeff Ross-Ibarra (rossibarra@gmail.com)			#
#                                                                         	#
#  This program is free software: you can redistribute it and/or modify   	#
#  it under the terms of the GNU General Public License as published by   	#
#  the Free Software Foundation, either version 3 of the License, or      	#
#  (at your option) any later version.                                     #
#                                             		                      	#
#  This program is distributed in the hope that it will be useful,        	#
#  but WITHOUT ANY WARRANTY; without even the implied warranty of         	#
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the				#
#  GNU General Public License <http://www.gnu.org/licenses/>					#
# 	for more details.																			#
##########################################################################

use strict;

my $hap;
my @haplotypes;
my $read=0; 

while ($hap = <stdin>){
	if( $hap=~m/^posit/ ){ $read=1; next; }
	 next unless $read;
    push (@haplotypes, $hap);
}

my @sites = split(//,$haplotypes[0]);#split the first haplotype into an array of sites 
my $N_sites = @sites - 1;                #for the purpose of counting how many sites there are#use chomp instead of subtracting one here.
my $N_haplotypes = @haplotypes; #the length of the array of haplotypes gives the num of haplotypes.
my $a_count=0;
my $b_count = 0; 
my $c_count = 0; 
my $d_count = 0;  
my $pos_A = 0;
my $pos_B = 0;
my $pos_C = 0;
my $pos_D = 0;
my $triplet_num =0;
my $quad_num = 0;
my $first_site = 0;
my $second_site = 0;
my $table_counter1 = 0;
my $table_counter2 = 0;
my @lookup = ();
my @maffer = ();
my $triplet_num=0;
my @first;
my @second;

for ($first_site = 0; $first_site < $N_sites - 1; $first_site++){
	for ($second_site = $first_site + 1; $second_site < $N_sites; $second_site++){
		for (my $i= 0; $i < $N_haplotypes; $i++){
			my @sites = split(//,$haplotypes[$i]);
			push(@first, $sites[$first_site]);
			push(@second, $sites[$second_site]);
		}
		$lookup[$first_site][$second_site] = &four_gamete(\@first,\@second)->[0];
		$maffer[$first_site][$second_site] = &four_gamete(\@first,\@second)->[1];
		@first = ();
		@second = ();
	}
}

# now the loops to find all of the quadruplets and then refer to the table for the four 
# gamete test results and determine if a pattern is found.
for ($pos_A=0; $pos_A < ($N_sites -2); $pos_A++) {
	my $initial_B = $pos_A + 1;
	for ($pos_B= $initial_B; $pos_B < ($N_sites-1); $pos_B++){	
		my $initial_C = $pos_B + 1;
		for ($pos_C= $initial_C; $pos_C < $N_sites; $pos_C++){
			#do triplets to get pattern a too
			$triplet_num++ if $maffer[$pos_A][$pos_C] && $maffer[$pos_B][$pos_C]; 
			$a_count++ if !$lookup[$pos_A][$pos_C] && $lookup[$pos_A][$pos_B] && $lookup[$pos_A][$pos_B] && $lookup[$pos_B][$pos_C];
			
			next if $pos_C == $N_sites;  # do this because otherwise D would be beyond last site
			
			my $initial_D = $pos_C + 1;
			for ($pos_D = $initial_D; $pos_D < $N_sites; $pos_D++){
				$quad_num++ if $maffer[$pos_A][$pos_C] && $maffer[$pos_C][$pos_D];
				$b_count++ if $lookup[$pos_A][$pos_D] && $lookup[$pos_B][$pos_C];
				$c_count++ if $lookup[$pos_A][$pos_D] && $lookup[$pos_B][$pos_C] && !$lookup[$pos_A][$pos_B] && !$lookup[$pos_C][$pos_D];
				$d_count++ if $lookup[$pos_A][$pos_B] && $lookup[$pos_C][$pos_D] && !$lookup[$pos_A][$pos_D] && !$lookup[$pos_B][$pos_C];
			}
		}
	}
}

print "triplets\tacount\tquads\tbcount\tcount\tdcount\n";
if( $triplet_num ){ print "$triplet_num\t", $a_count/$triplet_num; }
else{ print "0\t0"; }
if( $quad_num ){ print "\t$quad_num\t", $b_count/$quad_num, "\t", $c_count/$quad_num, "\t", $d_count/$quad_num, "\n"; }
else{ print "\t0\t0\t0\t0\n" }

sub four_gamete{
	my ($site1_ref, $site2_ref);
	($site1_ref,$site2_ref) = @_;
	my $sum_two = 0;
	my $sum_zero = 0;
	my $both_one = 0;
	my $arrangement1_sum1 = 0;
	my $arrangement2_sum1 = 0;
	my $four_gamete = 0;
	my @sum = ();
	my @maf=qw(0,0);
	my $min=0;
	for (my $j = 0; $j < $N_haplotypes; $j++){
		$sum[$j] = ($site1_ref->[$j]) + ($site2_ref->[$j]);
		$maf[0]+=$site1_ref->[$j];
		$maf[1]+=$site2_ref->[$j];
		if ($sum[$j] == 2){
			$sum_two++;
		}
		if ($sum[$j] == 1){
			if ($site1_ref->[$j]){
				$arrangement1_sum1++;
			}
			if ($site2_ref->[$j]){ 
				$arrangement2_sum1++;
			}
			if (($arrangement1_sum1) and ($arrangement2_sum1)){
				$both_one = 1;
			}
		}
		else{
			if ($sum[$j] == 0){
				$sum_zero++;
			}
		}
	}
	
	if (($both_one) and ($sum_two) and ($sum_zero)){
		$four_gamete = 1;
	}
	else {
		$four_gamete = 0;
	}
	
	$min++ if $maf[0] > 1 && $maf[1] > 1 && $maf[1] < $N_haplotypes-1 && $maf[0] < $N_haplotypes-1;
		
	my @result=( $four_gamete, $min );
	return(\@result);
}