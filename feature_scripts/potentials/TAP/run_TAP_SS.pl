#!/usr/bin/env perl
use warnings;
use strict;
use Math::Trig;

# input
my $f_pdb = $ARGV[0];
my $f_stride = $ARGV[1];
# output
my $f_out = $ARGV[2];

my $dir_work = $ARGV[3];
my $f_u2 = "$dir_work/potentials/TAP/parse_stride_ss_angle.pl";
require "$f_u2";
my $f_refer = "$dir_work/potentials/TAP/TAP_refer.9.txt";

###########################
# global configure values
my @aminos = ("ASP","PRO","LYS","ILE","TRP","CYS","GLY","PHE","GLN","SER","ASN","LEU","VAL","TYR","GLU","ARG","THR","ALA","MET","HIS");
my $interval = 9;
my $RT = 0.582;
my $sigma = 0.02;
my $penalty = 0;

##############################
# reference state information
my (%refers, %Ma);
for my $i (@aminos){ $Ma{$i} = 0; }
my $Ma_sum = 0;
open IN, "$f_refer" or die $!;
while(my $line = <IN>){
  chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;
  my @items = split(/\s+/, $line);
  $Ma_sum += $items[2]; 
  $Ma{$items[0]} += $items[2];
  my $k = "$items[0],$items[1]";
  $refers{$k} = $items[2];
}close IN;

open OUT, ">$f_out";
########################################
# secondary structure and torsion angles
my %stride;
if(! -e $f_stride || -z $f_stride){ die "Error, calculating TAP, no stride file."; }
parse_stride_ss_angle($f_stride, \%stride);

my $N = keys %stride;
if($N == 0){
  close OUT;
  exit;
}
     
my $asp_sum = 0;
for my $k1 (sort {$a<=>$b} keys %stride){
    if(@{$stride{$k1}} != 5){ next; }
    if(exists $stride{$k1-1} && exists $stride{$k1+1}){
        if(@{$stride{$k1-1}} != 5){ next; }
        if(@{$stride{$k1+1}} != 5){ next; }
        # phi1 and psi1
        my $c1 = get_class($stride{$k1-1}->[3], $interval);
        my $c2 = get_class($stride{$k1-1}->[4], $interval);
        # phi2 and psi2 : center
        my $c3 = get_class($stride{$k1}->[3], $interval);
        my $c4 = get_class($stride{$k1}->[4], $interval);
        # phi1 and psi1
        my $c5 = get_class($stride{$k1+1}->[3], $interval);
        my $c6 = get_class($stride{$k1+1}->[4], $interval);
        if($c1 == 0 || $c2 == 0 || $c3 == 0 || $c4 == 0 || $c5 == 0 || $c6 == 0){ next; } 
        my $ca1 = get_five_class($c1, $c3, $c5);
        my $ca2 = get_five_class($c2, $c4, $c6);
        my $ca3 = get_four_class($c1, $c3, $c5);
        my $ca4 = get_four_class($c2, $c4, $c6);       

        my $class3 = "$ca1,$ca2,$ca3,$ca4";
 
        my $cur_residue = $stride{$k1}->[1];
        my $f_observed_a = 0;
        my $f_reference = 0;
        my $ma;
        my $k1_ss;
        if(exists $stride{$k1}){
          $k1_ss = $stride{$k1}->[2];
        }else{
          next;
        }
        my $key_Ma = "$cur_residue";
        if(exists $Ma{$key_Ma}){
          $ma = $Ma{$key_Ma};
        }else{
          next;
        }
        my $rkey = "$cur_residue,$k1_ss,$class3";
        if(exists $refers{$rkey}){
          $f_observed_a = $refers{$rkey} / $ma;
        }
        my $asp = $penalty;
        my $ma_sum = 0;
        for my $tkey (keys %Ma){
          my $tk = "$tkey,$k1_ss,$class3";
          if(exists $refers{$tk}){
            $ma_sum += $refers{$tk};
          }
        }
        $f_reference = $ma_sum / $Ma_sum;
        if($f_observed_a && $f_reference){
          $asp = ($RT*log(1+$ma*$sigma)) - ($RT*log(1+($ma*$sigma*($f_observed_a/$f_reference))));
          $asp_sum += $asp;
        }else{
          $asp_sum += $penalty;
        }  
        print OUT "$k1\t$cur_residue\t$k1_ss\t$asp\n";
    }
}
close OUT;


sub get_class{
  my ($deg, $index) = @_;
  my $n = 0;
  my $bin;
  $bin = 360 / $index;
  for(my $i=-180; $i<180; $i+=$bin){
    my $j = $i+$bin;
    $n++;
    if($i == (180-$bin)){
      if($deg >= $i && $deg <= $j){ return $n;}
    }else{
      if($deg >= $i && $deg < $j){ return $n;}
    }
  }
  return 0;
}

sub get_five_class{
  my ($a1, $a2, $a3) = @_;
  if($a1 == $a2 && $a1 == $a3){ return 1; }
  if($a1 != $a2 && $a1 != $a3 && $a2 != $a3){ return 2; }
  if($a1 == $a2 && $a1 != $a3){ return 3; }
  if($a1 == $a3 && $a1 != $a2){ return 4; }
  if($a2 == $a3 && $a1 != $a2){ return 5; }
}

sub get_four_class{
  my ($a1, $a2, $a3) = @_;
  if($a2 == $a1 && $a2 == $a3){ return 1; }
  if($a2 != $a1 && $a2 != $a3){ return 2; }
  if($a2 == $a1 && $a2 != $a3){ return 3; }
  if($a2 == $a3 && $a2 != $a1){ return 4; }
}




