#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

# input
my $f_pdb = $ARGV[0];
my $f_ss = $ARGV[1];
# output
my $f_out = $ARGV[2];

my $dir_work = $ARGV[3];
my $f_u1 = "$dir_work/utils/ave_min_max.pl";
my $f_u2 = "$dir_work/potentials/RSAP/parse_stride.pl";
require "$f_u1";
require "$f_u2";
my $f_refer = "$dir_work/potentials/RSAP/RSAP_refer.txt";

###########################################
my $penalty = 0;
my $RT = 0.582;
my $sigma = 0.02;
my @aminos = ("ASP","PRO","LYS","ILE","TRP","CYS","GLY","PHE","GLN","SER","ASN","LEU","VAL","TYR","GLU","ARG","THR","ALA","MET","HIS");
my %aminoID;
for(my $i=0; $i<20; $i++){
  $aminoID{$aminos[$i]} = $i;
}

my @refers;
my %Ma;
my $Ma_sum = 0;
my $id = 0;
open IN, "$f_refer" or die $!;
while(my $line = <IN>){
  chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;
  my @items = split(/\s+/, $line);
  push(@refers, \@items);
  my $sum = ave_min_max(\@items, 4);
  $Ma_sum += $sum;
  $Ma{$aminos[$id]} = $sum;
  $id++;
}close IN;

open OUT, ">$f_out";
########################
# SS
my %ss;
if(! -e $f_ss || -z $f_ss){ die "Error, calculating TAP, no stride file."; }
parse_stride($f_ss, \%ss);
my $N = keys %ss;
if($N == 0){
  close OUT;
  exit;
}

for my $k1 (sort {$a<=>$b} keys %ss){
      my $amino1 = $ss{$k1}->[0];
      my $rid1 = $k1;
      my $ss1 = $ss{$k1}->[1];
      my $acc1 = $ss{$k1}->[2];
      my $class = get_class($acc1);

      my $asp = $penalty;
      my $f_o = 0;
      my $f_r = 0;
      my $ma = $Ma{$amino1};
        
      my $add = 0;
      if($ss1 eq 'C'){ $add = 0; }
      if($ss1 eq 'E'){ $add = 1; }
      if($ss1 eq 'H'){ $add = 2; }
      if($ma == 0){ next; }
      $f_o = $refers[$aminoID{$amino1}][$add*10+$class-1] / $ma;

      my $fr1 = 0;
      for(my $i=0; $i<20; $i++){
          $fr1 += $refers[$i][$add*10+$class-1];
      }
      $f_r = $fr1 / $Ma_sum; 
      if($f_r){ 
            $asp = ($RT*log(1+$ma*$sigma)) - ($RT*log(1+($ma*$sigma*($f_o/$f_r))));
      }
      print OUT "$rid1\t$amino1\t$asp\n";
}
close OUT;


sub get_class{
  my $ras = shift;
  if($ras >= 0 && $ras < 0.1){ return 1; }
  if($ras >= 0.1 && $ras < 0.2){ return 2; }
  if($ras >= 0.2 && $ras < 0.3){ return 3; }
  if($ras >= 0.3 && $ras < 0.4){ return 4; }
  if($ras >= 0.4 && $ras < 0.5){ return 5; }
  if($ras >= 0.5 && $ras < 0.6){ return 6; }
  if($ras >= 0.6 && $ras < 0.7){ return 7; }
  if($ras >= 0.7 && $ras < 0.8){ return 8; }
  if($ras >= 0.8 && $ras < 0.9){ return 9; }
  if($ras >= 0.9){ return 10; }
  return 0;
}



