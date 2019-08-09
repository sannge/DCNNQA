#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

# input
my $f_pdb = $ARGV[0];
# output
my $f_out = $ARGV[1];

my $dir_work = $ARGV[2];
my $f_refer = "$dir_work/potentials/SSDP/K_refer_CB.8.0_300_5";

################################
# global configuration values
my $aid = 2;
my $range1 = 0;
my $range2 = 300;
my $bin = 5;
my $separation = 8; # distance cutoff
my $penalty = 0;
my $RT = 0.582;
my $sigma = 0.02;

################################
# reference state information
my %refer_ijl;
my $Ma_sum = 0;
open IN, "$f_refer.ijk" or die $!;
while(my $line = <IN>){
  chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;
  my @items = split(/\s+/, $line);
  $refer_ijl{$items[0]} = $items[1];
  if($items[0] =~ /UNK/){ next; }
  $Ma_sum += $items[1];
}close IN;
my %refer_ij;
open IN, "$f_refer.ij" or die $!;
while(my $line = <IN>){
  chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;
  my @items = split(/\s+/, $line);
  $refer_ij{$items[0]} = $items[1];
}close IN;
my %refer_l;
open IN, "$f_refer.k" or die $!;
while(my $line = <IN>){
  chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;
  my @items = split(/\s+/, $line);
  $refer_l{$items[0]} = $items[1];
}close IN;

open OUT, ">$f_out";
my %coords;
get_coords($f_pdb, $aid, \%coords);
my $N = keys %coords;
if($N == 0){
  close OUT;
  exit;
}


my $asp_sum = 0;
for my $k1 (sort {$a<=>$b} keys %coords){
      my $amino1 = $coords{$k1}->[3];
      my $rid1 = $coords{$k1}->[4];
      for my $k2 (sort {$a<=>$b} keys %coords){
        if($k2 <= $k1){ next; }
        my $amino2 = $coords{$k2}->[3];
        my $rid2 = $coords{$k2}->[4];
        my $kij = abs($k1 - $k2); 
        my $dist = sqrt(($coords{$k1}->[0]-$coords{$k2}->[0])**2 + ($coords{$k1}->[1]-$coords{$k2}->[1])**2 + ($coords{$k1}->[2]-$coords{$k2}->[2])**2);
        if($dist <= $separation ){ # exist potential score
          my $asp = $penalty;
          my $class = get_class($range1, $range2, $bin, $kij);
          if($class == 0){ next; }
          my $key_ijl = "$amino1,$amino2,$class";
          my $key_ij = "$amino1,$amino2";
          my $key_l = $class;

          my $f_o = 0;
          my $f_r = 0;
          my $ma = 0;
          if(exists $refer_ij{$key_ij}){ $ma=$refer_ij{$key_ij}; }
          if($ma == 0){ next; }
 
          if(exists $refer_ijl{$key_ijl} && exists $refer_ij{$key_ij}){
            $f_o = $refer_ijl{$key_ijl} / $refer_ij{$key_ij};
          }
          if(exists $refer_l{$key_l}){
            if($refer_l{$key_l}){
              $f_r = $refer_l{$key_l} / $Ma_sum;
            }
          }
          if($f_r){ 
            $asp = ($RT*log(1+$ma*$sigma)) - ($RT*log(1+($ma*$sigma*($f_o/$f_r))));
          }
          $asp_sum += $asp;
          print OUT "$rid1\t$amino1\t$rid2\t$amino2\t$asp\n";
        }
      }
}
close OUT;

sub get_class{
  my ($range_1, $range_2, $bin, $d) = @_;
  my $n = 0;
  for(my $i=$range_1; $i<$range_2; $i+=$bin){
    my $j = $i+$bin;
    $n++;
    if($j >= $range_2){
      if($d>= $i && $d<=$j){ return $n;}
    }else{
      if($d>= $i && $d<$j){ return $n;}
    }
  }
  return 0;
}

sub get_coords{
  my ($f_pdb, $aid, $coords) = @_;
  # aid = 1, CA
  # aid = 2, CB
  my $atomID;
  if($aid == 1){ $atomID = "CA"; }
  if($aid == 2){ $atomID = "CB"; }
  open IN, "$f_pdb" or die $!;
  my $rid = 0;
  while(my $line = <IN>){
    chomp $line;
    if(substr($line, 0, 6) eq 'ATOM  '){
      my $atom_name = substr($line,12,4);
      $atom_name =~ s/^\s+//; $atom_name =~ s/\s+$//;
      if($atom_name eq $atomID){
        my ($id1, $id2, $x, $y, $z);
        $id1 = substr($line, 6, 5);
        $id2 = substr($line, 22, 4);
        $x = substr($line, 30, 8);
        $y = substr($line, 38, 8);
        $z = substr($line, 46, 8);
        my $residue = substr($line, 17, 3);
        $id1 =~ s/^\s+//; $id1 =~ s/\s+$//;
        $id2 =~ s/^\s+//; $id2 =~ s/\s+$//;
        $x =~ s/^\s+//; $x =~ s/\s+$//;
        $y =~ s/^\s+//; $y =~ s/\s+$//;
        $z =~ s/^\s+//; $z =~ s/\s+$//;
        $rid++;
        push(@{$coords->{$rid}}, $x);
        push(@{$coords->{$rid}}, $y);
        push(@{$coords->{$rid}}, $z);
        push(@{$coords->{$rid}}, $residue);
        push(@{$coords->{$rid}}, $id2);
      }
    }
  }
  close IN;
}

