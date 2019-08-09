#!/usr/bin/env perl
use warnings;
use strict;
use Math::Trig;

# input
my $f_pdb = $ARGV[0];
my $f_stride = $ARGV[1];
my $f_out = $ARGV[2];

my $dir_work = $ARGV[3];
my $f_u2 = "$dir_work/potentials/PAP/parse_stride.pl";
require "$f_u2";
my $f_refer = "$dir_work/potentials/PAP/PAP_SS_refer.30.txt";

##############################
# global configuration values
my $interval = 30;
my $penalty = 0;
my $RT = 0.582;
my $sigma = 0.02;
my @aminos = ("ASP","PRO","LYS","ILE","TRP","CYS","GLY","PHE","GLN","SER","ASN","LEU","VAL","TYR","GLU","ARG","THR","ALA","MET","HIS");

#####################################
# extract reference state information
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
#########################################
# get the secondary structure information
my %ss;
if(! -e $f_stride || -z $f_stride){ die "Error, when calculating PAP. no stride file."; }
parse_stride($f_stride, \%ss);

my %coords;
my $asp_sum = 0;
get_coords($f_pdb, \%coords);
my $N = keys %coords;
if($N == 0){
  close OUT;
  exit;
}

for my $k1 (sort {$a<=>$b} keys %coords){
    if(@{$coords{$k1}} != 3){ next; }
    if(exists $coords{$k1-1} && exists $coords{$k1+1}){
        if(@{$coords{$k1-1}} != 3){ next; }
        if(@{$coords{$k1+1}} != 3){ next; }
        # N
        my ($xn1, $yn1, $zn1) = ($coords{$k1-1}->[0][0], $coords{$k1-1}->[0][1], $coords{$k1-1}->[0][2]);
        my ($xn2, $yn2, $zn2) = ($coords{$k1}->[0][0], $coords{$k1}->[0][1], $coords{$k1}->[0][2]);
        my ($xn3, $yn3, $zn3) = ($coords{$k1+1}->[0][0], $coords{$k1+1}->[0][1], $coords{$k1+1}->[0][2]);
        my ($xn12, $yn12, $zn12)  = ($xn2-$xn1, $yn2-$yn1, $zn2-$zn1);
        my ($xn23, $yn23, $zn23)  = ($xn3-$xn2, $yn3-$yn2, $zn3-$zn2);
        my $cosN = ($xn12*$xn23 + $yn12*$yn23 + $zn12*$zn23) / (sqrt(($xn12**2)+($yn12**2)+($zn12**2))*sqrt(($xn23**2)+($yn23**2)+($zn23**2)));
        my $degreeN = acosd($cosN);
        my $classN = get_class($interval, $degreeN);
        if($classN == 0){ next; }
        # CA
        my ($xa1, $ya1, $za1) = ($coords{$k1-1}->[1][0], $coords{$k1-1}->[1][1], $coords{$k1-1}->[1][2]);
        my ($xa2, $ya2, $za2) = ($coords{$k1}->[1][0], $coords{$k1}->[1][1], $coords{$k1}->[1][2]);
        my ($xa3, $ya3, $za3) = ($coords{$k1+1}->[1][0], $coords{$k1+1}->[1][1], $coords{$k1+1}->[1][2]);
        my ($xa12, $ya12, $za12)  = ($xa2-$xa1, $ya2-$ya1, $za2-$za1);
        my ($xa23, $ya23, $za23)  = ($xa3-$xa2, $ya3-$ya2, $za3-$za2);
        my $cosA = ($xa12*$xa23 + $ya12*$ya23 + $za12*$za23) / (sqrt(($xa12**2)+($ya12**2)+($za12**2))*sqrt(($xa23**2)+($ya23**2)+($za23**2)));
        my $degreeA = acosd($cosA);
        my $classA = get_class($interval, $degreeA);
        if($classA == 0){ next; }
        # C
        my ($xc1, $yc1, $zc1) = ($coords{$k1-1}->[2][0], $coords{$k1-1}->[2][1], $coords{$k1-1}->[2][2]);
        my ($xc2, $yc2, $zc2) = ($coords{$k1}->[2][0], $coords{$k1}->[2][1], $coords{$k1}->[2][2]);
        my ($xc3, $yc3, $zc3) = ($coords{$k1+1}->[2][0], $coords{$k1+1}->[2][1], $coords{$k1+1}->[2][2]);
        my ($xc12, $yc12, $zc12)  = ($xc2-$xc1, $yc2-$yc1, $zc2-$zc1);
        my ($xc23, $yc23, $zc23)  = ($xc3-$xc2, $yc3-$yc2, $zc3-$zc2);
        my $cosC = ($xc12*$xc23 + $yc12*$yc23 + $zc12*$zc23) / (sqrt(($xc12**2)+($yc12**2)+($zc12**2))*sqrt(($xc23**2)+($yc23**2)+($zc23**2)));
        my $degreeC = acosd($cosC);
        my $classC = get_class($interval, $degreeC);
        if($classC == 0){ next; }

        my $class3 = "$classN,$classA,$classC"; 
        #my $tripeptide = "$amino_t1,$amino_t2,$amino_t3";
        my $cur_residue = $coords{$k1}->[0][3];
        my $f_observed_a = 0;
        my $f_reference = 0;
        my $ma;
        my $k1_2 = $coords{$k1}->[0][4];
        my $k1_ss;
        if(exists $ss{$k1_2}){
          $k1_ss = $ss{$k1_2};
        }else{
          next;
        }
        my $key_Ma = "$cur_residue";
        if(exists $Ma{$key_Ma}){
          $ma = $Ma{$key_Ma};
        }else{
          next;
        }
        #print "$target\t$model\t$k1\t$cur_residue..$ma\n"; 
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

sub acosd { acos($_[0])*180/pi }
sub get_class{
  my ($bin, $degree) = @_;
  my $class = 0;
  my $n = 0;
  for(my $i=0; $i<180; $i+=$bin){
    my $j = $i + $bin;
    $n++;
    if($j >= 180){
      if($degree >= $i && $degree <= $j){ $class = $n; }
    }else{
      if($degree >= $i && $degree < $j){ $class = $n; }
    }
  }
  return $class;
}

sub get_coords{
  my ($f_pdb, $coords) = @_;
  open IN, "$f_pdb" or die $!;
  my $cid = 0;
  my $cur_id2 = 0;
  while(my $line = <IN>){
    chomp $line;
    if(substr($line, 0, 6) eq 'ATOM  '){
      my ($id1, $id2, $x, $y, $z);
      $id1 = substr($line, 6, 5);
      $id2 = substr($line, 22, 4);
      $x = substr($line, 30, 8); $y = substr($line, 38, 8); $z = substr($line, 46, 8);
      my $residue = substr($line, 17, 3);
      $id1 =~ s/^\s+//; $id1 =~ s/\s+$//;
      $id2 =~ s/^\s+//; $id2 =~ s/\s+$//;
      $x =~ s/^\s+//; $x =~ s/\s+$//; $y =~ s/^\s+//; $y =~ s/\s+$//; $z =~ s/^\s+//; $z =~ s/\s+$//;
      my $atom_name = substr($line,12,4);
      $atom_name =~ s/^\s+//; $atom_name =~ s/\s+$//;

      if($atom_name eq "N" || $atom_name eq "CA" || $atom_name eq "C"){
        if($residue ~~ @aminos){
          if($id2 != $cur_id2){
            $cur_id2 = $id2;
            $cid++;
          }
          my @tmp;
          push(@tmp, $x); push(@tmp, $y); push(@tmp, $z);
          push(@tmp, $residue);
          push(@tmp, $id2);
          push(@{$coords->{$cid}}, \@tmp);
        }
      }
    }
  }
  close IN;
}



