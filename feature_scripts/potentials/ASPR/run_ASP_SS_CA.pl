#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

my $f_pdb = $ARGV[0];
my $f_ss = $ARGV[1];
my $f_out = $ARGV[2];

my $dir_work = $ARGV[3];
my $f_u1 = "$dir_work/utils/ave_min_max.pl";
my $f_u2 = "$dir_work/potentials/ASPR/parse_stride.pl";
require "$f_u1";
require "$f_u2";
my $f_refer = "$dir_work/potentials/ASPR/SS_ASP_refer_CA.11.50.2.25";

#####################
# configure values
my @aminos = ("ASP","PRO","LYS","ILE","TRP","CYS","GLY","PHE","GLN","SER","ASN","LEU","VAL","TYR","GLU","ARG","THR","ALA","MET","HIS");
my $aid = 1;
my $RT = 0.582;
my $sigma = 0.02;
my $range_distance = 11;
my $range_burial = 50;
my $res = 2;
my $penalty = 0; #
my $num_classes = $range_burial / $res;

my %refers;
my %refer_sum;
my $Ma_sum = 0;
open IN, "$f_refer" or die $!;
while(my $line = <IN>){
  chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;
  my @items = split(/\s+/, $line);
  my $amino = $items[0];
  shift @items;
  $refer_sum{$amino} = ave_min_max(\@items, 4);
  $Ma_sum += $refer_sum{$amino}; 
  $refers{$amino} = \@items;
}close IN;

open OUT, ">$f_out";
############
# ss
my %ss;
if(! -e $f_ss || -z $f_ss){ die "Error, calculating ASPR, no stride file."; }
parse_stride($f_ss, \%ss); 

my %coords;
get_coords($f_pdb, $aid, \%coords);
my $N = keys %coords;
if($N == 0){
  close OUT;
  exit;
}


my $asp_sum = 0;
for my $k1 (sort {$a<=>$b} keys %coords){
    my $num_atoms = 0;
    for my $k2 (keys %coords){
        if($k2 eq $k1){ next; }
        my $dist = sqrt(($coords{$k1}->[0]-$coords{$k2}->[0])**2 + ($coords{$k1}->[1]-$coords{$k2}->[1])**2 + ($coords{$k1}->[2]-$coords{$k2}->[2])**2);
        if($dist <= $range_distance){
          $num_atoms++;
        }
     }
     my $amino = $coords{$k1}->[3];
     if($amino eq "GLY" && $aid == 2){ next; }

     my $class = get_class($range_burial, $res, $num_atoms);
     my $asp = $penalty;
     if($class == 0){  
        $asp_sum += $penalty;
     }else{

         my $amino_ss;
         if(exists $ss{$coords{$k1}->[4]}){
            $amino_ss = $ss{$coords{$k1}->[4]};
         }else{
            next;
         }
        my $add = 0;
        if($amino_ss eq 'C'){ $add = 0; }
        if($amino_ss eq 'E'){ $add = 1; }
        if($amino_ss eq 'H'){ $add = 2; }
        my $key_refer = $add*$num_classes + $class - 1;

        my $ma = $refer_sum{$amino};
        if($ma == 0){ next; }

        my $f_observed_a = $refers{$amino}->[$key_refer]/$ma;
        my $f_reference;
        my $ma_sum = 0;
        for my $k3 (keys %refers){
          $ma_sum += $refers{$k3}->[$key_refer];
        }
        $f_reference = $ma_sum/$Ma_sum;
        
        #if($f_reference){
        if($f_observed_a && $f_reference){
          $asp = ($RT*log(1+$ma*$sigma)) - ($RT*log(1+($ma*$sigma*($f_observed_a/$f_reference))));
          $asp_sum += $asp;
        }else{
          $asp_sum += $penalty;
        }
     }
     print OUT "$k1\t$amino\t$asp\n";
}
close OUT;


sub get_class{
  my ($range, $res, $num) = @_;
  my $n=0;
  for(my $i=0; $i<$range; $i+=$res){
    my $j=$i+$res;
    $n++;
    if($j >= $range){
      if($num >= $i && $num <= $range){ return $n; }
    }else{
      if($num >= $i && $num < $j){ return $n; }
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
  my $cid = 0;
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
        if($residue ~~ @aminos){
          $cid++;
          push(@{$coords->{$cid}}, $x);
          push(@{$coords->{$cid}}, $y);
          push(@{$coords->{$cid}}, $z);
          push(@{$coords->{$cid}}, $residue);
          push(@{$coords->{$cid}}, $id2);
        }
      }
    }
  }
  close IN;
}

