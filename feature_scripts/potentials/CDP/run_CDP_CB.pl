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
my $f_u2 = "$dir_work/potentials/CDP/parse_stride.pl";
require "$f_u1";
require "$f_u2";
my $f_refer = "$dir_work/potentials/CDP/CDP_refer_CB.6.9";

#############################
# global configuration values
my $RT = 0.582;
my $sigma = 0.02;
my $aid = 2;
my $separation = 6;
my $distance = 9;
my $penalty = 0;
my @aminos = ("ASP","PRO","LYS","ILE","TRP","CYS","GLY","PHE","GLN","SER","ASN","LEU","VAL","TYR","GLU","ARG","THR","ALA","MET","HIS");
my %aminoID;
for(my $i=0; $i<20; $i++){
  $aminoID{$aminos[$i]} = $i;
}
#############################
# reference state information
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
if(! -e $f_ss || -z $f_ss){ die "Error, calculating CDP, no stride file."; }
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
    my $amino1 = $coords{$k1}->[3];
    my $rid1 = $coords{$k1}->[4];
    for my $k2 (sort {$a<=>$b} keys %coords){
        #if($k2 < $k1){ next; } 
        my $k12 = abs($k1 - $k2); 
        if($k12 < $separation){ next; }
        my $dist12 = sqrt(($coords{$k1}->[0]-$coords{$k2}->[0])**2 + ($coords{$k1}->[1]-$coords{$k2}->[1])**2 + ($coords{$k1}->[2]-$coords{$k2}->[2])**2);
        if($dist12 >= $distance){ next; }
        my $amino2 = $coords{$k2}->[3];
        my $rid2 = $coords{$k2}->[4];

        my $asp = $penalty;
        my $f_o = 0;
        my $f_r = 0;
        my $ma = $Ma{$amino1};
        
        if($ma == 0){ next; }

        my $amino2_ss;
        if(exists $ss{$coords{$k2}->[4]}){
          $amino2_ss = $ss{$coords{$k2}->[4]};
        }else{
          next;
        }
        my $add = 0;
        if($amino2_ss eq 'C'){ $add = 0; }
        if($amino2_ss eq 'E'){ $add = 1; }
        if($amino2_ss eq 'H'){ $add = 2; }
        
        $f_o = $refers[$aminoID{$amino1}][$add*20+$aminoID{$amino2}] / $ma;

        my $fr1 = 0;
        for(my $i=0; $i<20; $i++){
          $fr1 += $refers[$i][$add*20+$aminoID{$amino2}];
        }
        $f_r = $fr1 / $Ma_sum; 
        if($f_r){ 
            $asp = ($RT*log(1+$ma*$sigma)) - ($RT*log(1+($ma*$sigma*($f_o/$f_r))));
        }
        $asp_sum += $asp;
        print OUT "$rid1\t$amino1\t$rid2\t$amino2\t$asp\n";
      }
}
close OUT;


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

