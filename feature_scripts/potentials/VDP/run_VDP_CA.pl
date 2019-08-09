#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

# input
my $f_pdb = $ARGV[0];
my $f_ss = $ARGV[1];
my $f_vol = $ARGV[2];
# output
my $f_out = $ARGV[3];

my $dir_work = $ARGV[4];
my $f_u1 = "$dir_work/utils/ave_min_max.pl";
my $f_u2 = "$dir_work/potentials/VDP/parse_stride.pl";
require "$f_u1";
require "$f_u2";
my $f_refer = "$dir_work/potentials/VDP/VDP_CA_refer_10_30_2.txt";

#############################################
my $range1 = 10;
my $range2 = 30;
my $bin = 2;
my $num_class = ($range2-$range1)/$bin;
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
if(! -e $f_ss || -z $f_ss){ die "Error, calculating VDP, no stride file."; }
parse_stride($f_ss, \%ss);

my $N = keys %ss;
if($N == 0){ 
  close OUT; 
  next; 
}

if(! -e $f_vol || -z $f_vol){ die "Error, calculating VDP, no volume file."; }

open IN, "$f_vol" or die;
while(my $line = <IN>){
  chomp $line; $line =~ s/^\s+//; $line =~ s/\s+$//;
  my $atom_name = substr($line, 12, 4);
  $atom_name =~ s/^\s+//; $atom_name =~ s/\s+$//;
  my $res = substr($line, 17, 3);
  if($atom_name ne 'CA'){ next; }
  my $res_id = substr($line, 22, 4);
  $res_id =~ s/^\s+//; $res_id =~ s/\s+$//;

  if($res ~~ @aminos){
      my @items = split(/\s+/, $line);
      my $tmp_n = @items;
      my $vol = $items[$tmp_n-2];
      if($vol == -1){ next; }

      my $class = get_class($range1, $range2, $bin, $vol);
      if($class == 0){ next; } 
    
      if(exists $ss{$res_id}){
        my $ss1 = $ss{$res_id};
        my $asp = $penalty;
        my $f_o = 0;
        my $f_r = 0;
        my $ma = $Ma{$res};
        
        my $add = 0;
        if($ss1 eq 'C'){ $add = 0; }
        if($ss1 eq 'E'){ $add = 1; }
        if($ss1 eq 'H'){ $add = 2; }
        if($ma == 0){ next; }
        $f_o = $refers[$aminoID{$res}][$add*$num_class+$class-1] / $ma;

        my $fr1 = 0;
        for(my $i=0; $i<20; $i++){
          $fr1 += $refers[$i][$add*$num_class+$class-1];
        }
        $f_r = $fr1 / $Ma_sum; 
        if($f_r){ 
            $asp = ($RT*log(1+$ma*$sigma)) - ($RT*log(1+($ma*$sigma*($f_o/$f_r))));
        }
        print OUT "$res_id\t$res\t$asp\n";
      }
  }
}
close OUT;


sub get_class{
  my ($range1, $range2, $bin, $d) = @_;
  my $n = 0;
  for(my $i=$range1; $i<$range2; $i+=$bin){
    my $j = $i+$bin;
    $n++;
    if($j >= $range2){
      if($d>= $i && $d<=$j){ return $n;}
    }else{
      if($d>= $i && $d<$j){ return $n;}
    }
  }
  return 0;
}



