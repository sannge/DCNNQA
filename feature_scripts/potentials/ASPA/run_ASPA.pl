#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

# input
my $f_pdb = $ARGV[0];
# output 
my $f_out = $ARGV[1];

########################################
# global configuration values
my $RT = 0.582;
my $sigma = 0.02;
my $range_distance = 8;
my $range1 = 50;
my $range2 = 200;
my $res = 5;
my $penalty = 0;

my ($refer, $atomic_type) = get_refers_atomic_types();
my %refer_sum;
my $Ma_sum = 0;
for my $key (keys %{$refer}){
  my $tmp_sum = 0;
  for(my $i=0; $i<@{$refer->{$key}}; $i++){
    $tmp_sum += $refer->{$key}->[$i];
  }
  $refer_sum{$key} = $tmp_sum;
  $Ma_sum += $tmp_sum;
}

open OUT, ">$f_out";
my %coords;
get_coords($f_pdb, \%coords);
my $N = keys %coords;
if($N == 0){
  close OUT;
  exit;
}


my $score_sum = 0;
for my $k1 (sort {$a<=>$b} keys %coords){
      my $num_atoms = 0;
      for my $k2 (keys %coords){
        if($k2 eq $k1){ next; }
        my $dist = sqrt(($coords{$k1}->[0]-$coords{$k2}->[0])**2 + ($coords{$k1}->[1]-$coords{$k2}->[1])**2 + ($coords{$k1}->[2]-$coords{$k2}->[2])**2);
        if($dist <= $range_distance){
          $num_atoms++;
        }
      }
      my $amino = $coords{$k1}->[5];
      my $class = get_class($range1, $range2, $res, $num_atoms);
      my $score = $penalty;
      if($class == 0){  
        $score_sum += $penalty;
      }else{
        my $ma = $refer_sum{$amino};
        my $f_observed_a = $refer->{$amino}->[$class-1] / $ma;
        my $f_reference;
        my $ma_sum = 0;
        for my $k3 (keys %{$refer}){
          $ma_sum += $refer->{$k3}->[$class-1];
        }
        $f_reference = $ma_sum / $Ma_sum;
        
        if($f_observed_a && $f_reference){
          $score = ($RT*log(1+$ma*$sigma)) - ($RT*log(1+($ma*$sigma*($f_observed_a/$f_reference))));
          $score_sum += $score;
        }else{
          $score_sum += $penalty;
        }
      }
      print OUT "$k1\t$amino\t$score\n";
}
close OUT;

sub get_class{
  my ($range1, $range2, $res, $num) = @_;
  my $n=0;
  for(my $i=$range1; $i<$range2; $i+=$res){
    my $j=$i+$res;
    $n++;
    if($j >= $range2){
      if($num >= $i && $num <= $j){ return $n; }
    }else{
      if($num >= $i && $num < $j){ return $n; }
    }
  }
  return 0;
}

sub get_coords{
  my ($f_pdb, $coords) = @_;
  open IN, "$f_pdb" or die $!;
  my $rid = 0;
  my $cur_id = -1;
  my $rid2 = 0;
  while(my $line = <IN>){
    chomp $line;
    if(substr($line, 0, 6) eq 'ATOM  '){
      my $atom_name = substr($line,12,4);
      $atom_name =~ s/^\s+//; $atom_name =~ s/\s+$//;
      my $residue = substr($line, 17, 3);
      my $key = "$residue,$atom_name";
      if(exists $atomic_type->{$key}){
        my ($id1, $id2, $x, $y, $z);
        $id1 = substr($line, 6, 5);
        $id2 = substr($line, 22, 4);
        $x = substr($line, 30, 8);
        $y = substr($line, 38, 8);
        $z = substr($line, 46, 8);
        $id1 =~ s/^\s+//; $id1 =~ s/\s+$//;
        $id2 =~ s/^\s+//; $id2 =~ s/\s+$//;
        $x =~ s/^\s+//; $x =~ s/\s+$//;
        $y =~ s/^\s+//; $y =~ s/\s+$//;
        $z =~ s/^\s+//; $z =~ s/\s+$//;
        $rid++;
        if($id2 ne $cur_id){
          $cur_id = $id2;
          $rid2++;
        }
        push(@{$coords->{$rid}}, $x);
        push(@{$coords->{$rid}}, $y);
        push(@{$coords->{$rid}}, $z);
        push(@{$coords->{$rid}}, $residue);
        push(@{$coords->{$rid}}, $rid2);
        push(@{$coords->{$rid}}, $atomic_type->{$key});
      }
    }
  }
  close IN;
}

sub get_refers_atomic_types{
  my (%refers, %atomic_types);
  my @arr1 = qw(55995 79731 100410 121110 135403 141085 139104 137645 138877 143395 146698 136983 99917 54120 21612 6929 2223 808 364 206 149 101 74 54 49 42 34 45 38 41);
  my @arr2 = qw(9571 10432 10623 10026 9250 8790 8218 8315 8444 8669 9103 8813 7165 4313 1950 617 192 49 27 10 7 2 3 2 4 6 1 2 2 3);
  my @arr3 = qw(45430 67283 93022 118803 139375 151665 156060 156397 156388 156804 154790 136127 97901 54289 22420 7630 2460 875 362 212 146 75 90 69 47 52 44 44 46 42);
  my @arr4 = qw(54106 77116 102128 128362 148962 159736 162481 160706 158547 157073 153608 136816 100074 56211 24244 8152 2635 900 399 214 142 105 81 69 61 51 42 35 34 50);
  my @arr5 = qw(66689 82103 98605 115907 131194 142126 150575 155955 158076 157073 154541 136207 98862 56945 24854 8839 2710 894 443 218 133 99 73 68 61 50 46 30 35 57);
  my @arr6 = qw(31154 35423 38858 40873 43834 47623 52456 58995 69864 86791 111300 132305 123656 83489 39315 13562 3909 1317 509 215 134 86 63 35 28 25 28 23 24 32);
  my @arr7 = qw(6832 9204 11943 14308 16841 19884 23260 27328 31964 38878 48179 54848 48658 30708 12987 4342 1230 401 213 100 58 24 18 15 16 7 3 9 17 12);
  my @arr8 = qw(103613 117126 126144 130464 126646 120747 114721 109196 105669 102737 104975 105301 91576 59587 28431 10027 3135 981 439 234 157 79 63 45 40 44 36 23 23 36);
  my @arr9 = qw(887 932 1070 1150 1239 1428 1624 1759 2001 2365 3031 3676 3607 2876 1488 536 174 54 23 7 9 1 3 0 3 3 0 0 1 1);
  my @arr10 = qw(4306 5920 7373 8459 8433 8129 7318 6437 5532 4846 4201 3618 2823 1770 870 334 83 25 16 2 7 3 2 0 4 2 2 3 3 1);
  my @arr11 = qw(2584 3465 4285 5266 6399 7802 9817 11956 13824 15710 17405 19455 18977 14582 8224 3172 1032 306 120 42 28 13 9 7 7 6 7 9 2 7);
  my @arr12 = qw(16668 20000 23703 27896 32554 37474 42692 48341 55351 64362 75393 83344 78654 55798 29059 11016 3364 1024 334 156 87 46 45 21 27 36 29 26 14 24);
  my @arr13 = qw(373 482 637 852 1046 1327 1618 2014 2242 2573 2873 2988 2898 2316 1367 547 199 61 12 10 5 3 1 1 2 1 0 0 2 2);
  my @arr14 = qw(528 664 912 1117 1214 1538 1797 1981 2196 2286 2472 2637 2626 1934 1062 431 136 45 12 7 4 1 0 3 1 3 1 0 0 0);
  my @arr15 = qw(7547 8089 8329 8046 7348 6689 6439 6236 6074 6098 6161 6206 5787 4340 2420 902 317 99 36 27 14 9 6 3 1 0 0 0 6 2);
  my @arr16 = qw(13976 16070 16055 14871 14025 13352 12763 11984 11556 11281 11568 11654 10949 8360 4963 1960 634 200 69 26 22 13 8 10 5 2 3 7 2 3);
  my @arr17 = qw(5721 6795 7293 7483 7210 6710 6828 6549 6296 6112 6176 6417 5764 4078 2154 749 225 77 19 16 4 4 3 4 5 0 4 3 0 1);
  my @arr18 = qw(12156 11220 10083 9190 8307 7542 6571 5760 5388 4945 4709 4619 4045 3015 1550 563 173 56 33 7 10 4 1 5 2 4 3 0 2 3);
  my @arr19 = qw(361 364 456 583 723 903 1064 1360 1682 2295 2791 3211 3250 2450 1438 619 246 92 39 16 3 2 1 2 2 1 2 2 3 0);
  my @arr20 = qw(8119 6922 5442 4203 3235 2500 1866 1564 1191 995 874 681 484 312 131 39 14 7 4 1 1 1 0 0 0 0 1 1 0 0);
  my @arr21 = qw(8351 8633 8083 7072 6234 5388 4714 3826 3174 2619 2130 1710 1337 843 529 233 50 23 9 7 6 4 1 2 0 0 0 0 1 1);
  my @arr22 = qw(16120 15842 14383 12896 11151 9775 8420 7062 5974 4732 4042 3300 2485 1607 866 338 94 27 22 15 7 6 2 2 0 1 0 2 4 0);
  my @arr23 = qw(1989 2336 2630 2909 3014 3256 3433 3212 3150 2786 2649 2663 2270 1772 1071 469 173 40 13 11 2 4 1 2 3 0 2 0 0 3);
  my @arr24 = qw(2756 3118 3592 3936 4359 4605 4833 4979 4939 5040 5114 5031 4791 3618 2150 830 271 103 23 10 7 5 7 2 4 1 2 1 0 2);
  my @arr25 = qw(2508 2700 2810 2922 2814 2775 2633 2354 2309 2204 2132 1987 1783 1383 669 274 104 26 8 0 1 2 2 4 1 0 0 1 1 4);
  my @arr26 = qw(2705 2976 2939 2984 2886 2587 2548 2306 2126 2087 2053 1934 1807 1228 720 276 96 20 9 3 2 2 3 0 0 1 3 0 2 0);
  my @arr27 = qw(23587 22948 20452 17733 14814 12471 10296 8436 7204 6253 5634 4987 4059 2581 1353 529 176 57 32 10 9 8 5 2 2 3 0 3 2 1);
  my @arr28 = qw(44549 40950 35961 30334 26043 21606 17814 15155 13243 11814 10511 9368 7413 5063 2689 1070 353 125 46 27 13 14 6 9 5 4 1 2 3 3);
  my @arr29 = qw(1110 1402 1769 1934 2279 2545 3005 3376 4105 4833 5802 6951 6614 5033 2683 1071 368 133 47 22 8 6 7 8 2 2 1 0 3 3);
  my @arr30 = qw(1013 1078 1160 1208 1399 1480 1497 1665 1966 2232 2941 3664 3558 2420 1176 428 131 49 21 5 2 2 4 4 2 0 0 0 0 0);
  my @arr31 = qw(2113 2615 3296 3859 4258 4801 4968 4939 4815 4742 4869 4826 4390 3149 1695 597 176 43 19 13 5 3 3 3 3 1 0 2 0 1);
  my @arr32 = qw(5978 6755 6670 6817 6384 6067 5697 5133 4872 4443 4210 3865 3119 2149 1026 371 93 35 17 4 6 0 6 4 5 1 0 1 1 1);
  my @arr33 = qw(12616 12706 12036 10970 9858 8852 7821 7057 6299 5584 5149 4829 4478 3263 1805 684 188 64 31 14 7 5 7 3 1 4 0 3 0 2);
  my @arr34 = qw(12310 12062 11280 10340 9296 8393 7463 6763 6170 5586 5142 4777 4180 3262 1913 878 255 92 19 17 8 5 5 5 2 1 1 2 2 3);
  my @arr35 = qw(9523 8555 7363 6079 4783 3757 2970 2130 1649 1386 1097 892 688 362 192 88 19 9 3 1 2 0 0 0 0 0 1 1 1 1);
  my @arr36 = qw(7898 8444 8203 7906 6990 6224 5423 4510 3976 3284 2641 2189 1580 1129 638 279 85 36 14 11 5 5 3 2 3 0 0 0 0 1);
  my @arr37 = qw(7585 7856 8059 7898 7349 6763 6236 5449 4710 4142 3455 2904 2275 1484 815 409 118 48 14 10 5 5 2 6 1 6 3 0 0 0);
  my @arr38 = qw(2379 2735 2814 3113 3087 3147 3006 2906 2574 2361 2371 2319 2067 1622 1009 480 152 42 17 7 3 4 1 1 3 1 1 2 1 0);
  my @arr39 = qw(752 914 1074 1237 1493 1721 1815 1984 1941 2147 2394 2509 2374 1737 919 328 124 29 11 7 2 2 2 1 3 0 0 0 2 0);
  my @arr40 = qw(3016 3590 3779 4215 4173 4112 3923 4010 3949 4050 4070 4216 3809 2818 1426 495 164 31 16 7 6 3 3 2 1 3 2 0 0 0);
  $refers{"1"} = \@arr1; $refers{"2"} = \@arr2; $refers{"3"} = \@arr3; $refers{"4"} = \@arr4; $refers{"5"} = \@arr5; 
  $refers{"6"} = \@arr6; $refers{"7"} = \@arr7; $refers{"8"} = \@arr8; $refers{"9"} = \@arr9; $refers{"10"} = \@arr10; 
  $refers{"11"} = \@arr11; $refers{"12"} = \@arr12; $refers{"13"} = \@arr13; $refers{"14"} = \@arr14; $refers{"15"} = \@arr15; 
  $refers{"16"} = \@arr16; $refers{"17"} = \@arr17; $refers{"18"} = \@arr18; $refers{"19"} = \@arr19; $refers{"20"} = \@arr20; 
  $refers{"21"} = \@arr21; $refers{"22"} = \@arr22; $refers{"23"} = \@arr23; $refers{"24"} = \@arr24; $refers{"25"} = \@arr25; 
  $refers{"26"} = \@arr26; $refers{"27"} = \@arr27; $refers{"28"} = \@arr28; $refers{"29"} = \@arr29; $refers{"30"} = \@arr30; 
  $refers{"31"} = \@arr31; $refers{"32"} = \@arr32; $refers{"33"} = \@arr33; $refers{"34"} = \@arr34; $refers{"35"} = \@arr35; 
  $refers{"36"} = \@arr36; $refers{"37"} = \@arr37; $refers{"38"} = \@arr38; $refers{"39"} = \@arr39; $refers{"40"} = \@arr40; 
  $atomic_types{"GLY,N"} = 3; $atomic_types{"GLY,CA"} = 2; $atomic_types{"GLY,C"} = 4; $atomic_types{"GLY,O"} = 5; $atomic_types{"ALA,N"} = 3; $atomic_types{"ALA,CA"} = 1; 
  $atomic_types{"ALA,C"} = 4; $atomic_types{"ALA,O"} = 5; $atomic_types{"ALA,CB"} = 6; $atomic_types{"VAL,N"} = 3; $atomic_types{"VAL,CA"} = 1; $atomic_types{"VAL,C"} = 4; 
  $atomic_types{"VAL,O"} = 5; $atomic_types{"VAL,CB"} = 7; $atomic_types{"VAL,CG1"} = 6; $atomic_types{"VAL,CG2"} = 6; $atomic_types{"LEU,N"} = 3; $atomic_types{"LEU,CA"} = 1; 
  $atomic_types{"LEU,C"} = 4; $atomic_types{"LEU,O"} = 5; $atomic_types{"LEU,CB"} = 8; $atomic_types{"LEU,CG"} = 7; $atomic_types{"LEU,CD1"} = 6; $atomic_types{"LEU,CD2"} = 6; 
  $atomic_types{"ILE,N"} = 3; $atomic_types{"ILE,CA"} = 1; $atomic_types{"ILE,C"} = 4; $atomic_types{"ILE,O"} = 5; $atomic_types{"ILE,CB"} = 7; $atomic_types{"ILE,CG1"} = 8; 
  $atomic_types{"ILE,CG2"} = 6; $atomic_types{"ILE,CD1"} = 6; $atomic_types{"MET,N"} = 3; $atomic_types{"MET,CA"} = 1; $atomic_types{"MET,C"} = 4; $atomic_types{"MET,O"} = 5; 
  $atomic_types{"MET,CB"} = 8; $atomic_types{"MET,CG"} = 29; $atomic_types{"MET,SD"} = 9; $atomic_types{"MET,CE"} = 30; $atomic_types{"PRO,N"} = 10; $atomic_types{"PRO,CA"} = 1; 
  $atomic_types{"PRO,C"} = 4; $atomic_types{"PRO,O"} = 5; $atomic_types{"PRO,CB"} = 8; $atomic_types{"PRO,CG"} = 8; $atomic_types{"PRO,CD"} = 32; $atomic_types{"PHE,N"} = 3; 
  $atomic_types{"PHE,CA"} = 1; $atomic_types{"PHE,C"} = 4; $atomic_types{"PHE,O"} = 5; $atomic_types{"PHE,CB"} = 8; $atomic_types{"PHE,CG"} = 11; $atomic_types{"PHE,CD1"} = 12; 
  $atomic_types{"PHE,CD2"} = 12; $atomic_types{"PHE,CE1"} = 12; $atomic_types{"PHE,CE2"} = 12; $atomic_types{"PHE,CZ"} = 12; $atomic_types{"TRP,N"} = 3; $atomic_types{"TRP,CA"} = 1; 
  $atomic_types{"TRP,C"} = 4; $atomic_types{"TRP,O"} = 5; $atomic_types{"TRP,CB"} = 8; $atomic_types{"TRP,CG"} = 13; $atomic_types{"TRP,CD1"} = 24; $atomic_types{"TRP,CD2"} = 11; 
  $atomic_types{"TRP,NE1"} = 39; $atomic_types{"TRP,CE2"} = 14; $atomic_types{"TRP,CE3"} = 12; $atomic_types{"TRP,CZ2"} = 12; $atomic_types{"TRP,CZ3"} = 12; $atomic_types{"TRP,CH2"} = 12; 
  $atomic_types{"SER,N"} = 3; $atomic_types{"SER,CA"} = 1; $atomic_types{"SER,C"} = 4; $atomic_types{"SER,O"} = 5; $atomic_types{"SER,CB"} = 15; $atomic_types{"SER,OG"} = 16; 
  $atomic_types{"THR,N"} = 3; $atomic_types{"THR,CA"} = 1; $atomic_types{"THR,C"} = 4; $atomic_types{"THR,O"} = 5; $atomic_types{"THR,CB"} = 17; $atomic_types{"THR,OG1"} = 16; 
  $atomic_types{"THR,CG2"} = 6; $atomic_types{"ASN,N"} = 3; $atomic_types{"ASN,CA"} = 1; $atomic_types{"ASN,C"} = 4; $atomic_types{"ASN,O"} = 5; $atomic_types{"ASN,CB"} = 8; 
  $atomic_types{"ASN,CG"} = 33; $atomic_types{"ASN,OD1"} = 34; $atomic_types{"ASN,ND2"} = 18; $atomic_types{"GLN,N"} = 3; $atomic_types{"GLN,CA"} = 1; $atomic_types{"GLN,C"} = 4; 
  $atomic_types{"GLN,O"} = 5; $atomic_types{"GLN,CB"} = 8; $atomic_types{"GLN,CG"} = 8; $atomic_types{"GLN,CD"} = 33; $atomic_types{"GLN,OE1"} = 34; $atomic_types{"GLN,NE2"} = 18; 
  $atomic_types{"LYS,N"} = 3; $atomic_types{"LYS,CA"} = 1; $atomic_types{"LYS,C"} = 4; $atomic_types{"LYS,O"} = 5; $atomic_types{"LYS,CB"} = 8; $atomic_types{"LYS,CG"} = 8; 
  $atomic_types{"LYS,CD"} = 8; $atomic_types{"LYS,CE"} = 35; $atomic_types{"LYS,NZ"} = 20; $atomic_types{"TYR,N"} = 3; $atomic_types{"TYR,CA"} = 1; $atomic_types{"TYR,C"} = 4; 
  $atomic_types{"TYR,O"} = 5; $atomic_types{"TYR,CB"} = 8; $atomic_types{"TYR,CG"} = 11; $atomic_types{"TYR,CD1"} = 12; $atomic_types{"TYR,CD2"} = 12; $atomic_types{"TYR,CE1"} = 12; 
  $atomic_types{"TYR,CE2"} = 12; $atomic_types{"TYR,CZ"} = 31; $atomic_types{"TYR,OH"} = 40; $atomic_types{"CYS,N"} = 3; $atomic_types{"CYS,CA"} = 1; $atomic_types{"CYS,C"} = 4; 
  $atomic_types{"CYS,O"} = 5; $atomic_types{"CYS,CB"} = 29; $atomic_types{"CYS,SG"} = 19; $atomic_types{"GLU,N"} = 3; $atomic_types{"GLU,CA"} = 1; $atomic_types{"GLU,C"} = 4; 
  $atomic_types{"GLU,O"} = 5; $atomic_types{"GLU,CB"} = 8; $atomic_types{"GLU,CG"} = 8; $atomic_types{"GLU,CD"} = 27; $atomic_types{"GLU,OE1"} = 28; $atomic_types{"GLU,OE2"} = 28; 
  $atomic_types{"ASP,N"} = 3; $atomic_types{"ASP,CA"} = 1; $atomic_types{"ASP,C"} = 4; $atomic_types{"ASP,O"} = 5; $atomic_types{"ASP,CB"} = 8; $atomic_types{"ASP,CG"} = 27; 
  $atomic_types{"ASP,OD1"} = 28; $atomic_types{"ASP,OD2"} = 28; $atomic_types{"ARG,N"} = 3; $atomic_types{"ARG,CA"} = 1; $atomic_types{"ARG,C"} = 4; $atomic_types{"ARG,O"} = 5; 
  $atomic_types{"ARG,CB"} = 8; $atomic_types{"ARG,CG"} = 8; $atomic_types{"ARG,CD"} = 37; $atomic_types{"ARG,NE"} = 36; $atomic_types{"ARG,CZ"} = 21; $atomic_types{"ARG,NH1"} = 22; 
  $atomic_types{"ARG,NH2"} = 22; $atomic_types{"HIS,N"} = 3; $atomic_types{"HIS,CA"} = 1; $atomic_types{"HIS,C"} = 4; $atomic_types{"HIS,O"} = 5; $atomic_types{"HIS,CB"} = 8; 
  $atomic_types{"HIS,CG"} = 23; $atomic_types{"HIS,ND1"} = 38; $atomic_types{"HIS,CD2"} = 24; $atomic_types{"HIS,CE1"} = 26; $atomic_types{"HIS,NE2"} = 25;

  return (\%refers, \%atomic_types);
} 
