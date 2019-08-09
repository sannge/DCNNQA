#!/usr/bin/env perl
use warnings;
use strict;
use POSIX;

sub ave_min_max{
  my ($arr, $index) = @_;
  my $n = @{$arr};
  if($n <= 0){ print "Error, when calculating ave_min_max...\n"; exit; }
  my $val;
  if($index == 1){ # get the average value of the give array
    my $sum = 0;
    for(my $i=0; $i<$n; $i++){
      $sum += $arr->[$i];
    }
    $val = $sum/$n;
  }
  if($index == 2){ # get the min value of the give array
    my $min = $arr->[0];
    for(my $i=1; $i<$n; $i++){
      if($arr->[$i] < $min){
        $min = $arr->[$i];
      }
    }
    $val = $min;
  }
  if($index == 3){ # get the max value of the give array
    my $max = $arr->[0];
    for(my $i=1; $i<$n; $i++){
      if($arr->[$i] > $max){
        $max = $arr->[$i];
      }
    }
    $val = $max;
  }
  if($index == 4){ # get the sum value of the give array
    my $sum = 0;
    for(my $i=0; $i<$n; $i++){
      $sum += $arr->[$i];
    }
    $val = $sum;
  }

  return $val;
}

1;

