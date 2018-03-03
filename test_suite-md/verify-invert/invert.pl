#! /usr/bin/perl
#
# Invert all coordinates in an xyz file by a given amount.
#
use strict ;
use warnings ;

my $natoms = <> ;
print $natoms ;

$_ = <> ;
print $_ ;

for my $i ( 1 .. $natoms ) {
  $_ = <> ;
  my ( $ele, $x, $y, $z) = split(" ") ;
  $x = -$x ;
  $y = -$y ;
  $z = -$z ;
  printf("%s %21.14e %21.14e %21.14e\n", $ele, $x, $y, $z) ;
}

